#include <stdio.h>
#include <math.h>

#define PI 3.14159265359

struct point{
	double x;
	double y;
	double z;
};

struct res{
	double t1;
	double t2;
};

struct color{
	int r;
	int g;
	int b;
};

struct sphere{
	struct point center;
	double r;
	struct color color; //nesting structs?
};

struct light{
	int type; //0 = ambient, 1 = point, 2 = directional 
	double i;
	struct point center;
	struct point direction;
};

struct scene{
	struct sphere spheres[30]; 
	struct light lights[4];
	int nspheres;
	int nlights;
};

struct scene myScene; //need to figure out how to eliminate this dependancy!

struct color bgColor; //not needed
struct point camera; //not needed

struct color canvas[1920][1080]; //need to put this here because memory overflow
#pragma acc declare link(canvas) //only good for static vars

int screenWidth; //both not needed 
int screenHeight;

double viewWidth; //all 3 not needed 
double viewHeight;
double viewDist;

void loadScene(char* fileName){
	FILE *p;
    p = fopen(fileName, "r");

    fscanf(p, "%d", &screenWidth); fscanf(p, "%d", &screenHeight);

    fscanf(p, "%d", &bgColor.r); fscanf(p, "%d", &bgColor.g); fscanf(p, "%d", &bgColor.b);

    fscanf(p, "%lf", &camera.x); fscanf(p, "%lf", &camera.y); fscanf(p, "%lf", &camera.z);

    fscanf(p, "%lf", &viewWidth); fscanf(p, "%lf", &viewHeight); fscanf(p, "%lf", &viewDist);

    fscanf(p, "%d", &myScene.nspheres);
    for(int i = 0; i < myScene.nspheres; i++){
    	struct point c;
    	fscanf(p, "%lf", &c.x); fscanf(p, "%lf", &c.y); fscanf(p, "%lf", &c.z);
    	struct color sc;
    	fscanf(p, "%d", &sc.r); fscanf(p, "%d", &sc.g); fscanf(p, "%d", &sc.b);
    	struct sphere mySphere;
    	mySphere.center = c; mySphere.color = sc;
    	fscanf(p, "%lf", &mySphere.r); 
    	myScene.spheres[i] = mySphere;
    }

    fscanf(p, "%d", &myScene.nlights);
    for(int i = 0; i < myScene.nlights; i++){
    	struct light myLight;
    	fscanf(p, "%d", &myLight.type);
    	if(myLight.type == 1){
	    	struct point c;
	    	fscanf(p, "%lf", &c.x); fscanf(p, "%lf", &c.y); fscanf(p, "%lf", &c.z);
	    	myLight.center = c;
    	}
    	if(myLight.type == 2){
	    	struct point d;
	    	fscanf(p, "%lf", &d.x); fscanf(p, "%lf", &d.y); fscanf(p, "%lf", &d.z);
	    	myLight.direction = d;
    	}
    	fscanf(p, "%lf", &myLight.i);
    	myScene.lights[i] = myLight;
    }

    fclose(p);
}

//how to get function to return struct?
#pragma acc routine seq
struct point CanvasToViewport(int x, int y){
	struct point i;
	i.x = (float)x * 1.0 / 1920.0; i.y = (float)y * 0.6 / (float)1080.0; i.z = (float)1.0;

	return i;
}

#pragma acc routine seq
float dot(struct point v1, struct point v2){
	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

#pragma acc routine seq
struct res intersectRaySphere(struct point origin, struct point target, struct sphere collide){
	struct res ret;
	ret.t1 = INFINITY; ret.t2 = INFINITY;

	double r = collide.r;
	struct point oc;
	oc.x = origin.x - collide.center.x; oc.y = origin.y - collide.center.y; oc.z = origin.z - collide.center.z;

	float k1 = dot(target, target);
	float k2 = 2 * dot(oc, target);
	float k3 = dot(oc, oc) - r * r;

	float discrim = k2 * k2 - 4 * k1 * k3;
	if(discrim < 0) return ret;

	ret.t1 = (-k2 + sqrt(discrim)) / (2 * k1);
	ret.t2 = (-k2 - sqrt(discrim)) / (2 * k1);

	return ret;
}

#pragma acc routine seq
double length(struct point L){
	return sqrt((L.x * L.x) + (L.y * L.y) + (L.z * L.z));
}

#pragma acc routine seq
double computeLighting(struct point P, struct point N, struct scene myScene1){
	double i = 0.0;
	for(int j = 0; j < 3; j++){
		if(myScene1.lights[j].type == 0){
			i += myScene1.lights[j].i;
		} 
		else{
			struct point L;

			if(myScene1.lights[j].type == 1){
				L = myScene1.lights[j].center;
				L.x -= P.x; L.y -= P.y; L.z -= P.z;
			}
			else{
				L = myScene1.lights[j].direction;
			}

			double nL = dot(N, L);
			if(nL > 0) i += myScene1.lights[j].i * nL / (length(N) * length(L));
		}
	}

	return i;
}

#pragma acc routine seq
struct color traceSphere(struct point origin, struct point target, double t_min, double t_max, int day, struct scene myScene1){
	double closest_t = INFINITY;
	struct sphere closest_sphere; //how to do NULL in C?
	closest_sphere.r = 0;
	
	for(int i = 0; i < 3; i++){
		struct sphere temp = myScene1.spheres[i];

		//begin super specific animation-related code
		if(i == 1){
			double c = ((2.0 * PI) / 365.0) * (double)day;
			double xoff = sin(c) * 1.1;
			double zoff =  3 - (cos(c) * 1.1);

			//printf("%lf %lf\n", xoff, zoff);

			temp.center.x = xoff;
			temp.center.z = zoff;
		}

		if(i == 2){
			double c = ((2.0 * PI) / 365.0) * (double)day;
			double xoff = sin(c) * 1.1;
			double zoff =  3 - (cos(c) * 1.1);

			double c2 = ((2.0 * PI) / 27) * (double)day;
			double xoff2 = xoff + (sin(c2) * 0.02);
			double zoff2 =  zoff - (cos(c2) * 0.02);

			//printf("%lf %lf\n", xoff2, zoff2);

			temp.center.x = xoff2;
			temp.center.z = zoff2;
		}

		//now back to the regular program

		struct res inter = intersectRaySphere(origin, target, temp);

		if((inter.t1 >= t_min && inter.t1 <= t_max) && inter.t1 < closest_t){
			closest_t = inter.t1;
			closest_sphere = temp;
		}

		if((inter.t2 >= t_min && inter.t2 <= t_max) && inter.t2 < closest_t){
			closest_t = inter.t2;
			closest_sphere = temp;
		}
	}

	struct color bg;
	bg.r = 0; bg.g = 0; bg.b = 0;
	if(closest_sphere.r == 0) return bg;

	struct point p; p.x = origin.x + closest_t * target.x; p.y = origin.y + closest_t * target.y; p.z = origin.z + closest_t * target.z; 

	struct point n;
	n.x = p.x - closest_sphere.center.x; n.y = p.y - closest_sphere.center.y; n.z = p.z - closest_sphere.center.z;

	double ln = length(n);
	n.x /= ln; n.y /= ln; n.z /= ln;

	double ci = computeLighting(p, n, myScene1);

	struct color c;
	c.r = closest_sphere.color.r * ci; c.g = closest_sphere.color.g * ci; c.b = closest_sphere.color.b * ci;
	return c;
}

void save(struct color canvas[1920][1080], int id){ //pass in array in C
	char filename[50];
	sprintf(filename, "output/%d.ppm", id);

	FILE *ppm = fopen(filename,"w+");

	fprintf(ppm,"P3\n%u %u\n255\n", 1920, 1080);

	for(int y = 1080 - 1; y >= 0; y--){
		for(int x = 0; x < 1920; x++){
			fprintf(ppm, "%d %d %d ", 
          		(int)(canvas[x][y].r), 
          		(int)(canvas[x][y].g),
          		(int)(canvas[x][y].b)
          	);
		}
		fprintf(ppm, "\n");
	}

	fclose(ppm);
}

int main(int argc, char* argv[]){
	loadScene(argv[1]);

	//spin up canvas 
	

	int day = 0;
	for(day = 0; day <= 0; day++){

#pragma acc parallel loop collapse(2) copy(canvas[:][:]) copyin(day, myScene)
		for(int y = -(1080 / 2); y < 1080 / 2; y++){
			for(int x = -(1920 / 2); x < 1920 / 2; x++){
				struct point origin; origin.x = 0; origin.y =0; origin.z = 0;
				struct point target = CanvasToViewport(x, y);
				struct color res = traceSphere(origin, target, 1, INFINITY, day, myScene);
				//testing stuff
				struct color col; col.r = 255; col.g = 0; col.b = 0;
				
				canvas[x + (1920 / 2)][y + (1080 / 2)] = col;
			}
		}
printf("done\n");
		save(canvas, day);
	}

	return 0;
}



//Raytracer inspired by Gabriel Gambetta's Computer Graphics from Scratch

//todo: make file for loader, see if it works...
#include <stdio.h>
#include <math.h>

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
	struct sphere spheres[4]; 
	struct light lights[4];
	int nspheres;
	int nlights;
};

struct scene myScene;
struct color bgColor;
struct point camera;

struct color canvas[1024][768];

int screenWidth;
int screenHeight;

int viewWidth;
int viewHeight;
int viewDist;

void loadScene(char* fileName){
	FILE *p;
    p = fopen(fileName, "r");

    fscanf(p, "%d", &screenWidth); fscanf(p, "%d", &screenHeight);

    fscanf(p, "%d", &bgColor.r); fscanf(p, "%d", &bgColor.g); fscanf(p, "%d", &bgColor.b);

    fscanf(p, "%lf", &camera.x); fscanf(p, "%lf", &camera.y); fscanf(p, "%lf", &camera.z);

    fscanf(p, "%d", &viewWidth); fscanf(p, "%d", &viewHeight); fscanf(p, "%d", &viewDist);

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

void setup(){
	struct point sphereCenter;
	sphereCenter.x = 0; sphereCenter.y = 0; sphereCenter.z = 4;
	struct color sphereColor;
	sphereColor.r = 255; sphereColor.g = 0; sphereColor.b = 0;
	struct sphere mySphere;
	mySphere.center = sphereCenter; mySphere.color = sphereColor; mySphere.r = 1;
	myScene.spheres[0] = mySphere; myScene.nspheres = 1;

	struct light ambient;
	ambient.type = 0;
	ambient.i = 0.2;

	struct light point;
	point.type = 1;
	struct point lc1;
	lc1.x = 2; lc1.y = 1; lc1.z = 0;
	point.center = lc1;
	point.i = 0.6;

	struct light directional;
	directional.type = 2;
	struct point lc2;
	lc2.x = 1; lc2.y = 4; lc2.z = 4;
	directional.direction = lc2;
	directional.i = 0.2;

	myScene.nlights = 3;
	myScene.lights[0] = ambient; myScene.lights[1] = point; myScene.lights[2] = directional;

	bgColor.r = 0; bgColor.g = 0; bgColor.b = 0;

	camera.x = 0; camera.y = 0; camera.z = 0;
}

//how to get function to return struct?
struct point CanvasToViewport(int x, int y){
	struct point i;
	i.x = (float)x * (float)viewWidth / (float)screenWidth; i.y = (float)y * (float)viewHeight / (float)screenHeight; i.z = (float)viewDist;

	return i;
}

float dot(struct point v1, struct point v2){
	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

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

double length(struct point L){
	return sqrt((L.x * L.x) + (L.y * L.y) + (L.z * L.z));
}

double computeLighting(struct point P, struct point N){
	double i = 0.0;
	for(int j = 0; j < myScene.nlights; j++){
		if(myScene.lights[j].type == 0){
			i += myScene.lights[j].i;
		} 
		else{
			struct point L;

			if(myScene.lights[j].type == 1){
				L = myScene.lights[j].center;
				L.x -= P.x; L.y -= P.y; L.z -= P.z;
			}
			else{
				L = myScene.lights[j].direction;
			}

			double nL = dot(N, L);
			if(nL > 0) i += myScene.lights[j].i * nL / (length(N) * length(L));
		}
	}

	return i;
}

struct color traceSphere(struct point origin, struct point target, double t_min, double t_max){
	double closest_t = INFINITY;
	struct sphere closest_sphere; //how to do NULL in C?
	closest_sphere.r = 0;
	
	for(int i = 0; i < myScene.nspheres; i++){
		struct sphere temp = myScene.spheres[i];

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

	if(closest_sphere.r == 0) return bgColor;

	struct point p; p.x = origin.x + closest_t * target.x; p.y = origin.y + closest_t * target.y; p.z = origin.z + closest_t * target.z; 

	struct point n;
	n.x = p.x - closest_sphere.center.x; n.y = p.y - closest_sphere.center.y; n.z = p.z - closest_sphere.center.z;

	double ln = length(n);
	n.x /= ln; n.y /= ln; n.z /= ln;

	double ci = computeLighting(p, n);

	struct color c;
	c.r = closest_sphere.color.r * ci; c.g = closest_sphere.color.g * ci; c.b = closest_sphere.color.b * ci;
	return c;
}

void save(struct color canvas[screenWidth][screenHeight]){ //pass in array in C
	FILE *ppm = fopen("canvas.ppm","w+");
	fprintf(ppm,"P3\n%u %u\n255\n", screenWidth, screenHeight);

	for(int y = screenHeight - 1; y >= 0; y--){
		for(int x = 0; x < screenWidth; x++){
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
	//setup();

	struct point origin;
	origin.x = 0; origin.y = 0; origin.z = 0;
	for(int y = -(screenHeight / 2); y < screenHeight / 2; y++){
		for(int x = -(screenWidth / 2); x < screenWidth / 2; x++){
			struct point target = CanvasToViewport(x, y);

			canvas[x + (screenWidth / 2)][y + (screenHeight / 2)] = traceSphere(origin, target, 1, INFINITY);
		}
	}

	save(canvas);
	return 0;
}
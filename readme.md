rayBench is a simple raytracer designed to be run with different accelerator frameworks. It only does spheres
and basic lighting for now, I might add more things later. Inspired by Gabriel Gambetta's Computer Graphics from Scratch.

To run, do ./a.out (input filename)

outputs PPM file to canvas.ppm, can be viewed by most image viewers...

The input file is loosly organized like this...
displayX(int) displayY(int) (max 1024x768)

backgroundColor Red(int) backgroundColor Blue(int) backgroundColorGreen(int)

cameraX(double) cameraY(double) cameraZ(double)

viewWidth(double) viewHeight(double) viewDist(double) (usually 1 1 1)

nSpheres(int)

sphereX(double) sphereY(double) sphereZ(double)
sphereColor red(int) sphereColor green(int) sphereColor blue(int)
sphereRadius(double)

...

nLights(int)

lightType(int, 0=ambient, 1=point, 2=directional)
lightDirectionX(double) lightDirectionY(double) lightDirectionZ(double)
lightCenterX(double) lightCenterY(double) lightCenterZ(double)
lightIntensity(double)

...

Animation program to generate 30fps animation of earth revloving around sun with the moon is in animation folder. rayBench_animated4.c is the most recent version that can be compiled with OpenACC. It should output 365 frames, with each frame representing a day.

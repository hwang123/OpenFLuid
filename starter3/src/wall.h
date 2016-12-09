#ifndef WALL_H
#define WALL_H

#include <vector>
#include <vecmath.h>
#include <cstdint>

#include "vertexrecorder.h"

// helper for uniform distribution
float rand_uniform(float low, float hi);


struct GLProgram;
class Wall
{
public:

	Wall(Vector3f point, Vector3f normal);

    // getter method for the plane's point
	Vector3f getPoint() { return _point; };

    // getter method for the plane's normal
	Vector3f getNormal() { return _normal; };

    // void draw(GLProgram& gl, float len);

 protected:
      // member variables
    Vector3f _point;
    Vector3f _normal;
};

#endif

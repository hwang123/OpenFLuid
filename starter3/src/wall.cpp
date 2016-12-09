#include "wall.h"

#include <iostream>
using namespace std;


Wall::Wall(const Vector3f point, const Vector3f normal)
{
    // default constructor for Particle
    _point = point;
    _normal = normal.normalized();
}

// void Wall::draw(GLProgram& gl, float len) {

// 	float x = _point.x();
//     float y = _point.y();
//     float z = _point.z();

//     gl.updateModelMatrix(Matrix4f::translation(_point)); // update uniforms after mode change

//     // apply a rotation
//     // float angle1 = asin(_normal[0]) / (2 * M_PI) * 360.0;
//     // float angle2 = asin(_normal[1]) / (2 * M_PI) * 360.0;

//     // gl.rotate(-angle1, 0, 1, 0);
//     // gl.rotate(-angle2, 1, 0, 0);

//     const Vector3f color(0.1f, 0.9f, 0.4f);
//     VertexRecorder rec;

//     rec.record(Vector3f(x-len,z,y+len), color);
//     rec.record(Vector3f(x-len,z-len,y+len), color);
//     rec.record(Vector3f(x+len,z,y+len), color);
//     rec.record(Vector3f(x+len,z-len,y+len), color);
//     rec.record(Vector3f(x+-len,z,y-len), color);
//     rec.record(Vector3f(x+-len,z-len,y-len), color);
//     rec.record(Vector3f(x+len,z,y-len), color);
//     rec.record(Vector3f(x+len,z-len,y-len), color);

//     glLineWidth(5.0f);
//     drawQuadCustom(len*2, z-1);
//     rec.draw(GL_LINES);

// }
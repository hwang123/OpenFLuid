#ifndef PENDULUMSYSTEM_H
#define PENDULUMSYSTEM_H

#include <vector>

#include "particlesystem.h"

class PendulumSystem : public ParticleSystem
{
public:
    PendulumSystem();

    std::vector<Vector3f> evalF(std::vector<Vector3f> state) override;
    void draw(GLProgram&);
    Vector3f getSpringForce(Vector3f position1, Vector3f position2);

    // inherits 
    // std::vector<Vector3f> m_vVecState;
};

#endif

#ifndef FLUIDSYSTEM_H
#define FLUIDSYSTEM_H

#include <vector>

#include "particlesystem.h"
#include "particle.h"
#include "wall.h"


class FluidSystem : public ParticleSystem
{
    ///ADD MORE FUNCTION AND FIELDS HERE
public:
    FluidSystem();

    // evalF is called by the integrator at least once per time step
    std::vector<Particle> evalF(std::vector<Particle> state) override;

    // draw is called once per frame
    void draw(GLProgram& ctx);

    // inherits
    // std::vector<Vector3f> m_vVecState;

protected:
	std::vector<Wall> _walls;

};

#endif

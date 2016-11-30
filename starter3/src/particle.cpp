#include "particle.h"

#include <iostream>
using namespace std;


Particle::Particle()
{
    // default constructor for Particle
}

Particle::Particle(const Vector3f position)
{
    _position = position;
}

Particle::Particle(const Vector3f position, const Vector3f velocity)
{
    _position = position;
    _velocity = velocity;
}
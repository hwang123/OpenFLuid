#include "particle.h"

#include <iostream>
using namespace std;


Particle::Particle(int id)
{
    // default constructor for Particle
    _id = id;
}

Particle::Particle(int id, const Vector3f position)
{
	_id = id;
    _position = position;
}

Particle::Particle(int id, const Vector3f position, const Vector3f velocity)
{
	_id = id;
    _position = position;
    _velocity = velocity;
    _density = 998.29;
}
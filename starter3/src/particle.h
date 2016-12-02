#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <vecmath.h>
#include <cstdint>


// helper for uniform distribution
float rand_uniform(float low, float hi);

struct GLProgram;
class Particle
{
public:

	Particle(int id);
	Particle(int id, const Vector3f position);
	Particle(int id, const Vector3f position, const Vector3f velocity);

    // getter method for the particle's position
	Vector3f getPosition() { return _position; };

    // getter method for the particle's velocity
	Vector3f getVelocity() { return _velocity; };

    // getter method for the particle's neighbors
    std::vector<int> getNeighbors() { return _neighbors; };

    // setter method for the system's position
    void setPosition(const Vector3f  newPosition) { _position = newPosition; };

    // setter method for the system's position
    void setVelocity(const Vector3f  newVelocity) { _velocity = newVelocity; };
    
    // setter method for the system's position
    void setNeighbors(const std::vector<int> newNeighbors ) { _neighbors = newNeighbors; };


 protected:
      // member variables
    int _id;
    Vector3f _velocity;
    Vector3f _position;
    std::vector<int> _neighbors; 
};

#endif

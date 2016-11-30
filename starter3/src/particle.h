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

	Particle();
	Particle(const Vector3f position);
	Particle(const Vector3f position, const Vector3f velocity);

    // getter method for the particle's position
	Vector3f getPositon() { return position; };

    // getter method for the particle's velocity
	Vector3f getVelocity() { return velocity; };

    // getter method for the particle's neighbors
    std::vector<int> getNeighbors() { return neighbors; };

    // setter method for the system's position
    void setPosition(const Vector3f  newPosition) { position = newPosition; };

    // setter method for the system's position
    void setVelocity(const Vector3f  newVelocity) { velocity = newVelocity; };
    
    // setter method for the system's position
    void setNieghbors(const std::vector<int> newNeighbors ) { neighbors = newNeighbors; };


 protected:
    std::vector<Vector3f> m_vVecState;
};

private:
    // member variables
    int id;
    Vector3f velocity;
    Vector3f position;
    std::vector<int> neighbors; 
};
#endif

#include "pendulumsystem.h"

#include <cassert>
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
using namespace std;

// TODO adjust to number of particles.
const int NUM_PARTICLES = 4;
const float SPRING_CONSTANT = 50; // N/m
const float PARTICLE_MASS = 1.0; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;
const float SPRING_REST_LENGTH = .1; // m

PendulumSystem::PendulumSystem()
{

    // TODO 4.2 Add particles for simple pendulum
    // TODO 4.3 Extend to multiple particles

    // To add a bit of randomness, use e.g.
    // float f = rand_uniform(-0.5f, 0.5f);
    // in your initial conditions.
    
    m_vVecState.clear();
    // Add the fixed anchor point
    m_vVecState.push_back(Vector3f(0, 0, 0));
    m_vVecState.push_back(Vector3f(0, 0, 0));
    m_vVecState.push_back(Vector3f(0, -.01, 0));
    m_vVecState.push_back(Vector3f(.3 ,0, 0));
    m_vVecState.push_back(Vector3f(0, -.02, 0));
    m_vVecState.push_back(Vector3f(0 ,0, 0));
    m_vVecState.push_back(Vector3f(0, -.03, 0));
    m_vVecState.push_back(Vector3f(0 ,0, 0));
}


std::vector<Vector3f> PendulumSystem::evalF(std::vector<Vector3f> state)
{
    std::vector<Vector3f> f;
    // TODO 4.1: implement evalF
    //  - gravity
    //  - viscous drag
    //  - springs

    Vector3f zero = Vector3f(0,0,0);
    f.push_back(zero);
    f.push_back(zero);

    for (int i = 2; i < state.size(); i+=2){
        Vector3f position = state[i];
        Vector3f velocity = state[i+1];

        // F = mg
        float gForce = PARTICLE_MASS*GRAVITY;
        // only in the y-direction
        Vector3f gravityForce = Vector3f(0,-gForce, 0);

        // F = -k*v
        Vector3f dragForce = -1*DRAG_CONSTANT*velocity;

        Vector3f springForceParent = getSpringForce(position, state[i-2]);

        Vector3f springForceChild = Vector3f(0,0,0);
        if (i+2 < state.size()){
            springForceChild = getSpringForce(position, state[i+2]);
        }
        Vector3f totalForce = gravityForce+dragForce+springForceParent+springForceChild;
        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;

        f.push_back(velocity);
        f.push_back(acceleration);
    }

    return f;
}

Vector3f PendulumSystem::getSpringForce(Vector3f position1, Vector3f position2)
{       
        // −k*(|d| − r)* d/|d|

        // d = p1 - p0
        Vector3f displacement = position1 - position2;
        // |d|
        float absDisplacement = displacement.abs();
        // delta = |d| - r
        float delta = absDisplacement-SPRING_REST_LENGTH;
        Vector3f springForce = -1*SPRING_CONSTANT*(delta)*(displacement.normalized());

        return springForce;

}

// render the system (ie draw the particles)
void PendulumSystem::draw(GLProgram& gl)
{
    const Vector3f PENDULUM_COLOR(0.73f, 0.0f, 0.83f);
    gl.updateMaterial(PENDULUM_COLOR);

    // TODO 4.2, 4.3

    // example code. Replace with your own drawing  code
    Vector3f anchor = m_vVecState[0];
    gl.updateModelMatrix(Matrix4f::translation(anchor));
    drawSphere(0.075f, 10, 10);

    for (int i = 0; i < m_vVecState.size(); i+=2){
        Vector3f pos = m_vVecState[i];
        gl.updateModelMatrix(Matrix4f::translation(pos));
        drawSphere(0.075f, 10, 10);
    }
}

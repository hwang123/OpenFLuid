#include "fluidsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
#include <math.h>
#include "particle.h"
#include "wall.h"

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 5;

const float PARTICLE_RADIUS = .015;
const float PARTICLE_SPACING = .02;

const float H = .0457;
const float mu = 3.5;//0.001;

const float SPRING_CONSTANT = 5; // N/m
const float PARTICLE_MASS = 0.02;//.03; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;

float WPoly6(Vector3f R);
Vector3f WSpiky(Vector3f R);
float WViscosity(Vector3f R);

FluidSystem::FluidSystem()
{
    // back, front, left, right, bottom - respectively
    _walls.push_back(Wall(Vector3f(0,-0.5,-1), Vector3f(0,0,1)));
    _walls.push_back(Wall(Vector3f(0,-0.5,1), Vector3f(0,0,-1)));
    _walls.push_back(Wall(Vector3f(-1,-0.5,0), Vector3f(1,0,0)));
    _walls.push_back(Wall(Vector3f(1,-0.5,0), Vector3f(-1,0,0)));
    _walls.push_back(Wall(Vector3f(0,-1,0), Vector3f(0,1,0)));

    // Initialize m_vVecState with fluid particles. 
    // You can again use rand_uniform(lo, hi) to make things a bit more interesting
    m_vVecState.clear();
    int particleCount = 0;
    for (unsigned i = 0; i < N; i++){
        for (unsigned j = 0; j< N; j++){
            for (unsigned l = 0; l < N; l++){
                float x = i*PARTICLE_SPACING;
                float y = j*PARTICLE_SPACING;
                float z = l*PARTICLE_SPACING;
                // particles evenly spaced
                Vector3f position = Vector3f(x, y, z);
                // all particles stationary
                // Vector3f velocity = Vector3f( rand_uniform(-1,  1), 0, rand_uniform( -1, 1));
                Vector3f velocity = Vector3f(0);

                Particle particle = Particle(particleCount, position, velocity);
                m_vVecState.push_back(particle);
                particleCount += 1;
            }
        }
    }
}


std::vector<Particle> FluidSystem::evalF(std::vector<Particle> state)
{

    std::vector<Particle> f;
    // FLuid Particles undergo the following forces:
    // - gravity
    // - viscous drag
    // - pressure
    // ---Below Not Implemented---
    // - surface tension
    // - collision forces


    // F = mg  -> GRAVITY
    // ----------------------------------------
    float gForce = PARTICLE_MASS*GRAVITY;
    // only in the y-direction
    Vector3f gravityForce = Vector3f(0,-gForce, 0);
    // ----------------------------------------

    for (unsigned i = 0; i < m_vVecState.size(); i+=1){
        Particle& particle = state[i];
        Vector3f position = particle.getPosition();
        Vector3f velocity = particle.getVelocity();
        // compute density of every particle first
        float density_i = 0;
        for (unsigned j = 0; j < state.size(); j+=1) {
            if (j != i) {
                Particle particle_j = state[j];

                Vector3f delta = position - particle_j.getPosition();
                if (delta.absSquared() < H*H) {
                    density_i += PARTICLE_MASS*WPoly6(delta);
                    // cout << WPoly6(delta) << endl;
                }
            }
        }
        particle.density() = density_i;
        // cout << density_i << endl;

    }

    for (unsigned i = 0; i < m_vVecState.size(); i+=1){
        Particle particle = state[i];
        Vector3f position = particle.getPosition();
        Vector3f velocity = particle.getVelocity();
        float density = particle.density();

        // cout << density << endl;

        float pressure = particle.getPressure();


        // compute updated density and gradient of pressure
        // based on all other particles
        Vector3f f_pressure;
        Vector3f f_viscosity;
        for (unsigned j = 0; j < state.size(); j+=1) {
            if (j != i) {
                Particle particle_j = state[j];

                Vector3f delta = position - particle_j.getPosition() - Vector3f(2 * PARTICLE_RADIUS);
                if (delta.absSquared() < H*H) {
                    //  ---------------gradient of pressure computation-----------------

                    // Mueller value: (pi + pj) / 2roj
                    float p_factor = (pressure+particle_j.getPressure()) / (2*particle_j.getDensity()); 
                    // Lecture video value: pi/roi + pj/roj
                    // float p_factor = pressure/density + particle_j.getPressure()/particle_j.getDensity();
                    f_pressure += PARTICLE_MASS*p_factor*WSpiky(delta);
                    //  ---------------viscosity computation-----------------
                    float kernel_distance_viscosity = H-delta.abs();

                    Vector3f v_factor = (particle_j.getVelocity() - velocity) / particle_j.getDensity();

                    Vector3f viscosity_term = PARTICLE_MASS*WViscosity(delta)*v_factor;
                    // cout << "delta " << kernel_constant_pressure << endl;
                    // viscosity_term.print();
                    f_viscosity += viscosity_term;
                    // velocity.print();
                }
            }
        }

        // Total Force
        Vector3f totalForce = (gravityForce +(mu*f_viscosity) + f_pressure)/density;
        // totalForce.print();

        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;


        if (position.y() < -0.95){
            velocity = Vector3f(1,0,1)* velocity;
            acceleration = Vector3f(1,-.5,1)*acceleration;
        }

        if (position.x() < -0.95){
            velocity = Vector3f(-1,1,1)* velocity;
            acceleration = Vector3f(-1,1,1)*acceleration;
        }

        if (position.x() > 0.95){
            velocity = Vector3f(0, velocity.y(), velocity.z());
            acceleration = Vector3f(-1,1,1)*acceleration;
        }

        if (position.z() < -0.95){
            velocity = Vector3f(velocity.x(), velocity.y(), 0);
            acceleration = Vector3f(1,1,-1)*acceleration;
        }

        if (position.z() > 0.95){
            velocity = Vector3f(velocity.x(), velocity.y(), 0);
            acceleration = Vector3f(1,1,-1)*acceleration;
        }

        Particle newParticle = Particle(i, velocity, acceleration);
        f.push_back(newParticle);
    }

    return f;
}

float WPoly6(Vector3f R){
    float constant_term = 315.0 / (64.0 * M_PI * pow(H, 9));
    float kernel_distance_density = pow((H*H - R.absSquared()), 3);
    return kernel_distance_density*constant_term;
}

Vector3f WSpiky(Vector3f R){
    if (R.abs() < .1){
        return Vector3f(0,0,0);
    }
    float constant_term = -45.0 / (M_PI * pow(H, 6));
    Vector3f kernel_distance_pressure = pow((H - R.abs()), 2) * R.normalized();
    return constant_term * kernel_distance_pressure;
}

float WViscosity(Vector3f R){
    float constant_term = 45.0 / (M_PI * pow(H, 6));
    float kernel_distance_viscosity = H-R.abs();
    return constant_term * kernel_distance_viscosity;
}

void FluidSystem::draw(GLProgram& gl)
{
    //TODO 5: render the system 
    //         - ie draw the particles as little spheres
    //         - or draw the springs as little lines or cylinders
    //         - or draw wireframe mesh

    const Vector3f blue(0.0f, 0.0f, 1.0f);

    // EXAMPLE for how to render cloth particles.
    //  - you should replace this code.
    // EXAMPLE: This shows you how to render lines to debug the spring system.
    //
    //          You should replace this code.
    //
    //          Since lines don't have a clearly defined normal, we can't use
    //          a regular lighting model.
    //          GLprogram has a "color only" mode, where illumination
    //          is disabled, and you specify color directly as vertex attribute.
    //          Note: enableLighting/disableLighting invalidates uniforms,
    //          so you'll have to update the transformation/material parameters
    //          after a mode change.
    gl.disableLighting();
    gl.updateModelMatrix(Matrix4f::identity()); // update uniforms after mode change
    // drawBox(Vector3f(0,0,0), 1);

    // not working :(
    gl.updateMaterial(blue);

    for (unsigned i = 0; i < m_vVecState.size(); i+=1){
        Particle p = m_vVecState[i];
        Vector3f pos = p.getPosition();
        gl.updateModelMatrix(Matrix4f::translation(pos));
        drawSphere(PARTICLE_RADIUS, 5, 4);
    }

    gl.enableLighting(); // reset to default lighting model
}


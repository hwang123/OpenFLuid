#include "clothsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
#include "particle.h"

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 10;

const float PARTICLE_SPACING = .08;
const float PARTICLE_RADIUS = .03;

const float SPRING_CONSTANT = 5; // N/m
const float PARTICLE_MASS = .03; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;

ClothSystem::ClothSystem()
{
    // TODO 5. Initialize m_vVecState with cloth particles. 
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
                Vector3f velocity = Vector3f(0, 0, 0);

                Particle particle = Particle(particleCount, position, velocity);
                m_vVecState.push_back(particle);
                particleCount += 1;
            }
        }
    }
}


std::vector<Particle> ClothSystem::evalF(std::vector<Particle> state)
{

    std::vector<Particle> f;
    // TODO 5. implement evalF
    // - gravity
    // - viscous drag
    // - structural springs
    // - shear springs
    // - flexion springs
    // Particle p = Particle(4.0, Vector3f(0,0,0), Vector3f(0,0,0));


    // Gravity and Wind will be independent of the particle in question

    // F = mg  -> GRAVITY
    // ----------------------------------------
    float gForce = PARTICLE_MASS*GRAVITY*0;
    // only in the y-direction
    Vector3f gravityForce = Vector3f(0,-gForce, 0);
    // ----------------------------------------

    for (unsigned i = 0; i < state.size(); i+=1){
        Particle particle = state[i];
        Vector3f position = particle.getPosition();
        Vector3f velocity = particle.getVelocity();

        // Gravity
        Vector3f totalForce =  gravityForce;

        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;

        Particle newParticle = Particle(i, velocity, acceleration);
        f.push_back(newParticle);
    }

    return f;
}

void ClothSystem::draw(GLProgram& gl)
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


#include "clothsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
#include "particle.h"

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 10;

const int W = N;
const int H = N;
const float PARTICLE_SPACING = .1;

const float SPRING_CONSTANT = 5; // N/m
const float PARTICLE_MASS = .03; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;
const float STRUCTURAL_REST_LENGTH = .2; // m
const float SHEAR_REST_LENGTH = STRUCTURAL_REST_LENGTH*sqrt(2); // m
const float FLEXION_REST_LENGTH = STRUCTURAL_REST_LENGTH*2; // m


ClothSystem::ClothSystem()
{
    // TODO 5. Initialize m_vVecState with cloth particles. 
    // You can again use rand_uniform(lo, hi) to make things a bit more interesting
    m_vVecState.clear();
    for (unsigned i = 0; i < W; i++){
        for (unsigned j = 0; j<H; j++){
            for (unsigned l = 0; l < N; l++){
                float x = i*PARTICLE_SPACING;
                float y = j*PARTICLE_SPACING;
                float z = l*PARTICLE_SPACING;
                // particles evenly spaced
                Vector3f particle = Vector3f(x, y, z);
                // all particles stationary
                Vector3f velocity = Vector3f(0, 0, 0);

                m_vVecState.push_back(particle);
                m_vVecState.push_back(velocity);
            }
        }
    }
}


std::vector<Vector3f> ClothSystem::evalF(std::vector<Vector3f> state)
{

    std::vector<Vector3f> f;
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
    float gForce = PARTICLE_MASS*GRAVITY;
    // only in the y-direction
    Vector3f gravityForce = Vector3f(0,-gForce, 0);
    // ----------------------------------------

    for (unsigned i = 0; i < state.size(); i+=2){

        Vector3f position = state[i];
        Vector3f velocity = state[i+1];

        // Gravity, Drag, and Wind
        Vector3f totalForce =  gravityForce;

        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;

        f.push_back(velocity);
        f.push_back(acceleration);
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
    float w = 0.2f;
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

    for (int i = 0; i < m_vVecState.size(); i+=2){
        Vector3f pos = m_vVecState[i];
        gl.updateModelMatrix(Matrix4f::translation(pos));
        drawSphere(0.03f, 7, 4);
    }

    gl.enableLighting(); // reset to default lighting model
}


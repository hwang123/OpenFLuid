#include "fluidsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
#include <math.h>
#include "particle.h"

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 5;

const float PARTICLE_SPACING = .08;
const float PARTICLE_RADIUS = .03;
const int H = 3 * PARTICLE_RADIUS;


const float SPRING_CONSTANT = 5; // N/m
const float PARTICLE_MASS = .03; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;

FluidSystem::FluidSystem()
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


std::vector<Particle> FluidSystem::evalF(std::vector<Particle> state)
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
    float gForce = PARTICLE_MASS*GRAVITY;
    // only in the y-direction
    Vector3f gravityForce = Vector3f(0,-gForce, 0);
    // ----------------------------------------
    float kernel_constant_density = 315.0 / (64.0 * M_PI * pow(H, 9));
    float kernel_constant_pressure = -45.0 / (M_PI * pow(H, 6));

    for (unsigned i = 0; i < state.size(); i+=1){
        Particle particle = state[i];
        Vector3f position = particle.getPosition();
        Vector3f velocity = particle.getVelocity();
        float density = particle.getDensity();
        Vector3f pressure = particle.getPressure();

        // compute updated density and gradient of pressure
        // based on all other particles
        float density_i = 0;
        Vector3f grad_pressure;
        for (unsigned j = 0; j < state.size(); j+=1) {
            Particle particle_j = state[j];
            if (j != i) {
                Vector3f delta = position - particle_j.getPosition();

                // desnity computation
                float kernel_distance_density = pow((H*H - delta.absSquared()), 3);
                density_i += PARTICLE_MASS*kernel_constant_density*kernel_distance_density;

                // gradient of pressure computation
                Vector3f p_factor = pressure/(density*density) + particle_j.getPressure()/(particle_j.getDensity()*particle_j.getDensity());
                Vector3f kernel_distance_pressure = pow((H - delta.absSquared()), 2) * delta / delta.abs();

                grad_pressure += PARTICLE_MASS*p_factor*kernel_constant_pressure*kernel_distance_pressure;
            }


        }


        // Gravity
        Vector3f totalForce =  gravityForce;

        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;

        if (position.y() < -0.95){
            velocity = Vector3f(velocity.x(), 0, velocity.z());
        }

        Particle newParticle = Particle(i, velocity, acceleration);
        newParticle.setDensity(density_i);
        // newParticle.setPressure(density_i);
        f.push_back(newParticle);
    }

    return f;
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

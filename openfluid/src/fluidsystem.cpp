#include "fluidsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>
#include <math.h>
#include <algorithm>    // std::max
#include "particle.h"
#include "wall.h"

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 6;

const float PARTICLE_RADIUS = .015;
const float PARTICLE_SPACING = .0225;

const float H = 0.03;//.0457;
const float Hsquared = H*H;
const float Hfouth = Hsquared*Hsquared;
const float mu = .00081;//0.001;

const float SPRING_CONSTANT = 5; // N/m
const float PARTICLE_MASS = 0.02;//.03; // kg 
const float GRAVITY = 9.8; // m/s
const float DRAG_CONSTANT = .05;

const float WALL_K = 75.0; // wall spring constant
const float WALL_DAMPING = -0.0001; // wall damping constant


float WPoly6(Vector3f R);
Vector3f WSpiky(Vector3f R);
float WViscosity(Vector3f R);

FluidSystem::FluidSystem()
{
    // change box size - make sure this is the same as in main.cpp
    float len = 0.1;
    // back, front, left, right, bottom - respectively
    _walls.push_back(Wall(Vector3f(0,-len/2,-len), Vector3f(0,0,1)));
    _walls.push_back(Wall(Vector3f(0,-len/2,len), Vector3f(0,0,-1)));
    _walls.push_back(Wall(Vector3f(-len,-len/2,0), Vector3f(1,0,0)));
    _walls.push_back(Wall(Vector3f(len,-len/2,0), Vector3f(-1,0,0)));
    _walls.push_back(Wall(Vector3f(0,-len,0), Vector3f(0,1,0)));

    // Initialize m_vVecState with fluid particles. 
    // You can again use rand_uniform(lo, hi) to make things a bit more interesting
    m_vVecState.clear();
    int particleCount = 0;
    for (unsigned i = 0; i < N-1; i++){
        for (unsigned j = 0; j< N-1; j++){
            for (unsigned l = 0; l < N; l++){
                float x = -len + i*PARTICLE_SPACING;
                float y = 0.01 + -len + j*PARTICLE_SPACING;
                float z = -len + l*PARTICLE_SPACING;
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

    // for (unsigned i = 0; i < N; i++){
    //     for (unsigned j = 0; j< N; j++){
    //         float x = -len + i*PARTICLE_SPACING;
    //         float y = .1;
    //         float z = -len + j*PARTICLE_SPACING;
    //         // particles evenly spaced
    //         Vector3f position = Vector3f(x, y, z);
    //         // all particles stationary
    //         // Vector3f velocity = Vector3f( rand_uniform(-1,  1), 0, rand_uniform( -1, 1));
    //         Vector3f velocity = Vector3f(0);

    //         Particle particle = Particle(particleCount, position, velocity);
    //         m_vVecState.push_back(particle);
    //         particleCount += 1;
    //     }
    // }
}


std::vector<Particle> FluidSystem::evalF(std::vector<Particle> state)
{

    std::vector<Particle> f;
    // FLuid Particles undergo the following forces:
    // - gravity
    // - collision forces
    // - viscous drag
    // - pressure
    // ---Below Not Implemented---
    // - surface tension


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
        double density_i = 0;
        for (unsigned j = 0; j < state.size(); j+=1) {
            if (j != i) {
                Particle particle_j = state[j];

                Vector3f delta = position - particle_j.getPosition();
                if (delta.absSquared() < H*H) {
                    density_i += WPoly6(delta);
                    // cout << WPoly6(delta) << endl;
                }
            }
        }
        particle.density() = .001+PARTICLE_MASS*density_i;
        // cout <<particle.density() << endl;

        // if (particle.density() < 50){
        // particle.density() = 500;
        // }
        // cout << density_i << endl;

    }

    for (unsigned i = 0; i < m_vVecState.size(); i+=1){
        Particle particle = state[i];
        Vector3f position = particle.getPosition();
        Vector3f velocity = particle.getVelocity();
        double density = particle.density();
        float pressure = particle.getPressure();

        // compute updated density and gradient of pressure
        // based on all other particles
        Vector3f f_pressure;
        Vector3f f_viscosity;
        for (unsigned j = 0; j < state.size(); j+=1) {
            if (j != i) {
                Particle particle_j = state[j];

                Vector3f delta = position - particle_j.getPosition();// + Vector3f(2 * PARTICLE_RADIUS);
                // cout << "delta" << delta.x() << ","<< delta.y() << ","<< delta.z() << endl;
                if (delta.absSquared() < H*H) {
                    //  ---------------gradient of pressure computation-----------------

                    // sample implementation value: pi/roi2 + pj/roj2
                    float p_factor = (pressure/(density*density)) + (particle_j.getPressure()/(particle_j.density()*particle_j.density())); 

                    // Mueller value: (pi + pj) / 2roj
                    // float p_factor = (pressure+particle_j.getPressure()) / (2*particle_j.getDensity()); 

                    // Lecture video value: pi/roi + pj/roj
                    // float p_factor = pressure/density + particle_j.getPressure()/particle_j.getDensity();
                    // (WSpiky(delta)).print();
                    // delta.print();
                    f_pressure += p_factor*WSpiky(delta);
                    //  ---------------viscosity computation-----------------
                    float kernel_distance_viscosity = H-delta.abs();

                    Vector3f v_factor = (particle_j.getVelocity() - velocity);// / particle_j.getDensity();
                    // cout << WViscosity(delta) << endl;

                    Vector3f viscosity_term = (WViscosity(delta)/ particle_j.getDensity())*v_factor;
                    // cout << "delta " << kernel_constant_pressure << endl;
                    // viscosity_term.print();
                    f_viscosity += viscosity_term;
                    // velocity.print();
                        float dd = WViscosity(delta);

                    // cout << dd << endl;
                    // cout << particle_j.getDensity() << endl;

                    // if ((viscosity_term * mu * PARTICLE_MASS).y() > 20){
                    //     cout << j << endl;
                    //     cout << "alert" << endl;
                    //     cout << particle_j.getDensity() << endl;
                    //     cout << dd <<endl;
                    //     v_factor.print();
                    //     viscosity_term.print();
                    // }
                }
            }
        }

        f_pressure *= -.0051*PARTICLE_MASS;
        f_viscosity *= mu * PARTICLE_MASS;


        // Total Force
        // Vector3f totalForce = (gravityForce +(mu*f_viscosity) + f_pressure)/density;
        // f_pressure.print();
        // f_viscosity.print();
        // cout << density << endl;
        Vector3f f_collision = collisionForce(particle);

        // (f_viscosity/PARTICLE_MASS).print();

        // Total Force
        // f_collision = Vector3f(0, f_collision.y(), 0);
        // (f_viscosity*.001).print();
        Vector3f totalForce = gravityForce + f_collision + f_viscosity + f_pressure;
        // if (f_viscosity.abs() > 15){
        //     cout << "wtf" << endl;
        //            f_viscosity.print();
        // f_collision.print(); 
        // }

        // Vector3f totalForce = gravityForce + f_collision ;//+ f_pressure;
        // cout << "f_collision: " << f_collision.x() << ","<< f_collision.y() << ","<< f_collision.z() << endl;
        // cout << "f_pressure: " << f_pressure.x() << ","<< f_pressure.y() << ","<< f_pressure.z() << endl;
        // cout << "gravityForce: " << gravityForce.x() << ","<< gravityForce.y() << ","<< gravityForce.z() << endl;
        // totalForce.print();

        Vector3f acceleration = .2*totalForce / PARTICLE_MASS;

        Particle newParticle = Particle(i, velocity, acceleration, density);
        f.push_back(newParticle);
    }

    return f;
}

float WPoly6(Vector3f R){
    float constant_term = 315.0 / (64.0 * M_PI * Hfouth * Hfouth * H);
    float factor = Hsquared - R.absSquared();
    float kernel_distance_density = factor * factor * factor;
    return kernel_distance_density*constant_term;
}

Vector3f WSpiky(Vector3f R){
    if (R.abs() < .000001){
        return Vector3f(0,0,0);
    }
    float constant_term = -45.0 / (M_PI * Hfouth * Hsquared);
    float factor = H - R.abs();
    Vector3f kernel_distance_pressure = factor * factor * R.normalized();
    return constant_term * kernel_distance_pressure;
}

float WViscosity(Vector3f R){
    float constant_term = 45.0 / (M_PI * Hfouth * Hsquared);
    float kernel_distance_viscosity = H-R.abs();
    return constant_term * kernel_distance_viscosity /10.0;
}


Vector3f FluidSystem::collisionForce(Particle particle) {
    Vector3f f_collision;
    for (unsigned int i = 0; i < _walls.size(); i++) {

        Wall wall = _walls[i];
        Vector3f wallNormal = wall.getNormal();
        double d = Vector3f::dot(wallNormal, (wall.getPoint() - particle.getPosition())) +  0.008; // particle radius
        if (d > 0.0) {      
            f_collision += WALL_DAMPING * Vector3f::dot(wallNormal, particle.getVelocity()) * wallNormal + WALL_K * wallNormal * d;
        }
    }
    return f_collision;
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
        drawSphere(.009, 10, 4);
    }

    gl.enableLighting(); // reset to default lighting model
}


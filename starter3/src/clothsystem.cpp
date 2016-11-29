#include "clothsystem.h"
#include "camera.h"
#include "vertexrecorder.h"
#include <iostream>

using namespace std;

 // your system should at least contain 8x8 particles.
const int N = 11;

const int W = N;
const int H = N;
const float PARTICLE_SPACING = .2;

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
    for (int i = 0; i < W; i++){
        for (int j = 0; j<H; j++){
            // variables in case we want to change orientation
            float x = i*PARTICLE_SPACING;
            float y = j*PARTICLE_SPACING;
            // particles evenly spaced
            Vector3f particle = Vector3f(x, y, 2);
            // all particles stationary
            Vector3f velocity = Vector3f(0, 0, 0);

            m_vVecState.push_back(particle);
            m_vVecState.push_back(velocity);
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



    // Gravity and Wind will be independent of the particle in question

    // F = mg  -> GRAVITY
    // ----------------------------------------
    float gForce = PARTICLE_MASS*GRAVITY;
    // only in the y-direction
    Vector3f gravityForce = Vector3f(0,-gForce, 0);
    // ----------------------------------------

    // WIND FORCE
    // only in the z & x directions
    float w1 = (rand() % 20)/10.0 - 1;

    Vector3f windForce = Vector3f(0, 0, 2*w1);
    // ----------------------------------------

    for (int i = 0; i < state.size(); i+=2){

        // cout << "DEBUG" << endl;

        Vector3f position = state[i];
        Vector3f velocity = state[i+1];

        // F = -k*v -> VISCOUS DRAG 
        // ----------------------------------------
        Vector3f dragForce = -1*DRAG_CONSTANT*velocity;
        // ----------------------------------------

        Vector2f rowColumn = rowColumnOf(i/2);


        // STRUCTURAL SPRINGS
        // ----------------------------------------
        vector<Vector3f> structuralSpringForces;

        vector<int> structuralIndices = getStucturalNeighbors(rowColumn);

        for (int j = 0; j<structuralIndices.size(); j++){
            int structureIndex = structuralIndices[j];
            Vector3f structuralPosition = state[structureIndex];
            Vector3f structuralSpringForce = getSpringForce(position, structuralPosition, STRUCTURAL_REST_LENGTH);

            structuralSpringForces.push_back(structuralSpringForce);
        }
        // ----------------------------------------


        // SHEAR SPRINGS
        // ----------------------------------------
        vector<Vector3f> shearSpringForces;

        vector<int> shearIndices = getShearNeighbors(rowColumn);

        for (int j = 0; j<shearIndices.size(); j++){
            int shearIndex = shearIndices[j];
            Vector3f shearPosition = state[shearIndex];
            Vector3f shearSpringForce = getSpringForce(position, shearPosition, SHEAR_REST_LENGTH);

            shearSpringForces.push_back(shearSpringForce);
        }
        // ----------------------------------------

        // FLEXION SPRINGS
        // ----------------------------------------
        vector<Vector3f> flexionSpringForces;

        vector<int> flexionIndices = getFlexionNeighbors(rowColumn);

        for (int j = 0; j<flexionIndices.size(); j++){
            int flexionIndex = flexionIndices[j];
            Vector3f flexionPosition = state[flexionIndex];
            Vector3f flexionSpringForce = getSpringForce(position, flexionPosition, FLEXION_REST_LENGTH);

            flexionSpringForces.push_back(flexionSpringForce);
        }
        // ----------------------------------------


        // Now sum up all the forces that we have:

        // Gravity, Drag, and Wind
        Vector3f totalForce = gravityForce + dragForce + windForce;

        // As well as all the spring forces
        for (int j = 0; j<structuralSpringForces.size(); j++){
            totalForce = totalForce + structuralSpringForces[j];
        }
        for (int j = 0; j<shearSpringForces.size(); j++){
            totalForce = totalForce + shearSpringForces[j];
        }
        for (int j = 0; j<flexionSpringForces.size(); j++){
            totalForce = totalForce + flexionSpringForces[j];
        }

        Vector3f acceleration = (1.0/PARTICLE_MASS)*totalForce;

        // ugly, but yeah, this grabs the top corners of the cloth to keep them fixed
        if (i+2 == m_vVecState.size() or i+2 + (2*H*(W-1)) ==m_vVecState.size()){
            f.push_back(Vector3f(0));
            f.push_back(Vector3f(0));
        }
        else{
            f.push_back(velocity);
            f.push_back(acceleration);
        }
    }

    return f;
}

vector<int> ClothSystem::getStucturalNeighbors(Vector2f row_column)
{
    int row = row_column.x();
    int column = row_column.y();

    vector<int> returnIndices;

    if (column + 1 < W){
        int rightIndex = indexOf(row, column + 1);
        returnIndices.push_back(rightIndex);
    }

    if (column - 1 >= 0){
        int leftIndex = indexOf(row, column - 1);
        returnIndices.push_back(leftIndex);

    }

    if (row - 1 >= 0){
        int bottomIndex = indexOf(row - 1, column);
        returnIndices.push_back(bottomIndex);

    }

    if (row + 1 < H){
        int topIndex = indexOf(row + 1, column);
        returnIndices.push_back(topIndex);
    }

    return returnIndices;
}

vector<int> ClothSystem::getShearNeighbors(Vector2f row_column)
{
    int row = row_column.x();
    int column = row_column.y();

    vector<int> returnIndices;

    if (column + 1 < W){
        if (row - 1 >= 0){
            int topRight = indexOf(row-1, column + 1);
            returnIndices.push_back(topRight);


        }
        if (row + 1 < H){
            int bottomRight = indexOf(row+1, column + 1);
            returnIndices.push_back(bottomRight);
        }
    }

    if (column - 1 >= 0){
        if (row - 1 >= 0){
            int topLeft = indexOf(row-1, column - 1);
            returnIndices.push_back(topLeft);
        }
        if (row + 1 < H){
            int bottomLeft = indexOf(row+1, column - 1);
            returnIndices.push_back(bottomLeft);
        }
    }

    return returnIndices;
}

vector<int> ClothSystem::getFlexionNeighbors(Vector2f row_column)
{
    int row = row_column.x();
    int column = row_column.y();

    vector<int> returnIndices;

    if (column + 2 < W){
        int rightIndex = indexOf(row, column + 2);
        returnIndices.push_back(rightIndex);
    }

    if (column - 2 >= 0){
        int leftIndex = indexOf(row, column - 2);
        returnIndices.push_back(leftIndex);

    }

    if (row - 2 >= 0){
        int bottomIndex = indexOf(row - 2, column);
        returnIndices.push_back(bottomIndex);

    }

    if (row + 2 < H){
        int topIndex = indexOf(row + 2, column);
        returnIndices.push_back(topIndex);
    }

    return returnIndices;
}


Vector3f ClothSystem::getSpringForce(Vector3f position1, Vector3f position2, float restLength)
{       
        // −k*(|d| − r)* d/|d|

        // d = p1 - p0
        Vector3f displacement = position1 - position2;
        // |d|
        float absDisplacement = displacement.abs();
        // delta = |d| - r
        float delta = absDisplacement-restLength;
        Vector3f springForce = -1*SPRING_CONSTANT*(delta)*(displacement.normalized());

        return springForce;

}

int ClothSystem::indexOf(int row, int column)
{
    return 2*(column + row*W);
}

Vector2f ClothSystem::rowColumnOf(int index)
{
    int column = index%W;
    int row = index/W;
    Vector2f location = Vector2f(row, column);
    return location;
}

void ClothSystem::draw(GLProgram& gl)
{
    //TODO 5: render the system 
    //         - ie draw the particles as little spheres
    //         - or draw the springs as little lines or cylinders
    //         - or draw wireframe mesh

    const Vector3f CLOTH_COLOR(0.9f, 0.9f, 0.9f);
    gl.updateMaterial(CLOTH_COLOR);

    // EXAMPLE for how to render cloth particles.
    //  - you should replace this code.
    float w = 0.2f;
    Vector3f O(0.4f, 1, 0);


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
    VertexRecorder rec;

    for (int i = 2; i < m_vVecState.size()-2*W; i+=2){
        Vector3f pos = m_vVecState[i];
        // gl.updateModelMatrix(Matrix4f::translation(O + pos));

        Vector3f neighbor1 = m_vVecState[i-2];
        Vector3f neighbor2 = m_vVecState[i-2+2*W];

        Vector3f oneHalf = neighbor2-pos;
        Vector3f otherHalf = neighbor1-pos;

        Vector3f direction = Vector3f::cross(otherHalf, oneHalf);
        Vector3f normal = 1*direction.normalized();

        rec.record(O + pos, normal, CLOTH_COLOR);
        rec.record(O + neighbor2, normal, CLOTH_COLOR);
        rec.record(O + neighbor1, normal, CLOTH_COLOR);


        Vector3f neighbor3 = m_vVecState[i+2*W];

        Vector3f lastHalf = neighbor3-pos;

        Vector3f direction2 = Vector3f::cross(lastHalf, oneHalf);
        Vector3f normal2 = 1*direction2.normalized();

        rec.record(O + pos, normal2, CLOTH_COLOR);
        rec.record(O + neighbor3, normal2, CLOTH_COLOR);
        rec.record(O + neighbor2, normal2, CLOTH_COLOR);

        // drawSphere(0.075f, 10, 10);
    }

    rec.draw();

    // for (int i = 0; i < m_vVecState.size(); i+=2){
    //     Vector3f pos = m_vVecState[i];

    //     Vector2f dimensions = rowColumnOf(i/2);
    //     vector<int> indices = getStucturalNeighbors(dimensions);
    //     for (int j = 0; j < indices.size(); j++){
    //         int index = indices[j];
    //         Vector3f pos2 = m_vVecState[index];
    //         rec.record(O + pos, CLOTH_COLOR);
    //         rec.record(O + pos2, CLOTH_COLOR);
    //     } 

    // }
    // glLineWidth(3.0f);
    // rec.draw(GL_LINES);

    gl.enableLighting(); // reset to default lighting model
    // EXAMPLE END
}


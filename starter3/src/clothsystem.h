#ifndef CLOTHSYSTEM_H
#define CLOTHSYSTEM_H

#include <vector>

#include "particlesystem.h"

class ClothSystem : public ParticleSystem
{
    ///ADD MORE FUNCTION AND FIELDS HERE
public:
    ClothSystem();

    // evalF is called by the integrator at least once per time step
    std::vector<Vector3f> evalF(std::vector<Vector3f> state) override;

    // draw is called once per frame
    void draw(GLProgram& ctx);

    // inherits
    // std::vector<Vector3f> m_vVecState;

    int indexOf(int row, int column);
    std::vector<int> getStucturalNeighbors(Vector2f row_column);
    std::vector<int> getShearNeighbors(Vector2f row_column);
    std::vector<int> getFlexionNeighbors(Vector2f row_column);
    std::vector<Vector3f> getSpringForces(Vector2f row_column, Vector3f position, std::vector<Vector3f> state);
    Vector3f getSpringForce(Vector3f position1, Vector3f position2, float restLength);
    Vector2f rowColumnOf(int index);
};


#endif

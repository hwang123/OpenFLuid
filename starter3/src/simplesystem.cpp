#include "simplesystem.h"

#include "camera.h"
#include "vertexrecorder.h"


SimpleSystem::SimpleSystem()
{
    m_vVecState.clear();
    m_vVecState.push_back(Vector3f(1,0,0));
}

std::vector<Vector3f> SimpleSystem::evalF(std::vector<Vector3f> state)
{
    std::vector<Vector3f> f;

    Vector3f pos = state[0];
    Vector3f eval = Vector3f(-pos.y(), pos.x(), 0);
    f.push_back(eval);

    return f;
}

// render the system (ie draw the particles)
void SimpleSystem::draw(GLProgram& gl)
{

    // TODO 3.2: draw the particle. 
    //           we provide code that draws a static sphere.
    //           you should replace it with your own
    //           drawing code.
    //           In this assignment, you must manage two
    //           kinds of uniforms before you draw
    //            1. Update material uniforms (color)
    //            2. Update transform uniforms
    //           GLProgram is a helper object that has
    //           methods to set the uniform state.

    const Vector3f PARTICLE_COLOR(0.4f, 0.7f, 1.0f);
    gl.updateMaterial(PARTICLE_COLOR);
    Vector3f pos = m_vVecState[0];
    gl.updateModelMatrix(Matrix4f::translation(pos));
    drawSphere(0.075f, 10, 10);
}

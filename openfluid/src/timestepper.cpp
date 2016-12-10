#include "timestepper.h"

#include <cstdio>
#include <iostream>
#include "particle.h"

using namespace std;

void ForwardEuler::takeStep(ParticleSystem* particleSystem, float stepSize)
{
   // X(t + h) =X + hf(X, t)

	vector<Particle> new_state;

	vector<Particle> X = particleSystem->getState();
	vector<Particle> X_prime = particleSystem->evalF(X);

	Vector3f newPos = Vector3f();
	Vector3f newVel = Vector3f();

	for (int i = 0; i < X.size(); i++){
		Particle particle = X[i];

		newPos = particle.getPosition() + stepSize*X_prime[i].getPosition();
		newVel = particle.getVelocity() + stepSize*X_prime[i].getVelocity();

		Particle new_particle = Particle(i, newPos, newVel);

		new_state.push_back(new_particle);
	}

	particleSystem->setState(new_state);
}

void Trapezoidal::takeStep(ParticleSystem* particleSystem, float stepSize)
{
//  	// f0 =f(X, t)
// 	// f1 =f(X + hf0, t + h)
// 	// X(t + h) =X + (h/2)*(f0 + f1)

// 	vector<Vector3f> final_state;
// 	vector<Vector3f> intermediaryState;

// 	vector<Vector3f> X = particleSystem->getState();
// 	vector<Vector3f> f0 = particleSystem->evalF(X);

// 	for (int i =0; i<f0.size(); i++){
// 		Vector3f oldParticle = X[i];
// 		Vector3f particle = f0[i];
// 		// h*f0
// 		Vector3f step = stepSize*particle;
// 		// X + h*f0
// 		Vector3f finalParticle = oldParticle+step;
// 		intermediaryState.push_back(finalParticle);
	// }

// 	// f1 = X + h*f0
// 	vector<Vector3f> f1 = particleSystem->evalF(intermediaryState);

// 	for (int j = 0; j<f1.size(); j++){
// 		Vector3f startState = X[j];
// 		Vector3f particleF0 = f0[j];
// 		Vector3f particleF1 = f1[j];

// 		// (h/2)*(f0 + f1)
// 		Vector3f particleSum = particleF0+particleF1;
// 		particleSum = (stepSize/2.0)*particleSum;

// 		Vector3f finalParticle = startState+particleSum;

// 		final_state.push_back(finalParticle);
// 	}

// 	particleSystem->setState(final_state);
}


void RK4::takeStep(ParticleSystem* particleSystem, float stepSize)
{
	// LEAP FROG INTERGRATOR

	vector<Particle> new_state;

	vector<Particle> X = particleSystem->getState();
	vector<Particle> X_prime = particleSystem->evalF(X);

	Vector3f newPos = Vector3f();
	Vector3f newVel = Vector3f();

	for (int i = 0; i < X.size(); i++){
		Particle particle = X[i];

		newPos = particle.getPosition() + stepSize*X_prime[i].getPosition() + 0.5*stepSize*stepSize*X_prime[i].getVelocity();
		newVel = particle.getVelocity() + 0.5*stepSize*(X_prime[i].getVelocity()+X_prime[i].getVelocity());

		Particle new_particle = Particle(i, newPos, newVel);

		new_state.push_back(new_particle);
	}

	particleSystem->setState(new_state);

}

// 	vector<Vector3f> X = particleSystem->getState();
// 	vector<Vector3f> k1 = particleSystem->evalF(X);
// 	vector<Vector3f> finalState;

// 	vector<Vector3f> intermediaryState1;
// 	for (int i = 0; i<X.size(); i++){
// 		Vector3f obj = X[i];
// 		Vector3f k2_obj = obj+(k1[i]*(stepSize/2.0));
// 		intermediaryState1.push_back(k2_obj);
// 	}

// 	vector<Vector3f> k2 = particleSystem->evalF(intermediaryState1);

// 	vector<Vector3f> intermediaryState2;
// 	for (int i = 0; i<X.size(); i++){
// 		Vector3f obj = X[i];
// 		Vector3f k3_obj = obj+(k2[i]*(stepSize/2.0));
// 		intermediaryState2.push_back(k3_obj);
// 	}

// 	vector<Vector3f> k3 = particleSystem->evalF(intermediaryState2);

// 	vector<Vector3f> intermediaryState3;
// 	for (int i = 0; i<X.size(); i++){
// 		Vector3f obj = X[i];
// 		Vector3f k4_obj = obj+(k3[i]*stepSize);
// 		intermediaryState3.push_back(k4_obj);
// 	}	

// 	vector<Vector3f> k4 = particleSystem->evalF(intermediaryState3);

// 	for (int i = 0; i<X.size(); i++){
// 		Vector3f obj = X[i];
// 		Vector3f k1_part = k1[i];
// 		Vector3f k2_part = k2[i];
// 		Vector3f k3_part = k3[i];
// 		Vector3f k4_part = k4[i];

// 		Vector3f newStateObj = obj + (stepSize/6.0)*(k1_part+2*k2_part+2*k3_part+k4_part);

// 		finalState.push_back(newStateObj);

// 	}	

// 	particleSystem->setState(finalState);
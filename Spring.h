#pragma once

#include <Eigen/core>

using namespace Eigen;

class Spring {
	double initialLength; // in meters
	double springConstant; //K  in Newton/meter

	//Damper /// dont use this yet
	double dampingCoeffecient;

	/// dont need to store this
	//RowVector3d partical1Position;
	//RowVector3d partical2Position;

	int partical1Indice;
	int partical2Indice;



public:
	double invMass1;
	double invMass2;

	Spring(int indice1, int indice2, double invMass1, double invMass2, double initLength, double K, double C = 1) :
		partical1Indice(indice1), partical2Indice(indice2),
		invMass1(invMass1), invMass2(invMass2),
		initialLength(initLength), springConstant(K), dampingCoeffecient(C){

	}

	double getForce(VectorXd rawX) {
		return getForce(rawX[partical1Indice], rawX[partical2Indice]);
	}
	double getForce(double partical1Position, double partical2Position) {
		return -springConstant *(partical1Position - partical2Position - initialLength);
	}

	double getImpulse(double partical1Position, double partical2Position, double timeStep) {
		return getForce(partical1Position, partical2Position) * timeStep;
	}

	double getImpulse(VectorXd rawX, double timeStep) {
		return getForce(rawX) * timeStep;
	}

/*	RowVector3d dampSpringVelocity(RowVector3d speed1, RowVector3d speed2) {
		return -dampingCoeffecient*(speed1 - speed2);
	}*/

	/// don't use this yet
	double dampSpringForce(double speed1, double speed2) {
		return -dampingCoeffecient*(speed1 - speed2);
	}

	double dampSpringImpulse(VectorXd rawVel, double timeStep) {
		return dampSpringForce(rawVel[partical1Indice], rawVel[partical2Indice]) * timeStep;
	}

	int getParticleIndice1() {
		return partical1Indice;
	}

	int getParticleIndice2() {
		return partical2Indice;
	}
};
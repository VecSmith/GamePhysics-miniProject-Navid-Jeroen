#pragma once

#include <Eigen/core>

using namespace Eigen;

extern double SpringConstant;

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
		return -SpringConstant *(partical1Position - partical2Position - initialLength);
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
	double dampSpringForce(double impulse) {
		double vel1 = impulse * invMass1;
		double vel2 = -impulse * invMass2;
		return -dampingCoeffecient*(vel1 - vel2);
	}

	/// don't use this yet
	double dampSpringForce(double speed1, double speed2) {
		return -dampingCoeffecient*(speed1 - speed2);
	}

	double dampSpringImpulse(VectorXd rawVel, double timeStep) {
		return dampSpringForce(rawVel[partical1Indice], rawVel[partical2Indice]) * timeStep;
	}

	double getTotalImpulse(VectorXd rawX, double timeStep) {
		// mass * acceleration + C * velocity + K * position = 0
		// ma + Cv + Kx = 0;
		// m(d^2 x)/(dt^2) + c dx/dt + k x = 0
		//  homogeneous second order differential equation
		// x = e^lambda
		// m lambda^2 + c lambda + k = 0
		// abc formula
		// lamda = (-c + or - sqrt(c^2 - 4mk) ) / (2m)

		/// uhm 2 masses?
	/*	double useless1 = dampingCoeffecient*dampingCoeffecient - 4 / invMass1* springConstant;
		if (useless1 > 0) {
			cout << "hey" << endl;
		}
		double useless = sqrt(dampingCoeffecient*dampingCoeffecient - 4 / invMass1* springConstant);
		double lambdaPlus = (-dampingCoeffecient + sqrt(dampingCoeffecient*dampingCoeffecient - 4 / invMass1* springConstant)) * (2 * invMass1);
		double lambdaMinus = (-dampingCoeffecient - sqrt(dampingCoeffecient*dampingCoeffecient - 4 / invMass1* springConstant)) * (2 * invMass1);

		
		///shouldnt just take lambdaPlus
		
		double x = exp(lambdaPlus);
		// v_avg = \Delta s / \Delta t
		double velocity = x / timeStep;
		return velocity / invMass1;*/
	/*	double mass1 = 1 / invMass1;
		double mass2 = 1 / invMass2;
		double mass = (mass1 * mass2) / (mass1 + mass2);

		double angular = sqrt(SpringConstant/mass);
		double dampingRatio = dampingCoeffecient / (2 * sqrt(mass*SpringConstant));

		double omega = angular * dampingRatio;
		double alpha = angular * sqrt(1.0 - dampingRatio * dampingRatio);

		double expValue = exp(-omega * timeStep);
		double cosValue = cos(alpha * timeStep);
		double sinValue = sin(alpha * timeStep);

		double pos1 = rawX[partical1Indice];
		double pos2 = rawX[partical2Indice];

		double length = pos1 - pos2;
		double pos = (initialLength + expValue*length) / alpha;
		expValue*((c1*omega - c2*alpha)*cosValue +	(c1*alpha + c2*omega)*sinValue);*/
	}

	int getParticleIndice1() {
		return partical1Indice;
	}

	int getParticleIndice2() {
		return partical2Indice;
	}
};
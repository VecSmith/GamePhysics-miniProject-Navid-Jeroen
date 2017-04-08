//
//  constraints.h
//  practical2
//
//  Created by Amir Vaxman on 14/03/2017.
//
//

#ifndef constraints_h
#define constraints_h

#include <Eigen/core>

using namespace Eigen;
using namespace std;

extern double RigidityAllowance;

//constraint types
typedef enum ConstraintType{ATTACHMENT, RIGIDITY, COLLISION, BARRIER, ATTACHMENTSTATIC} ConstraintType;

class Constraint{
public:

    VectorXi particleIndices;  //list of participating indices
    double currValue;                       //the current value of the constraint
    VectorXd currGradient;               //practicleIndices-sized.
    MatrixXd invMassMatrix;              //M^{-1} matrix
    VectorXd radii;

    double refValue;                    //reference value to compare against. Can be rest length, dihedral angle, barrier, etc.
    double stiffness;
    ConstraintType constraintType;  //the type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!

    Constraint(const ConstraintType _constraintType, const VectorXi& _particleIndices, const VectorXd& _radii, const VectorXd& invMasses, const double _refValue, const double _stiffness):constraintType(_constraintType), refValue(_refValue), stiffness(_stiffness)
    {
        currValue=0.0;
        particleIndices=_particleIndices;
        currGradient=VectorXd::Zero(particleIndices.size());
        invMassMatrix=invMasses.asDiagonal();
        radii=_radii;
    }

    Constraint(const ConstraintType _constraintType, const int& particleIndex, const double& radius, const double& invMass, const double _refValue, const double _stiffness):constraintType(_constraintType), refValue(_refValue), stiffness(_stiffness)
    {
        currValue=0.0;
        particleIndices.resize(1); particleIndices(0)=particleIndex;
        currGradient=VectorXd::Zero(particleIndices.size());
        invMassMatrix.resize(1,1); invMassMatrix(0,0)=invMass;
        radii.resize(1); radii(0)=radius;
    }

    ~Constraint(){}

    //updating the value and the gradient vector with given values in xyzxyzxyz format
    void updateValueGradient(const VectorXd& currPos){
        switch (constraintType){

            case ATTACHMENT:{
				RowVector3d ParticleCenter1 = RowVector3d(currPos(0), currPos(1), currPos(2));
				RowVector3d ParticleCenter2 = RowVector3d(currPos(3), currPos(4), currPos(5));
				RowVector3d ConnectorVector = ParticleCenter1 - ParticleCenter2;
				currValue = ConnectorVector.norm() - refValue - (RigidityAllowance / refValue);
				ConnectorVector = ConnectorVector.normalized();
				currGradient(0) = ConnectorVector(0);
				currGradient(1) = ConnectorVector(1);
				currGradient(2) = ConnectorVector(2);
				currGradient(3) = -ConnectorVector(0);
				currGradient(4) = -ConnectorVector(1);
				currGradient(5) = -ConnectorVector(2);

                break;
            }

            case RIGIDITY:
            {

                RowVector3d ParticleCenter1 = RowVector3d( currPos(0), currPos(1), currPos(2) );
                RowVector3d ParticleCenter2 = RowVector3d( currPos(3), currPos(4), currPos(5) );
                RowVector3d ConnectorVector = ParticleCenter1 - ParticleCenter2;
                double Range = ( RigidityAllowance != 0 ) ? ( refValue / RigidityAllowance ) : (0);
                currValue = ConnectorVector.norm() - refValue - Range;
                ConnectorVector = ConnectorVector.normalized();
                currGradient(0) = ConnectorVector(0);
                currGradient(1) = ConnectorVector(1);
                currGradient(2) = ConnectorVector(2);
                currGradient(3) = -ConnectorVector(0);
                currGradient(4) = -ConnectorVector(1);
                currGradient(5) = -ConnectorVector(2);

                break;
            }

            case COLLISION:
            {
                RowVector3d ParticleCenter1 = RowVector3d( currPos(0), currPos(1), currPos(2) );
                RowVector3d ParticleCenter2 = RowVector3d( currPos(3), currPos(4), currPos(5) );
                RowVector3d ConnectorVector = ParticleCenter1 - ParticleCenter2;
                currValue = ConnectorVector.norm() - ( radii(0) + radii(3) );
                ConnectorVector = ConnectorVector.normalized();
                currGradient(0) = -ConnectorVector(0) * currValue;
                currGradient(1) = -ConnectorVector(1) * currValue;
                currGradient(2) = -ConnectorVector(2) * currValue;
                currGradient(3) = ConnectorVector(0) * currValue;
                currGradient(4) = ConnectorVector(1) * currValue;
                currGradient(5) = ConnectorVector(2) * currValue;
                break;
            }

            case BARRIER:
            {
                currValue = currPos(0) - refValue;
                currGradient(0) = 1.0;
				int debug = 0;
                break;
            }

			/*case SPRING:
			{
				double length = (currPos(0) - currPos(1));
				currValue = length - refValue; //length - initial length
				// if currValue < 0 than the init and max length hasn't been reached
				//		how about a min value?
				// if currValue == 0 max length
				// if currvalue > 0 bigger than max length so needs to be fixed
				currGradient(0) = 1.0;
				currGradient(1) = -1.0;

				/// shouldnt be here
				double K = stiffness;
				double force = -K * currValue;
				currGradient(0) = 1;
				currGradient(1) = -1; // --K

				break;
			}*/

			case ATTACHMENTSTATIC:
			{
				//use this as an example on how to fill these fields
				currValue = (currPos(0) - currPos(1)) - refValue;
				currGradient(0) = 1.0;
				currGradient(1) = -1.0;
				break;
			}
        }
    }


    //computes the position differences to resolve the constraint
    void resolveConstraint(const VectorXd& currPos, VectorXd& posDiffs){
        updateValueGradient(currPos);
        if (((constraintType==COLLISION)||(constraintType==BARRIER))&&(currValue>=0.0)){
            //if constraint is anyhow valid, nothing happens
            posDiffs=VectorXd::Zero(particleIndices.size());
            return;
        }
		if (currValue == 0.0) {
			return;
		}

        //compute posDiffs so that C(currPos+posDiffs) ~= C(currPos)+grad(C)*posdiffs=0, using the lagrange multiplier s.t.
        //Lagrange multiplier lambda holds posdiffs=lambda*invMassMatrix*grad(C) as taught in class
        //don't forget to call updateValueGradient() to get the most update values

		// C(currPos) + lambda * gradient(C)^transposed * inverseMass * gradient(C) = 0
		// lambda = -C(currPos) / (gradient(C)^transposed * inverseMass * gradient(C))
		// I think C(currPos) equals currValue after an  updateValueGradient(currPos); call so than
		//currValue + lambda * currGradient.transpose() * invMassMatrix * currGradient;

		double lambda = -currValue / (currGradient.transpose() * invMassMatrix * currGradient);
		posDiffs = lambda * invMassMatrix * currGradient * stiffness;

		/*for (int i = 0; i < particleIndices.size(); i++) {
			double result = currValue + currGradient(i)*posDiffs(i); // = 0
			updateValueGradient(currPos + posDiffs);
 			assert(result > -0.001 && result < 0.001);
			assert(currValue > -0.001 && currValue < 0.001);
		}*/

    }
};



#endif /* constraints_h */

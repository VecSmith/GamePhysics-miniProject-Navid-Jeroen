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

//constraint types
typedef enum ConstraintType{ATTACHMENT, RIGIDITY, COLLISION, BARRIER} ConstraintType;

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
                //use this as an example on how to fill these fields
                currValue=(currPos(0)-currPos(1))-refValue;
                currGradient(0)=1.0;
                currGradient(1)=-1.0;
                break;
            }

            case RIGIDITY:{

				//TODO within an object
                double DeltaX = currPos(0) - currPos(3);
                double DeltaY = currPos(1) - currPos(4);
                double DeltaZ = currPos(2) - currPos(5);
                double EdgeLength = sqrt( ( DeltaX * DeltaX ) + ( DeltaY * DeltaY ) + ( DeltaZ * DeltaZ ) );
                currValue = EdgeLength - refValue;
                currGradient(0) = DeltaX/EdgeLength;
                currGradient(1) = DeltaY/EdgeLength;
                currGradient(2) = DeltaZ/EdgeLength;
                currGradient(3) = -DeltaX/EdgeLength;
                currGradient(4) = -DeltaY/EdgeLength;
                currGradient(5) = -DeltaZ/EdgeLength;
                if (currValue < 0)
                    int j=0;
                break;
            }

            case COLLISION:{

				//TODO particles of different objects colliding

                break;
            }

            case BARRIER:{

				//TODO for the platform
				//currValue = barrierLocation - randomOtherThing - referenceDistance; // why a double when we should have xyz values how should we know where the barrier counts?
				//currValue =  - refValue; //

                currValue = currPos(0)-refValue;
                currGradient(0)=1.0;
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

		if (constraintType == RIGIDITY && abs(currValue) < 0.01)
		{
			posDiffs = VectorXd::Zero(particleIndices.size());
			return;
		}

        //compute posDiffs so that C(currPos+posDiffs) ~= C(currPos)+grad(C)*posdiffs=0, using the lagrange multiplier s.t.
        //Lagrange multiplier lambda holds posdiffs=lambda*invMassMatrix*grad(C) as taught in class
        //don't forget to call updateValueGradient() to get the most update values
        /*******************
         TODO
         *******************/
		//
		//C(currPos) + lambda * gradient(C)^transposed * inverseMass * gradient(C) = 0
		// lambda = -C(currPos) / (gradient(C)^transposed * inverseMass * gradient(C))
		// I think C(currPos) equals currValue after an  updateValueGradient(currPos); call so than

		//currValue + lambda * currGradient.transpose() * invMassMatrix * currGradient;
		double lambda = -currValue / (currGradient.transpose() * invMassMatrix * currGradient);
		posDiffs = lambda * invMassMatrix * currGradient;
		//currValue + currGradient*posDiff; // = 0

		//posdiffs=lambda*invMassMatrix*grad(C)
    }
};



#endif /* constraints_h */

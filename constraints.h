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
        violated=false;
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
                /*******************
                 TODO
                 *******************/
                break;
            }
                
            case COLLISION:{
                /*******************
                 TODO
                 *******************/
                
                break;
            }
                
            case BARRIER:{
                /*******************
                 TODO
                 *******************/
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
        
        //compute posDiffs so that C(currPos+posDiffs) ~= C(currPos)+grad(C)*posdiffs=0, using the lagrange multiplier s.t.
        //Lagrange multiplier lambda holds posdiffs=lambda*invMassMatrix*grad(C) as taught in class
        //don't forget to call updateValueGradient() to get the most update values
        /*******************
         TODO
         *******************/
    }
};



#endif /* constraints_h */

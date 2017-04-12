#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include "constraints.h"
#include "auxfunctions.h"
#include <igl/per_vertex_normals.h>
#include <igl/edge_topology.h>
#include "volInt.h"

#include "Spring.h"

using namespace Eigen;
using namespace std;

//NOTE default mathematical and world values
#define PI 3.14159265358979323846264338327950
#define GravityAcceleration 9.80665

extern double TirePressure;
extern double RigidityAllowance;
double TireCenterX = 0,TireCenterY = 12.5,TireCenterZ = 0;

void ConstructEFi(const MatrixXi& FE, const MatrixXi& EF, MatrixXi& EFi, MatrixXd& FESigns)
{

    EFi=MatrixXi::Constant(EF.rows(), 2,-1);
    FESigns=MatrixXd::Zero(FE.rows(),FE.cols());
    for (int i=0;i<EF.rows();i++)
        for (int k=0;k<2;k++){
            if (EF(i,k)==-1)
                continue;

            for (int j=0;j<3;j++)
                if (FE(EF(i,k),j)==i)
                    EFi(i,k)=j;
        }


    //doing edge signs
    for (int i=0;i<EF.rows();i++){
        if (EFi(i,0)!=-1) FESigns(EF(i,0),EFi(i,0))=1.0;
        if (EFi(i,1)!=-1) FESigns(EF(i,1),EFi(i,1))=-1.0;
    }


}


//the class the contains each individual rigid objects and their functionality
class Mesh{
public:

    //geometry
    MatrixXd origX;   //original particle positions - never change this!
    MatrixXd prevX;   //the previous time step positions
    MatrixXd currX;   //current particle positions
    MatrixXi T;       //the triangles of the original mesh
    MatrixXd currNormals;
    VectorXd invMasses;   //inverse masses of particles, computed (autmatically) as 1.0/(density * particle area)
    VectorXd radii;      //radii of particles
    int rawOffset;  //the raw index offset of the x value of V from the beginning of the 3*|V| particles


    //kinematics
    bool isFixed;  //is the object immobile
    bool isTire = false;
    double rigidity;  //how much the mesh is really rigid
    MatrixXd currVel;     //velocities per particle. Exactly in the size of origV
    vector<Constraint> meshConstraints;  //the systematic constraints in the mesh (i.e., rigidity)
    MatrixXd currImpulses;


    //Quick-reject checking collision between mesh bounding boxes.
    //Does not need updating

    bool isBoxCollide(const Mesh& m2){
        RowVector3d XMin1=currX.colwise().minCoeff();
        RowVector3d XMax1=currX.colwise().maxCoeff();
        RowVector3d XMin2=m2.currX.colwise().minCoeff();
        RowVector3d XMax2=m2.currX.colwise().maxCoeff();


        //New code 20/Mar/2017!!!
        double rmax1=radii.maxCoeff();
        double rmax2=m2.radii.maxCoeff();
        XMin1.array()-=rmax1;
        XMax1.array()+=rmax1;
        XMin2.array()-=rmax2;
        XMax2.array()+=rmax2;
        //end of new code

        //checking all axes for non-intersection of the dimensional interval
        for (int i=0;i<3;i++)
            if ((XMax1(i)<XMin2(i))||(XMax2(i)<XMin1(i)))
                return false;

        return true;  //all dimensional intervals are overlapping = intersection

    }


    //this function creates all collision constraints between two particle meshes
	void CreateCollisionConstraint(const Mesh& m, vector<Constraint>& collConstraints) {


		//collision between bounding boxes
		if (!isBoxCollide(m))
			return;

		if (isFixed && m.isFixed)
			return;
        //checking collision between every two particles
        //This assumes that prevX are not intersecting. That could potentially cause artifacts
        for (int i=0;i<currX.rows();i++){
            for (int j=0;j<m.currX.rows();j++){
                if ((currX.row(i)-m.currX.row(j)).norm()>radii(i)+m.radii(j))  //mutual distance longer than sum of radii
                    continue;  //no collision

                //cout<<"collision between particle at "<<currX.row(i)<<" with radius "<<radii(i)<<" and particle at "<<m.currX.row(j)<<" with radius "<<m.radii(j)<<endl;

                //naive constraint
               // if ( !isFixed || !m.isFixed )
                {
                    VectorXi particleIndices(6); particleIndices<<rawOffset+3*i, rawOffset+3*i+1, rawOffset+3*i+2, m.rawOffset+3*j, m.rawOffset+3*j+1, m.rawOffset+3*j+2;
                    VectorXd rawInvMasses(6); rawInvMasses<<invMasses(i), invMasses(i),invMasses(i),m.invMasses(j),m.invMasses(j),m.invMasses(j);
                    VectorXd rawRadii(6); rawRadii <<radii(i), radii(i), radii(i), m.radii(j), m.radii(j), m.radii(j);
                    Constraint c(COLLISION,  particleIndices, rawRadii,rawInvMasses, 0.0, 1.0, false);
                    collConstraints.push_back(c);
                }
            }
        }
    }

    //Updating the velocities currVel of the particles, from currImpulses and the external forces
    //You need to modify this to integrate from acceleration in the field (basically gravity)
    void integrateVelocity(double timeStep){

        if (isFixed)
            return;



		RowVector3d GravityEffect = RowVector3d( 0, ( -GravityAcceleration * timeStep), 0 );
		for ( int rowCounter = 0; rowCounter < currVel.rows(); rowCounter++ )
		{


			//currImpulses.row(rowCounter) += GravityEffect / invMasses(rowCounter); //alternative way with impulses instead of directly
			currVel.row(rowCounter) += GravityEffect;

		}

        for (int ImpulseCounter = 0; ImpulseCounter <  currImpulses.rows(); ImpulseCounter++)
        {
			if (!currImpulses.row(ImpulseCounter).isZero()) {
				//impulse = F * timestep
				// F = m * a
				// v+ = v- + a*timestep
				// impulse / timestep / mass * timestep + v- = v+
				// v+ = v- + impulse * invMass
				//cout << "before imp" << currVel.row(rowCounter) << endl;

				//cout << "after imp" << currVel.row(rowCounter) << endl;

				/*cout << "CRCoeff * posDiff / timeStep / invMass" << impulse;
				dV = impulse / timeStep * invMasses(ImpulseCounter) * timeStep;
				cout << "dV" << dV;*/

				//cout << ImpulseCounter << ": " << currImpulses.row(ImpulseCounter) << endl;

				currVel.row(ImpulseCounter) += currImpulses.row(ImpulseCounter) * invMasses(ImpulseCounter);

				//cout << "imp" << currImpulses.row(ImpulseCounter) << endl;
			}

        }



        currImpulses.setZero();
    }


    //Update the current position currX
    void integratePosition(double timeStep){
        if (isFixed)
            return;  //a fixed object is immobile

		//NOTE this is called after integrating velocities so it uses v(t+dt)
        /*
		for ( int rowCounter = 0; rowCounter < currVel.rows(); rowCounter++ )
		{
	//		cout << "before imp" << currX.row(rowCounter) << endl;
			currX.row(rowCounter) += currVel.row(rowCounter) * timeStep;
	//		cout << "ater imp" << currX.row(rowCounter) << endl;
	//		cout << "imp" << currVel.row(rowCounter) << endl;
		}*/
		//cout << "currXold" << currX << endl;
		//cout << "currVEL" << currVel << endl;
		currX += currVel * timeStep;
		//cout << "currX" << currX << endl;

        igl::per_vertex_normals(currX, T, currNormals);
    }


    //updating the current velocities to match the positional changes
    void projectVelocities(double timeStep){

        //NOTE this method is called after constraints are resolved to update the new velocities. lecture 7 slide 6 the yellow thingy
		/*MatrixXd currVel1 = currVel;
        for ( int rowCounter = 0; rowCounter < currVel.rows(); rowCounter++ )
        {
            currVel1(rowCounter) = (currX(rowCounter) -prevX(rowCounter) ) / timeStep;
        }*/
		currVel = (currX - prevX) / timeStep;
		//cout << "new" << currVel << endl << "old" << currVel1 << endl;

        prevX=currX;
    }


    //the full integration for the time step (velocity + position)
    void integrate(double timeStep){
        integrateVelocity(timeStep);
        integratePosition(timeStep);
    }


    Mesh(const MatrixXd& _X, const MatrixXi& _T, const int _rawOffset, const double density, const double _rigidity, const bool _isFixed, const RowVector3d& userCOM, const RowVector4d& userOrientation, const bool isTire_){
        origX=_X;
        T=_T;
        isFixed=_isFixed;
        isTire = isTire_;
        rawOffset=_rawOffset;
        rigidity=_rigidity;
        currVel=MatrixXd::Zero(origX.rows(),3);
        currImpulses=MatrixXd::Zero(origX.rows(),3);

        RowVector3d naturalCOM;  //by the geometry of the object
        Matrix3d invIT; //it is not used in this practical
        double mass; //as well

        //initializes the origianl geometry (COM + IT) of the object
        getCOMandInvIT(origX, T, density, mass, naturalCOM, invIT);

        origX.rowwise() -= naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)

        currX.resize(origX.rows(), 3);
        for (int i=0;i<currX.rows();i++)
            currX.row(i)<<QRot(origX.row(i), userOrientation)+userCOM;

        prevX=currX;

        //dynamics initialization
        VectorXd A;
        igl::doublearea(currX,T,A);
        VectorXd massV=VectorXd::Zero(currX.rows());
        for (int i=0;i<T.rows();i++){
            for (int j=0;j<3;j++)
                massV(T(i,j))+=A(i);
        }

        massV*=density/3.0;
        //massV.setOnes();
        invMasses=1.0/massV.array();

        //New ode 30/mar/2017
        if (isFixed)
            invMasses.setZero();
        //end new code

        //radii are the maximum half-edge lengths
        radii=VectorXd::Zero(currX.rows());
        for (int i=0;i<T.rows();i++){
            for (int j=0;j<3;j++){
                double edgeLength=(currX.row(T(i,j))-currX.row(T(i,(j+1)%3))).norm();
                radii(T(i,j))=std::max(radii(T(i,j)), edgeLength/2.0);
                radii(T(i,(j+1)%3))=std::max(radii(T(i,(j+1)%3)), edgeLength/2.0);
            }
        }

        igl::per_vertex_normals(currX, T, currNormals);

        MatrixXi EV;
        MatrixXi FE;
        MatrixXi EF;
        MatrixXi EFi;
        MatrixXd FESigns;

        igl::edge_topology(currX,T,EV, FE,EF);
        ConstructEFi(FE, EF, EFi, FESigns);

        for (int i=0;i<EV.rows();i++){
            int f=EF(i,0);
            int g=EF(i,1);

            //from the side i->k
            int v[4];
            v[0]=EV(i,0);
            v[1]=EV(i,1);
            v[2]=T(g,(EFi(i,1)+2)%3);
            v[3]=T(f,(EFi(i,0)+2)%3);
            for (int j=0;j<4;j++){
                for (int k=j+1;k<4 && !isFixed;k++){
                    VectorXi particleIndices(6); particleIndices << 3*v[j], 3*v[j]+1, 3*v[j]+2, 3*v[k], 3*v[k]+1, 3*v[k]+2;
                    particleIndices.array()+=rawOffset;
                    VectorXd rawRadii(6); rawRadii << radii(v[j]), radii(v[j]), radii(v[j]), radii(v[k]), radii(v[k]), radii(v[k]);
                    VectorXd rawInvMasses(6); rawInvMasses << invMasses(v[j]), invMasses(v[j]), invMasses(v[j]), invMasses(v[k]), invMasses(v[k]), invMasses(v[k]);
                    double edgeLength = (currX.row(v[j]) - currX.row(v[k])).norm();
                    meshConstraints.push_back(Constraint(RIGIDITY, particleIndices, rawRadii, rawInvMasses, edgeLength, 1.0, isTire));
                }
            }

        }

    }

    ~Mesh(){}
};



//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
    double currTime;

    VectorXd rawX;
    VectorXd rawVel;
    VectorXd rawImpulses;
    MatrixXi T;

	MatrixXi particleModelT;
	MatrixXd particleModelX;

    vector<Mesh> meshes;
    double platWidth, platHeight;  //used to create the platform constraint

    vector<Constraint> interMeshConstraints;   //constraints between meshes (mostly attachments read from the file

	vector<Spring> springs; // sad but have to do it like this for now.


    //updates from global raw indices back into mesh current positions.
    void updateMeshValues(){
        for (int i=0;i<meshes.size();i++){
            for (int j=0;j<meshes[i].currX.rows();j++){
                //cout<<"meshes[i].currX.row(j) before:"<<meshes[i].currX.row(j)<<endl;
                meshes[i].currX.row(j)<<rawX.segment(meshes[i].rawOffset+3*j,3).transpose();
                //cout<<"meshes[i].currX.row(j) after:"<<meshes[i].currX.row(j)<<endl;
                meshes[i].currVel.row(j)<<rawVel.segment(meshes[i].rawOffset+3*j,3).transpose();
                meshes[i].currImpulses.row(j)<<rawImpulses.segment(meshes[i].rawOffset+3*j,3).transpose();
			//	cout << meshes[i].currImpulses.row(j) << endl;
            }
        }
    }

    //update from mesh current positions into global raw indices.
    void updateRawValues(){
        for (int i=0;i<meshes.size();i++){
            for (int j=0;j<meshes[i].currX.rows();j++){
                rawX.segment(meshes[i].rawOffset+3*j,3)<<meshes[i].currX.row(j).transpose();
                rawVel.segment(meshes[i].rawOffset+3*j,3)<<meshes[i].currVel.row(j).transpose();
                rawImpulses.segment(meshes[i].rawOffset+3*j,3)<<meshes[i].currImpulses.row(j).transpose();
            }
        }

    }


    //adding an objects. You do not need to update this generally
    void addMesh(const MatrixXd& meshX, const MatrixXi& meshT, const double density, const double rigidity, const bool isFixed, const RowVector3d& COM, const RowVector4d orientation, const bool isTire){

        Mesh m(meshX,meshT, rawX.size(), density, rigidity, isFixed, COM, orientation, isTire);
        meshes.push_back(m);
        int oldTsize=T.rows();
        T.conservativeResize(T.rows()+meshT.rows(),3);
        T.block(oldTsize,0,meshT.rows(),3)=meshT.array()+rawX.size()/3;  //to offset T to global index
        rawX.conservativeResize(rawX.size()+meshX.size());
        rawVel.conservativeResize(rawVel.size()+meshX.size());
        rawImpulses.conservativeResize(rawImpulses.size()+meshX.size());
        updateRawValues();

        //cout<<"rawVel: "<<rawVel<<endl;
    }


    /*********************************************************************
    This function handles a single position-based time step
     1. Integrating velocities, and position
     2. detecting collisions and generating constraints
     3. Resolving constraints iteratively until the system is valid
     4. updating velocities to match positions
     *********************************************************************/
    void updateScene(double timeStep, double CRCoeff, MatrixXd& fullX, MatrixXi& fullT, const double tolerance, const int maxIterations, const MatrixXd& platX, const MatrixXi& platT){

        //1. integrating velocity, position and orientation from forces and previous states
        for (int i=0;i<meshes.size();i++)
            meshes[i].integrate(timeStep);

		// first get the curVel to rawVel
		updateRawValues();

		// damp spring velocity from previous stuff before frame
		/*if (timeStep > 0) {
			for (Spring s : springs) {
				double impulse = s.dampSpringImpulse(rawVel, timeStep); //impulse but needs to be fixed before the frame is done so need to update velocity according to this impulse right now
				double impulse2 = s.dampSpringForce(rawVel[s.getParticleIndice2()], rawVel[s.getParticleIndice1()]) * timeStep;

				//double F = impulse / timeStep;
				//RowVector3d a = F * invMasses(ImpulseCounter);
				//RowVector3d dV = a * timeStep;
				if (impulse > 0) {
				//	cout << "old1" << rawVel[s.getParticleIndice1()] << endl;
				//	cout << "old2" << rawVel[s.getParticleIndice2()] << endl;
					rawVel[s.getParticleIndice1()] += impulse / timeStep * s.invMass1 * timeStep; // can remove timeStep
					rawVel[s.getParticleIndice2()] += impulse2 / timeStep * s.invMass2 * timeStep;
				//	cout << "new1" << rawVel[s.getParticleIndice1()] << endl;
				//	cout << "new2" << rawVel[s.getParticleIndice2()] << endl;
				//	cout << "particles:" << s.getParticleIndice1() << endl << s.getParticleIndice2() << endl;
				}
			}
		}*/
		// move rawVel back to curVel
		/*for (int i = 0; i < meshes.size(); i++) {
			//cout << "old currVel" << endl;
			//cout << meshes[i].currVel << endl;
		}
		updateMeshValues(); // need to update the meshValues so the dampening doesn't get removed
		for (int i = 0; i < meshes.size(); i++) {
			//cout << "new currVel" << endl;
			//cout << meshes[i].currVel << endl;
		}
		// also need to update the position to wait first uhm idk
		// cant use intergrate since that gravity gets added twice so
		for (int i = 0; i < meshes.size(); i++) {
			meshes[i].integratePosition(timeStep);
		}*/

        //cout<<"raw positions: "<<rawV<<endl;

        vector<Constraint> fullConstraints;
        vector<Constraint> collConstraints;
        vector<Constraint> rigidityConstraints;
        vector<Constraint> barrierConstraints;

        //2. detecting collisions and generating constraints the are aggragated in collConstraints
        for (int i=0;i<meshes.size();i++)
			for (int j = i + 1; j < meshes.size(); j++) {
				meshes[i].CreateCollisionConstraint(meshes[j], collConstraints);
			}

        //NOTE creating Platform Barrier constraints
        for (int MeshCounter = 0; MeshCounter < meshes.size(); MeshCounter++)
        {
            Mesh Object = meshes[MeshCounter];

            //NOTE RIGIDITY, inter-mesh constraints
            if ( !Object.isFixed )
            {
                rigidityConstraints.insert( rigidityConstraints.end(), Object.meshConstraints.begin(), Object.meshConstraints.end() );
                //NOTE BARRIER, loop through the objects particles
                for ( int ParticleCounter = 0; ParticleCounter < Object.currX.rows() ; ParticleCounter++ )
                {
    				RowVector3d ParticlePosition = Object.currX.row(ParticleCounter);
    				double ParticleRadius = Object.radii(ParticleCounter);

                    //HACK assuming the first object is the tire
                    /*
                    if ( MeshCounter == 0 )
                    {
                        TireCenterX = ( TireCenterX + ParticlePosition(0) ) / 2;
                        TireCenterY = ( TireCenterX + ParticlePosition(1) ) / 2;
                        TireCenterZ = ( TireCenterX + ParticlePosition(2) ) / 2;
                    }
                    */
                    VectorXi particleIndices(1); particleIndices << Object.rawOffset + (ParticleCounter * 3);
    				VectorXd rawInvMasses(1); rawInvMasses << Object.invMasses(ParticleCounter);
                    VectorXd rawRadii(1); rawRadii << ParticleRadius;

                    bool XLimitation = abs(ParticlePosition[0]) - ParticleRadius < (platWidth / 2);
                    bool ZLimitation = abs(ParticlePosition[2]) - ParticleRadius < (platWidth / 2);

                    //NOTE check whether object need the barrier constraint or not
                    if ( XLimitation && ZLimitation )
                    {
                        Constraint platformBarrier(BARRIER, particleIndices, rawRadii, rawInvMasses, platHeight/2, 1.0, false);
						barrierConstraints.push_back(platformBarrier);
                    }
                }
            }
        }

        //NOTE lecture 7 slide 11 onwards

        //aggregating mesh and inter-mesh constraints

		fullConstraints.insert(fullConstraints.end(), collConstraints.begin(), collConstraints.end());
        fullConstraints.insert(fullConstraints.end(), rigidityConstraints.begin(), rigidityConstraints.end());
		fullConstraints.insert(fullConstraints.end(), interMeshConstraints.begin(), interMeshConstraints.end());
        fullConstraints.insert(fullConstraints.end(), barrierConstraints.begin(), barrierConstraints.end());

        //3. Resolving constraints iteratively until the system is valid (all constraints are below "tolerance" , or passed maxIteration*fullConstraints.size() iterations
        //add proper impulses to rawImpulses for the corrections (CRCoeff*posDiff/timeStep). Don't do that on the initialization step.
		bool done = false;
		for (int iteration = 0; !done && iteration < maxIterations; iteration++)
        {
			done = true;
            double Impulse = -1;
			for (Constraint c : fullConstraints) {
				VectorXd posDiffs;

                if ( c.constraintType == BARRIER )
                {
					double impulse = -1;
					for ( int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++ )
                    {
						int indices = (c.particleIndices[ParticleIndex]) + 1;
						VectorXd test(1); test << rawX[indices];

						c.resolveConstraint(test, posDiffs);

						// norm should always be positive
						if ( posDiffs.norm() > tolerance )
                        {

							test += posDiffs;
                            done = false;
                            rawX[(c.particleIndices[ParticleIndex]) + 1] = test(0);
                            if ( timeStep > 0.0 && iteration > 0)
                            {
                             //   double Radius = ( posDiffs(0) > 0 ) ? ( c.radii(0) ) : ( -c.radii(0) ) ;
                              //  rawImpulses[(c.particleIndices[ParticleIndex]) + 1] += ( ( CRCoeff * ( posDiffs(0) * ( 1 +  (2 / Radius) ) ) ) / timeStep );
							    double invMass = c.invMassMatrix.row(ParticleIndex)[ParticleIndex];
                                rawImpulses[(c.particleIndices[ParticleIndex]) + 1] += ( (  CRCoeff  * posDiffs(0) ) / timeStep ) / invMass;

								/*double a = CRCoeff*posDiffs(ParticleIndex) / timeStep;
								rawImpulses(indices) = a * 1;
								if (impulse == -1 || impulse == 0) {
									impulse = a;
								}

								if (impulse > 0) {

									int tempOffset = 0;
									int tempOffset2 = 0;
									// find the mesh where this particle belongs too
									for (int i = 0; i < meshes.size(); i++) {
										tempOffset2 = meshes[i].rawOffset;
										if (tempOffset2 <= indices) {
											if (tempOffset2 > indices) {
												break;
											}
											else {
												tempOffset = i;
											}
										}

									}
									// update all particles of this mesh
									for (int i = 0; i < meshes[tempOffset].currX.rows(); i++) {
										rawImpulses(meshes[tempOffset].rawOffset + i * 3 + 1) += impulse;
									}
								}*/
                            }



						}
                    }
				}

				if ( c.constraintType == ATTACHMENT)
                {
					VectorXd AllParticles(6);
					AllParticles << rawX[(c.particleIndices[0])], rawX[(c.particleIndices[1])], rawX[(c.particleIndices[2])],
						rawX[(c.particleIndices[3])], rawX[(c.particleIndices[4])], rawX[(c.particleIndices[5])];
					c.updateValueGradient(AllParticles);
					if (abs(c.currValue) > tolerance + (RigidityAllowance / c.refValue))
					{
						done = false;
						c.resolveConstraint(AllParticles, posDiffs);
						for (int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++)
						{
							rawX[(c.particleIndices[ParticleIndex])] += posDiffs(ParticleIndex);
							if (timeStep > 0.0 && iteration > 0)
							{
								double invMass = c.invMassMatrix.row(ParticleIndex)[ParticleIndex];
								rawImpulses[(c.particleIndices[ParticleIndex])] += (CRCoeff * posDiffs(ParticleIndex) / timeStep) / invMass;
							}
						}
					}
				}

				if  ( c.constraintType == RIGIDITY )
                {
                    VectorXd AllParticles(6);
                    AllParticles << rawX[ ( c.particleIndices[0] )], rawX[ ( c.particleIndices[1] ) ], rawX[ ( c.particleIndices[2] ) ],
                                    rawX[ ( c.particleIndices[3] )], rawX[ ( c.particleIndices[4] ) ], rawX[ ( c.particleIndices[5] ) ];
                    c.updateValueGradient( AllParticles );
                    double Range = c.isTire ? (c.refValue * ( 1 - (TirePressure / MaxTirePressure) )) : (0);
                    if ( c.currValue > 0 || c.currValue < - tolerance - Range )
                    {
                        done = false;
                        c.resolveConstraint( AllParticles, posDiffs );
                        for (int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++)
                        {
                            rawX[(c.particleIndices[ParticleIndex])] += posDiffs(ParticleIndex);
                            if ( timeStep > 0.0 && iteration > 0 )
                            {
								double invMass = c.invMassMatrix.row(ParticleIndex)[ParticleIndex];
                                rawImpulses[(c.particleIndices[ParticleIndex])] += ( CRCoeff * posDiffs(ParticleIndex) / timeStep ) / invMass;
                            }
                        }
                    }
                }

				if (c.constraintType == ATTACHMENTSTATIC)
				{
					VectorXd CurrentParticlePositions(c.particleIndices.size());
					for (int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++) {
						CurrentParticlePositions(ParticleIndex) = rawX[(c.particleIndices[ParticleIndex])];
					}
					c.updateValueGradient(CurrentParticlePositions);//cant do else since c isnt used next time :(
					if (abs(c.currValue) > tolerance) // needs to happen here since resolve doesnt have the tolerance otherwise update could be removed
					{
						done = false;
						c.resolveConstraint(CurrentParticlePositions, posDiffs);
						for (int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++)
						{
							rawX[(c.particleIndices[ParticleIndex])] += posDiffs(ParticleIndex);
							if (timeStep > 0.0 && (ParticleIndex > 2) && posDiffs(ParticleIndex) > tolerance)
							{
								rawImpulses[(c.particleIndices[ParticleIndex])] += ((CRCoeff * posDiffs(ParticleIndex)) / timeStep);
							}
						}

					}
				}

                /*if ( c.constraintType == COLLISION )
                {
                    /*
                    VectorXd AllParticles(6);
                    AllParticles << rawX[ ( c.particleIndices[0] )], rawX[ ( c.particleIndices[1] ) ], rawX[ ( c.particleIndices[2] ) ],
                                    rawX[ ( c.particleIndices[3] )], rawX[ ( c.particleIndices[4] ) ], rawX[ ( c.particleIndices[5] ) ];
                    c.updateValueGradient( AllParticles );
                    if ( abs(c.currValue) > tolerance)
                    {
                        done = false;
                        c.resolveConstraint( AllParticles, posDiffs );
                        for (int ParticleIndex = 0; ParticleIndex < c.particleIndices.size(); ParticleIndex++)
                        {
                            rawX[(c.particleIndices[ParticleIndex])] += posDiffs(ParticleIndex);
                            if ( timeStep > 0.0 && iteration > 0)
                            {
								double invMass = c.invMassMatrix.row(ParticleIndex)[ParticleIndex];
                                rawImpulses[(c.particleIndices[ParticleIndex])] += ( ( CRCoeff ) * posDiffs(ParticleIndex) ) / timeStep  / invMass;
                            }
                        }
                    }
                    */
                }

			for (Spring s : springs) {
				double impulse = s.getImpulse(rawX, timeStep);
				rawImpulses[s.getParticleIndice1()] += impulse;
				rawImpulses[s.getParticleIndice2()] += -impulse;	
			}
		}

        fullConstraints.clear();
        collConstraints.clear();
        barrierConstraints.clear();
        rigidityConstraints.clear();

        //4. updating velocities to match positions (the position-based step)
        updateMeshValues();
        if (timeStep>tolerance)  //just to allow initialization with t=0.0 in the beginning of the run.
            for (int i=0;i<meshes.size();i++)
                meshes[i].projectVelocities(timeStep);


        updateRawValues();

        //Updating visualization variables
        currTime+=timeStep;
        fullX.conservativeResize(rawX.size()/3+platX.rows(),3);
        for (int i=0;i<rawX.size()/3;i++)
            fullX.row(i)=rawX.segment(3*i,3).transpose();

        fullX.block(rawX.size()/3, 0, platX.rows(), 3)=platX;

        fullT.conservativeResize(T.rows()+platT.rows(),3);
        fullT<<T, platT.array()+ rawX.size() / 3;
    }

    //loading a scene from the scene .txt files
    //you do not need to update this function
    bool loadScene(const std::string dataFolder, const std::string sceneFileName, const double _platWidth, const double _platHeight, VectorXi& attachM1, VectorXi& attachV1, VectorXi& attachM2, VectorXi& attachV2, VectorXi& springM1, VectorXi& springV1, VectorXi& springM2, VectorXi& springV2){

        platWidth=_platWidth;
        platHeight=_platHeight;
        ifstream sceneFileHandle;
        sceneFileHandle.open(dataFolder+std::string("/")+sceneFileName);
        if (!sceneFileHandle.is_open())
            return false;
		int numofObjects, numofConstraints, numOfSprings;
        cout << '1' << endl;
        currTime=0;
        sceneFileHandle>>numofObjects>>numofConstraints>>numOfSprings;
		RowVector3d COMMesh1;
		RowVector4d orientationMesh1;
        for (int i=0;i<numofObjects;i++){
            MatrixXi objT;
            MatrixXd objX;
            std::string OFFFileName;
            bool isFixed;
            double density, rigidity;
            RowVector3d COM;
            RowVector4d orientation;
            sceneFileHandle>>OFFFileName>>density>>rigidity>>isFixed>>COM(0)>>COM(1)>>COM(2)>>orientation(0)>>orientation(1)>>orientation(2)>>orientation(3);
            orientation.normalize();
            igl::readOFF(dataFolder+std::string("/")+OFFFileName,objX,objT);
			bool isTire = false;
			if (i == 0)
				isTire = true;
            addMesh(objX,objT,density, rigidity, isFixed, COM, orientation, isTire);
        }

        //reading and adding inter-mesh attachment constraints
        attachM1.resize(numofConstraints);
        attachV1.resize(numofConstraints);
        attachM2.resize(numofConstraints); 
        attachV2.resize(numofConstraints);
		for (int i = 0; i < numofConstraints; i++) {
			sceneFileHandle >> attachM1(i) >> attachV1(i) >> attachM2(i) >> attachV2(i);

			int rawIndice1 = meshes[attachM1(i)].rawOffset + 3 * attachV1(i);
			int rawIndice2 = meshes[attachM2(i)].rawOffset + 3 * attachV2(i);
			VectorXi particleIndices(6); particleIndices << rawIndice1, rawIndice1 + 1, rawIndice1 + 2, rawIndice2, rawIndice2 + 1, rawIndice2 + 2;
			double radii1 = meshes[attachM1(i)].radii(attachV1(i));
			double radii2 = meshes[attachM2(i)].radii(attachV2(i));
			VectorXd rawRadii(6); rawRadii << radii1, radii1, radii1, radii2, radii2, radii2;
			double invMass1 = meshes[attachM1(i)].invMasses(attachV1(i));
			double invMass2 = meshes[attachM2(i)].invMasses(attachV2(i));
			RowVector3d pos1 = RowVector3d(rawX[rawIndice1], rawX[rawIndice1 + 1], rawX[rawIndice1 + 2]);
			RowVector3d pos2 = RowVector3d(rawX[rawIndice2], rawX[rawIndice2 + 1], rawX[rawIndice2 + 2]);
			VectorXd rawInvMasses(6); rawInvMasses << invMass1, invMass1, invMass1, invMass2, invMass2, invMass2;

			double edgeLength = (pos1 - pos2).norm();

			cout << edgeLength << endl;

			interMeshConstraints.push_back(Constraint(ATTACHMENT, particleIndices, rawRadii, rawInvMasses, edgeLength, 1.0, false));
		}
		springM1.resize(numOfSprings);
		springV1.resize(numOfSprings);
		springM2.resize(numOfSprings);
		springV2.resize(numOfSprings);
		for (int i = 0; i<numOfSprings; i++) {
			
			sceneFileHandle >> springM1(i) >> springV1(i) >> springM2(i) >> springV2(i);

			int rawIndice1 = meshes[springM1(i)].rawOffset + 3 * springV1(i);
			int rawIndice2 = meshes[springM2(i)].rawOffset + 3 * springV2(i);
			VectorXi particleIndices(6); particleIndices << rawIndice1, rawIndice1 + 1, rawIndice1 + 2, rawIndice2, rawIndice2 + 1, rawIndice2 + 2;
			double radii1 = meshes[springM1(i)].radii(springV1(i));
			double radii2 = meshes[springM2(i)].radii(springV2(i));
			VectorXd rawRadii(6); rawRadii << radii1, radii1, radii1, radii2, radii2, radii2;
			double invMass1 = meshes[springM1(i)].invMasses(springV1(i));
			double invMass2 = meshes[springM2(i)].invMasses(springV2(i));
			RowVector3d pos1 = RowVector3d(rawX[rawIndice1], rawX[rawIndice1 + 1], rawX[rawIndice1 + 2]);
			RowVector3d pos2 = RowVector3d(rawX[rawIndice2], rawX[rawIndice2 + 1], rawX[rawIndice2 + 2]);
			VectorXd rawInvMasses(6); rawInvMasses << invMass1, invMass1, invMass1, invMass2, invMass2, invMass2;

			double edgeLength = (pos1 - pos2).norm();

			//interMeshConstraints.push_back(Constraint(ATTACHMENT, particleIndices, rawRadii, rawInvMasses, edgeLength, 1.0, false));

			double K = 50000;
			double dampingCoeffecient = 25000;
			// + 1 since y axis
			springs.push_back(Spring(rawIndice1+1, rawIndice2+1, invMass1, invMass2, edgeLength, K, dampingCoeffecient));
		}
        return true;
    }



    Scene(){}
    ~Scene(){}
};




#endif

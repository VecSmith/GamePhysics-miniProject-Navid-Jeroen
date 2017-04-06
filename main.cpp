#include <igl/viewer/Viewer.h>
#include "scene.h"

Eigen::MatrixXd X;
Eigen::MatrixXi T;

double currTime = 0;

//initial values
double timeStep = 0.02;
double CRCoeff= 1.0;
double tolerance=10e-3;
int maxConstraintIterations=100;

MatrixXd platX;
MatrixXi platT;
VectorXi attachM1, attachM2, attachV1, attachV2;
double platWidth=100.0;
double platHeight=10.0;
double RigidityAllowance = 0.0;

Scene scene;

void createPlatform(const double platWidth, const double platHeight)
{
    platX.resize(8,3);
    platT.resize(12,3);
    platX<<-platWidth/2.0,platHeight/2.0,-platWidth/2.0,
            -platWidth/2.0,platHeight/2.0,platWidth/2.0,
            platWidth/2.0,platHeight/2.0,platWidth/2.0,
            platWidth/2.0,platHeight/2.0, -platWidth/2.0,
            -platWidth/2.0,-platHeight/2.0,-platWidth/2.0,
            -platWidth/2.0,-platHeight/2.0,platWidth/2.0,
            platWidth/2.0,-platHeight/2.0,platWidth/2.0,
            platWidth/2.0,-platHeight/2.0, -platWidth/2.0;
    platT<<0,1,2,
            2,3,0,
            6,5,4,
            4,7,6,
            1,0,5,
            0,4,5,
            2,1,6,
            1,5,6,
            3,2,7,
            2,6,7,
            0,3,4,
            3,7,4;
}

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == ' ')
    {
        viewer.core.is_animating = !viewer.core.is_animating;
        return true;
    }

    if (key == 'S')
    {
        if (!viewer.core.is_animating){
            scene.updateScene(timeStep, CRCoeff, X,T, tolerance,maxConstraintIterations, platX, platT);
            currTime+=timeStep;
            std::cout <<"currTime: "<<currTime<<std::endl;
            return true;
        }
    }
    return false;
}


bool pre_draw(igl::viewer::Viewer &viewer)
{
    using namespace Eigen;
    using namespace std;

    if (viewer.core.is_animating){
        scene.updateScene(timeStep, CRCoeff, X,T, tolerance,maxConstraintIterations, platX, platT);
        currTime+=timeStep;
        //cout <<"currTime: "<<currTime<<endl;
    }
    viewer.data.clear();
    viewer.data.set_mesh(X,T);
    MatrixXd P1, P2;
    P1.resize(attachM1.rows(),3);
    P2.resize(attachM2.rows(),3);
    for (int i=0;i<attachM1.rows();i++){
        P1.row(i)<<scene.meshes[attachM1(i)].currX.row(attachV1(i));
        P2.row(i)<<scene.meshes[attachM2(i)].currX.row(attachV2(i));
    }
    viewer.data.add_edges(P1,P2,Eigen::RowVector3d(255.0,0.0,0.0));

    return false;
}



int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;


    // Load scene
    if (argc<2){
        cout<<"Please provide name of scene file!"<<endl;
		argv[1] = "singleboxscene.txt";
        //return 0;
    }
    cout<<"scene file: "<<std::string(argv[1])<<endl;


    //create platform
    createPlatform(100.0, 10.0);

    scene.loadScene(std::string(DATA_PATH), std::string(argv[1]),platWidth, platHeight, attachM1, attachV1, attachM2, attachV2);

    //scene.addMesh(platV, platF, 10000.0, true, platCOM, platOrientation);*/
    scene.updateScene(0.0, CRCoeff, X,T, tolerance,maxConstraintIterations, platX, platT);


    cout<<"initial X: "<<X<<endl;
    cout<<"initial T: "<<T<<endl;

    // Viewer Settings
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(X,T);
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 50.;


    //Adding options to GUI
    //To add new options, just modify the code below in a similar manner.
    viewer.callback_init = [&](igl::viewer::Viewer& viewer)
    {
        //algorithmic options
        viewer.ngui->addGroup("Simulation");
        viewer.ngui->addVariable("CR Coeff",CRCoeff);
        viewer.ngui->addVariable<double>("Time Step",[&](double val) {
            viewer.core.animation_max_fps = (((int)1.0/val));
            timeStep=val;
        },[&]() {
            return timeStep;
        });

        viewer.ngui->addVariable<double>("Rigidity Range",[&](double val_0) {
            RigidityAllowance=val_0;
        },[&]() {
            return RigidityAllowance;
        });

        // call to generate menu
        viewer.screen->performLayout();
        return false;
    };


    cout<<"Press [space] to toggle continuous simulation" << endl;
    cout<<"Press 'S' to advance time step-by-step"<<endl;
    viewer.launch();
}

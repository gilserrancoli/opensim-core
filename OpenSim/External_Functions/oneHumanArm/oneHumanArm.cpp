/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Gil Serrancoli
    Contributor: Joris Gillis, Antoine Falisse, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Simulation/SimbodyEngine/Body.h>
#include "SimTKcommon/internal/recorder.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace OpenSim;

/*  The function F describes the OpenSim model and, implicitly, the skeleton
    dynamics. F takes as inputs joint positions and velocities (states x),
    joint accelerations (controls u), a platform perturbation value (p), and
    returns the joint torques. F is templatized using type T. F(x,u,p)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 10;     // # degrees of freedom
constexpr int NX = ndof*2;  // # states
constexpr int NU = ndof;    // # controls
constexpr int NR = ndof;    // # residual torques

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

	// OpenSim model: create components
    /// Model
	Model* model;
    /// Bodies
	OpenSim::Body* pelvis;
	OpenSim::Body* torso;
	OpenSim::Body* humerus_r;
	OpenSim::Body* ulna_r;
	OpenSim::Body* radius_r;
	OpenSim::Body* hand_r;
	
    /// Joints
	OpenSim::CustomJoint* ground_pelvis;
	OpenSim::WeldJoint* back;
	OpenSim::CustomJoint* acromial_r;
	OpenSim::CustomJoint* elbow_r;
	OpenSim::CustomJoint* radioulnar_r;
	OpenSim::CustomJoint* radius_hand_r;



	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	pelvis = new OpenSim::Body("pelvis", 11.777, Vec3(-0.0707, 0, 0),
            Inertia(0.1028, 0.0871, 0.0579, 0, 0, 0));
	torso = new OpenSim::Body("torso", 26.8266, Vec3(-0.03, 0.32, 0),
		Inertia(1.4745, 0.7555, 1.4314, 0, 0, 0));
	humerus_r = new OpenSim::Body("humerus_r", 2.0325, Vec3(0, -0.164502, 0),
            Inertia(0.011946, 0.004121, 0.013409, 0, 0, 0));
	ulna_r = new OpenSim::Body("ulna_r", 0.6075, Vec3(0, -0.120525, 0),
            Inertia(0.002962, 0.000618, 0.003213, 0, 0, 0));
	radius_r = new OpenSim::Body("radius_r", 0.6075, Vec3(0, -0.120525, 0),
		Inertia(0.002962, 0.000618, 0.003213, 0, 0, 0));
	hand_r = new OpenSim::Body("hand_r", 0.4575, Vec3(0, -0.068095, 0),
		Inertia(0.000892, 0.000547, 0.00134, 0, 0, 0));
	

    /// Joint specifications
	SpatialTransform st_ground_pelvis;
	st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
	st_ground_pelvis[0].setFunction(new LinearFunction());
	st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
	st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
	st_ground_pelvis[1].setFunction(new LinearFunction());
	st_ground_pelvis[1].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_rotation", 1, 1));
	st_ground_pelvis[2].setFunction(new LinearFunction());
	st_ground_pelvis[2].setAxis(Vec3(0, 1, 0));
	ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0),
		    *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
	back = new WeldJoint("back", *pelvis, Vec3(-0.1007, 0.0815, 0),
		Vec3(0), *torso, Vec3(0), Vec3(0));
	
	SpatialTransform st_acromial_r;
	st_acromial_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
	st_acromial_r[0].setFunction(new LinearFunction());
	st_acromial_r[0].setAxis(Vec3(0, 0, 1));
	st_acromial_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
	st_acromial_r[1].setFunction(new LinearFunction());
	st_acromial_r[1].setAxis(Vec3(1, 0, 0));
	st_acromial_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
	st_acromial_r[2].setFunction(new LinearFunction());
	st_acromial_r[2].setAxis(Vec3(0, 1, 0));
	acromial_r = new CustomJoint("acromial_r", *torso, Vec3(0.003155, 0.3715, 0.17),
            Vec3(0), *humerus_r, Vec3(0), Vec3(0),st_acromial_r);
	
	SpatialTransform st_elbow_r;
	st_elbow_r[2].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
	st_elbow_r[2].setFunction(new LinearFunction());
	st_elbow_r[2].setAxis(Vec3(0.226047, 0.022269, 0.973862));
	elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.013144, -0.286273, -0.009595), Vec3(0),
            *ulna_r, Vec3(0), Vec3(0), st_elbow_r);
	
	SpatialTransform st_radioulnar_r;
	st_radioulnar_r[2].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
	st_radioulnar_r[2].setFunction(new LinearFunction());
	st_radioulnar_r[2].setAxis(Vec3(0.056398, 0.998406, 0.001952));
	radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.006727, -0.013007, 0.026083),
		Vec3(0), *radius_r, Vec3(0), Vec3(0), st_radioulnar_r);
	
	SpatialTransform st_radius_hand_r;
	st_radius_hand_r[0].setCoordinateNames(OpenSim::Array<std::string>("wrist_flex_r", 1, 1));
	st_radius_hand_r[0].setFunction(new LinearFunction());
	st_radius_hand_r[0].setAxis(Vec3(0, 0, 1));
	st_radius_hand_r[1].setCoordinateNames(OpenSim::Array<std::string>("wrist_dev_r", 1, 1));
	st_radius_hand_r[1].setFunction(new LinearFunction());
	st_radius_hand_r[1].setAxis(Vec3(1, 0, 0));
	st_radius_hand_r[2].setAxis(Vec3(0, 1, 0));
	radius_hand_r = new CustomJoint("radius_hand_r", *radius_r, Vec3(-0.008797, -0.235841, 0.01361), Vec3(0), *hand_r, Vec3(0), Vec3(0), st_radius_hand_r);
	
    /// Add bodies and joints to model
	model->addBody(pelvis);
	model->addBody(torso);
	model->addBody(humerus_r);
	model->addBody(ulna_r);
	model->addBody(radius_r);
	model->addBody(hand_r);
	model->addJoint(ground_pelvis);
	model->addJoint(back);
	model->addJoint(acromial_r);
	model->addJoint(elbow_r);
	model->addJoint(radioulnar_r);
	model->addJoint(radius_hand_r);

	// Initialize system and state
	State* state;
	state = new State(model->initSystem());

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
	T ua[NU]; /// joint accelerations (Qdotdots) - controls

	// Assign inputs to model variables
    /// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
	for (int i = 0; i < NU; ++i) ua[i] = u[i];

    // Set state variables and realize
	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// Compute residual forces
    /// appliedMobilityForces (# mobilities)
	Vector appliedMobilityForces(ndof);
	appliedMobilityForces.setToZero();
    /// appliedBodyForces (# bodies + ground)
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1;
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Set gravity
	Vec3 gravity(0);
	gravity[1] = -9.81;
	/// Set platform perturbation
	Vec3 Vec_gravity = Vec3(0, gravity[1], 0);
	/// Add to model
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
			Vec_gravity*model->getBodySet().get(i).getMass(),
            appliedBodyForces);
	}
    /// knownUdot
	Vector knownUdot(ndof);
	knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
	Vector residualMobilityForces(ndof);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

	// Extract results
    /// Residual forces
	for (int i = 0; i < NR; ++i) res[0][i] =
        value<T>(residualMobilityForces[i]);

	return 0;

}

/* In main(), the Recorder is used to save the expression graph of function F.
This expression graph is saved as a MATLAB function named foo.m. From this
function, a c-code can be generated via CasADi and then compiled as a dll. This
dll is then imported in MATLAB as an external function. With this workflow,
CasADi can use algorithmic differentiation to differentiate the function F.
*/
int main() {

	Recorder x[NX];
	Recorder u[NU];
	Recorder tau[NR];

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

	Recorder::stop_recording();

	return 0;

}

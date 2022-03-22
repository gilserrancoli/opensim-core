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
constexpr int ndof = 7;     // # degrees of freedom
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
	OpenSim::Body* femur_r;
	OpenSim::Body* tibia_r;
	OpenSim::Body* talus_r;
	OpenSim::Body* calcn_r;
	OpenSim::Body* toes_r;
	OpenSim::Body* torso;
    /// Joints
	OpenSim::CustomJoint* ground_pelvis;
	OpenSim::PinJoint* hip_r;
	OpenSim::CustomJoint* knee_r;
	OpenSim::PinJoint* ankle_r;
	OpenSim::WeldJoint* subtalar_r;
	OpenSim::WeldJoint* mtp_r;
	OpenSim::PinJoint* back;



	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	pelvis = new OpenSim::Body("pelvis", 11.777, Vec3(-0.0707, 0, 0),
            Inertia(0.1028, 0.0871, 0.0579, 0, 0, 0));
	femur_r = new OpenSim::Body("femur_r", 9.3014, Vec3(0, -0.17, 0),
            Inertia(0.1339, 0.0351, 0.1412, 0, 0, 0));
	tibia_r = new OpenSim::Body("tibia_r", 3.7075, Vec3(0, -0.1867, 0),
            Inertia(0.0504, 0.0051, 0.0511, 0, 0, 0));
	talus_r = new OpenSim::Body("talur_r", 0.1, Vec3(0),
		Inertia(0.001));
	calcn_r = new OpenSim::Body("calcn_r", 1.25, Vec3(0.1, 0.03, 0),
		Inertia(0.0014, 0.0039, 0.0041, 0, 0, 0));
	toes_r = new OpenSim::Body("toes_r", 0.2166, Vec3(0.0346, 0.006, -0.0175),
		Inertia(0.0001, 0.0002, 0.0001, 0, 0, 0));
	torso = new OpenSim::Body("torso", 34.2366, Vec3(-0.03, 0.32, 0),
		Inertia(1.4745, 0.7555, 1.4314, 0, 0, 0));

    /// Joint specifications
	SpatialTransform st_ground_pelvis;
	st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
	st_ground_pelvis[2].setFunction(new LinearFunction());
	st_ground_pelvis[2].setAxis(Vec3(0, 0, 1));
	st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
	st_ground_pelvis[3].setFunction(new LinearFunction());
	st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
	st_ground_pelvis[4].setFunction(new LinearFunction());
	st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));


	ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0),
		    *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
	hip_r = new PinJoint("hip_r", *pelvis, Vec3(-0.0707, -0.0661, 0.0835),
            Vec3(0), *femur_r, Vec3(0), Vec3(0));
	SpatialTransform st_knee_r;
	st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[2].setFunction(new LinearFunction());
	st_knee_r[2].setAxis(Vec3(0, 0, 1));
	//st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("", 1, 1));
	st_knee_r[4].setFunction(new Constant(-0.41));
	st_knee_r[4].setAxis(Vec3(0, 1, 0));

	knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0, 0, 0), Vec3(0),
            *tibia_r, Vec3(0), Vec3(0), st_knee_r);
	ankle_r = new PinJoint("ankle_r", *tibia_r, Vec3(0, -0.43, 0),
		Vec3(0), *talus_r, Vec3(0), Vec3(0));
	subtalar_r = new WeldJoint("subtalar_r", *talus_r, Vec3(-0.04877, -0.04195, 0.00792), Vec3(0), *calcn_r, Vec3(0), Vec3(0));
	mtp_r = new WeldJoint("mtp_r", *calcn_r, Vec3(0.1788, -0.002, 0.00108), Vec3(0), *toes_r, Vec3(0), Vec3(0));
	back = new PinJoint("back", *pelvis, Vec3(-0.1007, 0.0815, 0),
		Vec3(0), *torso, Vec3(0), Vec3(0));
    /// Add bodies and joints to model
	model->addBody(pelvis);
	model->addBody(femur_r);
	model->addBody(tibia_r);
	model->addBody(talus_r);
	model->addBody(calcn_r);
	model->addBody(toes_r);
	model->addBody(torso);
	model->addJoint(ground_pelvis);
	model->addJoint(hip_r);
	model->addJoint(knee_r);
	model->addJoint(ankle_r);
	model->addJoint(subtalar_r);
	model->addJoint(mtp_r);
	model->addJoint(back);

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

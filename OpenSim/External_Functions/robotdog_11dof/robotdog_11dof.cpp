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
constexpr int n_in = 3;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 11;     // # degrees of freedom
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
	OpenSim::Body* torso;
	OpenSim::Body* humer1;
	OpenSim::Body* ulna1;
	OpenSim::Body* humer2;
	OpenSim::Body* ulna2;
	OpenSim::Body* femur3;
	OpenSim::Body* tibia3;
	OpenSim::Body* femur4;
	OpenSim::Body* tibia4;
	
    /// Joints
	OpenSim::CustomJoint* torso_ground;
	OpenSim::PinJoint* shoulder1;
	OpenSim::PinJoint* elbow1;
	OpenSim::PinJoint* shoulder2;
	OpenSim::PinJoint* elbow2;
	OpenSim::PinJoint* hip3;
	OpenSim::PinJoint* knee3;
	OpenSim::PinJoint* hip4;
	OpenSim::PinJoint* knee4;

	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	torso = new OpenSim::Body("torso", 20.0, Vec3(0),
            Inertia(6.87, 0.1, 6.87, 0, 0, 0));
	humer1 = new OpenSim::Body("humer1", 8.0, Vec3(0),
			Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	ulna1 = new OpenSim::Body("ulna1", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	humer2 = new OpenSim::Body("humer2", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	ulna2 = new OpenSim::Body("ulna2", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	femur3 = new OpenSim::Body("femur3", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	tibia3 = new OpenSim::Body("tibia3", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	femur4 = new OpenSim::Body("femur4", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	tibia4 = new OpenSim::Body("tibia4", 8.0, Vec3(0),
		Inertia(1.04, 0.04, 1.04, 0, 0, 0));
	
    /// Joint specifications
	SpatialTransform st_pin1;
	st_pin1[0].setCoordinateNames(OpenSim::Array<std::string>("r1_z", 1, 1));
	st_pin1[0].setFunction(new LinearFunction());
	st_pin1[0].setAxis(Vec3(0, 0, 1));
	st_pin1[1].setAxis(Vec3(1, 0, 0));
	st_pin1[2].setAxis(Vec3(0, 1, 0));
	st_pin1[3].setCoordinateNames(OpenSim::Array<std::string>("t1_x", 1, 1));
	st_pin1[3].setFunction(new LinearFunction());
	st_pin1[3].setAxis(Vec3(1, 0, 0));
	st_pin1[4].setCoordinateNames(OpenSim::Array<std::string>("t1_y", 1, 1));
	st_pin1[4].setFunction(new LinearFunction());
	st_pin1[4].setAxis(Vec3(0, 1, 0));
	torso_ground = new CustomJoint("torso_ground", model->getGround(), Vec3(0,0,0),
		Vec3(0), *torso, Vec3(0, -0.3,0), Vec3(0,0,0),st_pin1);

	shoulder1 = new PinJoint("shoulder1", *torso, Vec3(0, -0.3, 0),
		Vec3(0), *humer1, Vec3(0, 0.3, 0), Vec3(0));
	elbow1 = new PinJoint("elbow1", *humer1, Vec3(0, -0.3, 0),
		Vec3(0), *ulna1, Vec3(0, 0.3, 0), Vec3(0));
	shoulder2 = new PinJoint("shoulder2", *torso, Vec3(0, -0.3, 0),
		Vec3(0), *humer2, Vec3(0, 0.3, 0), Vec3(0));
	elbow2 = new PinJoint("elbow2", *humer2, Vec3(0, -0.3, 0),
		Vec3(0), *ulna2, Vec3(0, 0.3, 0), Vec3(0));
	hip3 = new PinJoint("hip3", *torso, Vec3(0, 0.5, 0),
		Vec3(0), *femur3, Vec3(0, 0.3, 0), Vec3(0));
	knee3 = new PinJoint("knee3", *femur3, Vec3(0, -0.3, 0),
		Vec3(0), *tibia3, Vec3(0, 0.3, 0), Vec3(0));
	hip4 = new PinJoint("hip4", *torso, Vec3(0, 0.5, 0),
		Vec3(0), *femur4, Vec3(0, 0.3, 0), Vec3(0));
	knee4 = new PinJoint("knee4", *femur4, Vec3(0, -0.3, 0),
		Vec3(0), *tibia4, Vec3(0, 0.3, 0), Vec3(0));
	
    /// Add bodies and joints to model
	model->addBody(torso);
	model->addBody(humer1);
	model->addBody(ulna1);
	model->addBody(humer2);
	model->addBody(ulna2);
	model->addBody(femur3);
	model->addBody(tibia3);
	model->addBody(femur4);
	model->addBody(tibia4);
	
	model->addJoint(torso_ground);
	model->addJoint(shoulder1);
	model->addJoint(elbow1);
	model->addJoint(shoulder2);
	model->addJoint(elbow2);
	model->addJoint(hip3);
	model->addJoint(knee3);
	model->addJoint(hip4);
	model->addJoint(knee4);
	
	// Initialize system and state
	State* state;
	state = new State(model->initSystem());

	std::cout << model->getStateVariableNames() << std::endl;

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);
	std::vector<T> CF(arg[2], arg[2] + 8);

    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
	T ua[NU]; /// joint accelerations (Qdotdots) - controls
	Vec3 F1(0);
	Vec3 F2(0);
	Vec3 F3(0);
	Vec3 F4(0);
	F1[0] = CF[0]; //leg1
	F1[1] = CF[1];
	F2[0] = CF[2]; //leg2
	F2[1] = CF[3];
	F3[0] = CF[4]; //leg 3
	F3[1] = CF[5];
	F4[0] = CF[6]; //leg 4
	F4[1] = CF[7];

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
	std::cout << nbodies << std::endl;
	appliedBodyForces.setToZero();
	/// Set gravity
	Vec3 gravity(0);
	gravity[1] = -9.81;
	/// Set gravity vector
	Vec3 Vec_gravity = Vec3(0, gravity[1], 0);
	/// Add to model
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
			Vec_gravity*model->getBodySet().get(i).getMass(),
            appliedBodyForces);
	}
	/// Add contact forces
	model->getMatterSubsystem().addInStationForce(*state,
		model->getBodySet().get("ulna1").getMobilizedBodyIndex(),
		Vec3(0,-0.3,0),F1,appliedBodyForces);
	model->getMatterSubsystem().addInStationForce(*state,
		model->getBodySet().get("ulna2").getMobilizedBodyIndex(),
		Vec3(0, -0.3, 0), F2, appliedBodyForces);
	model->getMatterSubsystem().addInStationForce(*state,
		model->getBodySet().get("tibia3").getMobilizedBodyIndex(),
		Vec3(0, -0.3, 0), F3, appliedBodyForces);
	model->getMatterSubsystem().addInStationForce(*state,
		model->getBodySet().get("tibia4").getMobilizedBodyIndex(),
		Vec3(0, -0.3, 0), F4, appliedBodyForces);

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
	Recorder CF[8];
	Recorder tau[NR];

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;
	for (int i = 0; i < 8;  ++i) CF[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u, CF };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

	Recorder::stop_recording();

	return 0;

}

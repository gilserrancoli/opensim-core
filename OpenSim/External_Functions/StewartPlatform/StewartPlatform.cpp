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
constexpr int ndof = 6;     // # degrees of freedom
constexpr int NX = ndof*2;  // # states
constexpr int NU = ndof;    // # controls
constexpr int NP = 3;		// # number of output points
constexpr int NR = ndof + NP*3;    // # residual torques + output point coordinates


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
	OpenSim::Body* platform;
	
    /// Joints
	OpenSim::CustomJoint* joint1;

	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	platform = new OpenSim::Body("platform", 508.9, Vec3(0),
            Inertia(412.3, 824.4, 412.3, 0, 0, 0));
	

    /// Joint specifications
	SpatialTransform st_joint1;
	st_joint1[0].setCoordinateNames(OpenSim::Array<std::string>("r1_z", 1, 1));
	st_joint1[0].setFunction(new LinearFunction());
	st_joint1[0].setAxis(Vec3(0, 0, 1));
	st_joint1[1].setCoordinateNames(OpenSim::Array<std::string>("r2_x", 1, 1));
	st_joint1[1].setFunction(new LinearFunction());
	st_joint1[1].setAxis(Vec3(1, 0, 0));
	st_joint1[2].setCoordinateNames(OpenSim::Array<std::string>("r3_y", 1, 1));
	st_joint1[2].setFunction(new LinearFunction());
	st_joint1[2].setAxis(Vec3(0, 1, 0));
	st_joint1[3].setCoordinateNames(OpenSim::Array<std::string>("t4_x", 1, 1));
	st_joint1[3].setFunction(new LinearFunction());
	st_joint1[3].setAxis(Vec3(1, 0, 0));
	st_joint1[4].setCoordinateNames(OpenSim::Array<std::string>("t5_y", 1, 1));
	st_joint1[4].setFunction(new LinearFunction());
	st_joint1[4].setAxis(Vec3(0, 1, 0));
	st_joint1[5].setCoordinateNames(OpenSim::Array<std::string>("t6_z", 1, 1));
	st_joint1[5].setFunction(new LinearFunction());
	st_joint1[5].setAxis(Vec3(0, 0, 1));
	joint1 = new CustomJoint("joint1", model->getGround(), Vec3(0,0,0),
		Vec3(0), *platform, Vec3(0,0,0), Vec3(0,0,0), st_joint1);
	
    /// Add bodies and joints to model
	model->addBody(platform);
	
	model->addJoint(joint1);
	
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

	//Report pos values in absolute frame
	Vec3 Q1_loc = Vec3(1.55, 0, 0.9);
	Vec3 Q2_loc = Vec3(-1.55, 0, 0.9);
	Vec3 Q3_loc = Vec3(0, 0, 1.8);
	Vec3 Q1_abs=platform->findStationLocationInGround(*state, Q1_loc);
	Vec3 Q2_abs = platform->findStationLocationInGround(*state, Q2_loc);
	Vec3 Q3_abs = platform->findStationLocationInGround(*state, Q3_loc);

	// Extract results
    /// Residual forces
	for (int i = 0; i < ndof; ++i) res[0][i] =
        value<T>(residualMobilityForces[i]);
	for (int i = 0; i < 3; ++i) res[0][i + ndof] =
		value<T>(Q1_abs[i]);
	for (int i = 0; i < 3; ++i) res[0][i + ndof + 3] =
		value<T>(Q2_abs[i]);
	for (int i = 0; i < 3; ++i) res[0][i + ndof + 6] =
		value<T>(Q3_abs[i]);

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

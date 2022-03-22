/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Gil Serrancoli
    Contributor: Joris Gillis, Antoine Falisse, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
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
constexpr int ndof = 3;     // # degrees of freedom
constexpr int NX = ndof*2;  // # states
constexpr int NU = ndof;    // # controls
constexpr int NP = 4;       // # parameter values
constexpr int NR = ndof;    // # residual torques

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);
	std::vector<T> p(arg[2], arg[2] + NP);

	T param[NP]; /// parameters to change

	param[0] = p[0];
	param[1] = p[1];
	param[2] = p[2];
	param[3] = p[3];

	T m1 = param[0];
	T s1 = param[1];
	T s2 = param[2];
	T s3 = param[3];
	// OpenSim model: create components
    /// Model
	Model* model;
    /// Bodies
	OpenSim::Body* pend1;
	OpenSim::Body* pend2;
	OpenSim::Body* pend3;
	OpenSim::Body* ball1;
	OpenSim::Body* ball2;
    /// Joints
	OpenSim::PinJoint* pin1;
	OpenSim::PinJoint* pin2;
	OpenSim::PinJoint* pin3;
	OpenSim::WeldJoint* w1;
	OpenSim::WeldJoint* w2;

	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	pend1 = new OpenSim::Body("pend1", 10.0, Vec3(0),
            Inertia(0.01, 13.3, 13.3));
	pend2 = new OpenSim::Body("pend2", 5.0, Vec3(0),
            Inertia(0.01,1.7,1.7));
	pend3 = new OpenSim::Body("pend3", 5.0, Vec3(0),
            Inertia(0.01,1.7,1.7));
	ball1 = new OpenSim::Body("ball1", m1, Vec3(0),
		Inertia((2/5)*m1*(pow(0.5,2))));
	ball2 = new OpenSim::Body("ball2", 10000, Vec3(0),
		Inertia(1000));
    /// Joint specifications
	pin1 = new PinJoint("pin1", model->getGround(), Vec3(0), Vec3(0),
		    *pend1, Vec3(0), Vec3(0));
	pin2 = new PinJoint("pin2", *pend1, Vec3(0, -s1, 0),
            Vec3(0), *pend2, Vec3(0,-s2/2,0), Vec3(0));
	pin3 = new PinJoint("pin3", *pend1, Vec3(0, s3, 0), Vec3(0),
            *pend3, Vec3(0,-1,0), Vec3(0));
	w1 = new WeldJoint("w1", *pend2, Vec3(0, s2/2, 0), Vec3(0),
		*ball1, Vec3(0), Vec3(0));
	w2 = new WeldJoint("w2", *pend3, Vec3(0, 1, 0), Vec3(0),
		*ball2, Vec3(0), Vec3(0));
    /// Add bodies and joints to model
	model->addBody(pend1);
	model->addBody(pend2);
	model->addBody(pend3);
	model->addBody(ball1);
	model->addBody(ball2);
	model->addJoint(pin1);
	model->addJoint(pin2);
	model->addJoint(pin3);
	model->addJoint(w1);
	model->addJoint(w2);

	// Initialize system and state
	State* state;
	state = new State(model->initSystem());

	

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
	Recorder p[NP];
	Recorder tau[NR];

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;
	for (int i = 0; i < NP; ++i) p[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u, p };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

	Recorder::stop_recording();

	return 0;

}

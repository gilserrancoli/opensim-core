/*  This code describes the OpenSim model and the skeleton dynamics of a Rat Hindlimb model
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
#include <OpenSim/Common/MultivariatePolynomialFunction.h>
#include <OpenSim/Common/MultiplierFunction.h>
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
constexpr int ndof = 14;     // # degrees of freedom
constexpr int NX = ndof*2;  // # states
constexpr int NU = ndof;    // # controls
constexpr int NR = ndof;    // # residual torques

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// OpenSim and Simbody use different indices for the states/controls when the
// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
// and legs with pelvis as origin).
// The two following functions allow getting the indices from one reference
// system to the other. These functions are inspired from
// createSystemYIndexMap() in Moco.
// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
// reference system. Note that we only care about the order here so we divide
// by 2 because the states include both Qs and Qdots.
SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
	auto s = model.getWorkingState();
	const auto svNames = model.getStateVariableNames();
	SimTK::Array_<int> idxOSInSimbody(s.getNQ());
	s.updQ() = 0;
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		s.updQ()[iy] = SimTK::NaN;
		const auto svValues = model.getStateVariableValues(s);
		for (int isv = 0; isv < svNames.size(); ++isv) {
			if (SimTK::isNaN(svValues[isv])) {
				s.updQ()[iy] = 0;
				idxOSInSimbody[iy] = isv / 2;
				break;
			}
		}
	}
	return idxOSInSimbody;
}
// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
// reference system.
SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
	auto idxOSInSimbody = getIndicesOSInSimbody(model);
	auto s = model.getWorkingState();
	SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
			if (idxOSInSimbody[iyy] == iy) {
				idxSimbodyInOS[iy] = iyy;
				break;
			}
		}
	}
	return idxSimbodyInOS;
}

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

	// OpenSim model: create components
    /// Model
	Model* model;
    /// Bodies
	OpenSim::Body* spine;
	OpenSim::Body* pelvis;
	OpenSim::Body* femur;
	OpenSim::Body* tibia;
	OpenSim::Body* foot;
	
    /// Joints
	OpenSim::CustomJoint* ground_spine;
	OpenSim::CustomJoint* sacroiliac;
	OpenSim::CustomJoint* hip;
	OpenSim::CustomJoint* knee;
	OpenSim::CustomJoint* ankle;

	// OpenSim model: initialize components
    /// Model
	model = new OpenSim::Model();
    /// Body specifications
	spine = new OpenSim::Body("spine", 0.1145, Vec3(0),
            Inertia(3.7987e-6, 3.7987e-6, 3.7987e-6, 0, 0, 0));
	pelvis = new OpenSim::Body("pelvis", 0.01908, Vec3(0),
			Inertia(3.1655e-6, 2.7699e-6, 2.7699e-6, 0, 0, 0));
	femur = new OpenSim::Body("femur", 0.01326, Vec3(0,-0.01445,0),
            Inertia(9.6e-7, 6.0e-8, 9.6e-7, 0, 0, 0));
	tibia = new OpenSim::Body("tibia", 0.00906, Vec3(0,0.01346,0),
            Inertia(1.4e-6, 3.0e-8, 1.4e-6, 0, 0, 0));
	foot = new OpenSim::Body("foot", 0.00131, Vec3(0.007026,-0.005018,0),
		Inertia(2.6e-8, 4.4e-8, 4.7e-8, 7.6e-9, 0, 0));
	

    /// Joint specifications
	SpatialTransform st_ground_spine;
	st_ground_spine[0].setCoordinateNames(OpenSim::Array<std::string>("sacrum_pitch",1,1));
	st_ground_spine[0].setFunction(new LinearFunction());
	st_ground_spine[0].setAxis(Vec3(0, 0, 1));
	st_ground_spine[1].setCoordinateNames(OpenSim::Array<std::string>("sacrum_roll", 1, 1));
	st_ground_spine[1].setFunction(new LinearFunction());
	st_ground_spine[1].setAxis(Vec3(1, 0, 0));
	st_ground_spine[2].setCoordinateNames(OpenSim::Array<std::string>("sacrum_yaw", 1, 1));
	st_ground_spine[2].setFunction(new LinearFunction());
	st_ground_spine[2].setAxis(Vec3(0, 1, 0));
	st_ground_spine[3].setCoordinateNames(OpenSim::Array<std::string>("sacrum_x", 1, 1));
	st_ground_spine[3].setFunction(new LinearFunction());
	st_ground_spine[3].setAxis(Vec3(1, 0, 0));
	st_ground_spine[4].setCoordinateNames(OpenSim::Array<std::string>("sacrum_y", 1, 1));
	st_ground_spine[4].setFunction(new LinearFunction());
	st_ground_spine[4].setAxis(Vec3(0, 1, 0));
	st_ground_spine[5].setCoordinateNames(OpenSim::Array<std::string>("sacrum_z", 1, 1));
	st_ground_spine[5].setFunction(new LinearFunction());
	st_ground_spine[5].setAxis(Vec3(0, 0, 1));
	ground_spine = new CustomJoint("ground_spine", model->getGround(), Vec3(0), Vec3(0),
		    *spine, Vec3(0,0,0), Vec3(0,0,0),st_ground_spine);


	SpatialTransform st_sacroiliac;
	st_sacroiliac[0].setCoordinateNames(OpenSim::Array<std::string>("sacroiliac_flx", 1, 1));
	st_sacroiliac[0].setFunction(new Constant(3.7*Pi/180));
	st_sacroiliac[0].setAxis(Vec3(0, 0, 1));
	st_sacroiliac[1].setAxis(Vec3(1, 0, 0));
	st_sacroiliac[2].setAxis(Vec3(0, 1, 0));
	sacroiliac = new CustomJoint("sacroiliac", *spine, Vec3(0.004553, -0.002732, 0.009107),
		Vec3(0), *pelvis, Vec3(0,0,0), Vec3(0,0,0),st_sacroiliac);

	SpatialTransform st_hip;
	st_hip[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flx", 1, 1));
	st_hip[0].setFunction(new LinearFunction());
	st_hip[0].setAxis(Vec3(0, 0, 1));
	st_hip[1].setCoordinateNames(OpenSim::Array<std::string>("hip_add", 1, 1));
	st_hip[1].setFunction(new LinearFunction());
	st_hip[1].setAxis(Vec3(1, 0, 0));
	st_hip[2].setCoordinateNames(OpenSim::Array<std::string>("hip_int", 1, 1));
	st_hip[2].setFunction(new LinearFunction());
	st_hip[2].setAxis(Vec3(0, 1, 0));
	hip = new CustomJoint("hip", *pelvis, Vec3(0),
		Vec3(0), *femur, Vec3(0), Vec3(0),st_hip);

	SpatialTransform st_knee;
	st_knee[0].setAxis(Vec3(0, 0, 1));
	st_knee[0].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	st_knee[0].setFunction(new LinearFunction());
	st_knee[1].setAxis(Vec3(1, 0, 0));
	st_knee[1].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	st_knee[1].setFunction(new LinearFunction(0.0667,0.0874)); // This is Y=0.0667X+0.0874
	st_knee[2].setAxis(Vec3(0, 1, 0));
	st_knee[2].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	st_knee[2].setFunction(new Constant(0.2618)); // This is Y=0.2618 rad
	st_knee[3].setAxis(Vec3(1, 0, 0));
	st_knee[3].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	osim_double_adouble coefficients_knee_tx[3] = {0.0021, 0.0549, 0.0171}; //This is Y=0.0171X^2+0.0549X+0.0021
	Vector coefficients_tx_vec(3);
	for (int i = 0; i < 3; ++i) coefficients_tx_vec[i] = coefficients_knee_tx[i];
	auto* multivar_knee_tx = new MultivariatePolynomialFunction();
	multivar_knee_tx->setDimension(1);
	multivar_knee_tx->setOrder(2);
	multivar_knee_tx->setCoefficients(coefficients_tx_vec);
	st_knee[3].setFunction(new MultiplierFunction(multivar_knee_tx, 0.9637));
	st_knee[4].setAxis(Vec3(0, 1, 0));
	st_knee[4].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	osim_double_adouble coefficients_knee_ty[4] = {-0.0721,0.0005,0.0248,0.0053}; //This is Y=0.0053X^3+0.0248X^2+0.0005X-0.0721;
	Vector coefficients_ty_vec(4);
	for (int i = 0; i < 4; ++i) coefficients_ty_vec[i] = coefficients_knee_ty[i];
	auto* multivar_knee_ty = new MultivariatePolynomialFunction();
	multivar_knee_ty->setDimension(1);
	multivar_knee_ty->setOrder(3);
	multivar_knee_ty->setCoefficients(coefficients_ty_vec);
	st_knee[4].setFunction(new MultiplierFunction(multivar_knee_ty, 0.9637));
	st_knee[5].setAxis(Vec3(0, 0, 1));
	st_knee[5].setCoordinateNames(OpenSim::Array<std::string>("knee_flx", 1, 1));
	osim_double_adouble coefficients_knee_tz[4] = {-0.0039,0.0,0.0022,0.0005}; //This is Y=0.0005X^3+0.0022*X^2-0.0039;
	Vector coefficients_tz_vec(4);
	for (int i = 0; i < 4; ++i) coefficients_tz_vec[i] = coefficients_knee_tz[i];
	auto* multivar_knee_tz = new MultivariatePolynomialFunction();
	multivar_knee_tz->setDimension(1);
	multivar_knee_tz->setOrder(3);
	multivar_knee_tz->setCoefficients(coefficients_tz_vec);
	st_knee[5].setFunction(new MultiplierFunction(multivar_knee_tz, 0.9637));
	knee = new CustomJoint("knee", *femur, Vec3(0, 0, 0),
		Vec3(0), *tibia, Vec3(0,0,0), Vec3(0, 0, 0), st_knee);

	SpatialTransform st_ankle;
	st_ankle[0].setAxis(Vec3(0, 0, 1));
	st_ankle[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_flx", 1, 1));
	st_ankle[0].setFunction(new LinearFunction());
	st_ankle[1].setAxis(Vec3(1, 0, 0));
	st_ankle[1].setCoordinateNames(OpenSim::Array<std::string>("ankle_add", 1, 1));
	st_ankle[1].setFunction(new LinearFunction());
	st_ankle[2].setAxis(Vec3(0, 1, 0));
	st_ankle[2].setCoordinateNames(OpenSim::Array<std::string>("ankle_int", 1, 1));
	st_ankle[2].setFunction(new LinearFunction());
	ankle = new CustomJoint("ankle", *tibia, Vec3(0, 0, 0),
		Vec3(0), *foot, Vec3(0, 0, 0), Vec3(0),st_ankle);
	
    /// Add bodies and joints to model
	model->addBody(spine);
	model->addBody(pelvis);
	model->addBody(femur);
	model->addBody(tibia);
	model->addBody(foot);
	model->addJoint(ground_spine);
	model->addJoint(sacroiliac);
	model->addJoint(hip);
	model->addJoint(knee);
	model->addJoint(ankle);
	
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
	 /// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];
	std::cout << indicesOSInSimbody << std::endl;
	std::cout << model->getStateVariableNames() << std::endl;

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
	
	/// Add to model
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
			gravity *model->getBodySet().get(i).getMass(),
            appliedBodyForces);
	}
    /// knownUdot
	Vector knownUdot(ndof);
	knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
	Vector residualMobilityForces(ndof);
	residualMobilityForces.setToZero();

	std::cout << "Applied mob forces " << appliedMobilityForces.size() << std::endl;


	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

	// Extract results
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    /// Residual forces
	for (int i = 0; i < NR; ++i) res[0][i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);

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

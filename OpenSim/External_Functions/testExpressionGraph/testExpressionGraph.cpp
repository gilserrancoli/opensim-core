/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
#include "SimTKcommon/internal/recorder.h"
#include <initializer_list>
#include <unordered_map>

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
    joint accelerations (controls u), contact forces as well as several
    parameters related to the contact models (parameters p), and returns the
    joint torques as well as several variables for use in the
    optimal control problems. F is templatized using type T. F(x,u,p)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 1;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int NX = 2;      // # n inputs
constexpr int NR =1;    // # outputs

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    // create components
    Real x1;
    Real x2;
    
    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);

    x1 = x[0];
    x2 = x[1];

    Real y = cos(x2) - x1 * x2;
     
    // Extract results
  
    res[0][0] = value<T>(y);

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
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}

/*  This code describes the OpenSim model, it is based on the code, where
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

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace SimTK;
using namespace OpenSim;

/*  The function F describes the OpenSim model. F takes as inputs joint
    positions and velocities (states x), and returns several variables for use
    when computing the contact forces in the optimal control problems. F is
    templatized using type T. F(x)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 1;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 34;        // # degrees of freedom (excluding locked)
constexpr int NX = ndof*2;      // # states
constexpr int NR = 84;          // # output variables 84

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
    OpenSim::Model* model;
    /// Bodies
    OpenSim::Body* pelvis;
    OpenSim::Body* femur_r;
    OpenSim::Body* femur_l;
    OpenSim::Body* tibia_r;
    OpenSim::Body* tibia_l;
    OpenSim::Body* femoral_component;
    OpenSim::Body* tibial_tray;
    OpenSim::Body* talus_r;
    OpenSim::Body* talus_l;
    OpenSim::Body* calcn_r;
    OpenSim::Body* calcn_l;
    OpenSim::Body* toes_r;
    OpenSim::Body* toes_l;
    OpenSim::Body* torso;
    OpenSim::Body* humerus_r;
    OpenSim::Body* humerus_l;
    OpenSim::Body* ulna_r;
    OpenSim::Body* ulna_l;
    OpenSim::Body* radius_r;
    OpenSim::Body* radius_l;
    OpenSim::Body* hand_r;
    OpenSim::Body* hand_l;
    /// Joints
    OpenSim::CustomJoint* ground_pelvis;
    OpenSim::CustomJoint* hip_r;
    OpenSim::CustomJoint* hip_l;
    OpenSim::CustomJoint* knee_r;
    OpenSim::CustomJoint* knee_l;
    OpenSim::CustomJoint* ankle_r;
    OpenSim::CustomJoint* ankle_l;
    OpenSim::CustomJoint* subtalar_r;
    OpenSim::CustomJoint* subtalar_l;
    OpenSim::WeldJoint* mtp_r;
    OpenSim::WeldJoint* mtp_l;
    OpenSim::CustomJoint* back;
    OpenSim::CustomJoint* shoulder_r;
    OpenSim::CustomJoint* shoulder_l;
    OpenSim::CustomJoint* elbow_r;
    OpenSim::CustomJoint* elbow_l;
    OpenSim::CustomJoint* radioulnar_r;
    OpenSim::CustomJoint* radioulnar_l;
    OpenSim::WeldJoint* radius_hand_r;
    OpenSim::WeldJoint* radius_hand_l;
    OpenSim::WeldJoint* femoral_component_weld;
    OpenSim::WeldJoint* tibial_tray_weld;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 9.53459239861278, Vec3(-0.0672, 0, 0), Inertia(0.0832132493453181, 0.0705026541156487, 0.0469103262792837, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 7.53034185009555, Vec3(0, -0.17, 0), Inertia(0.108360110411211, 0.0284388137872461, 0.114303913935877, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 7.53034185009555, Vec3(0, -0.17, 0), Inertia(0.108360110411211, 0.0284388137872461, 0.114303913935877, 0, 0, 0));
    femoral_component = new OpenSim::Body("femoral_component", 0.00914431311487012, Vec3(0, 0, 0), Inertia(9.14431311487012e-005, 9.14431311487012e-005, 9.14431311487012e-005, 0, 0, 0));
    tibial_tray = new OpenSim::Body("tibial_tray", 0.00914431311487012, Vec3(0, 0, 0), Inertia(9.14431311487012e-005, 9.14431311487012e-005, 9.14431311487012e-005, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 3.00162077995612, Vec3(0, -0.2054, 0), Inertia(0.0407836364923208, 0.00411494090169155, 0.0413322952792129, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 3.00162077995612, Vec3(0, -0.2054, 0), Inertia(0.0407836364923208, 0.00411494090169155, 0.0413322952792129, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.0809271710666006, Vec3(0, 0, 0), Inertia(0.000809600905938141, 0.000809600905938141, 0.000809600905938141, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.0809271710666006, Vec3(0, 0, 0), Inertia(0.000809600905938141, 0.000809600905938141, 0.000809600905938141, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 1.01200113242268, Vec3(0.1, 0.03, 0), Inertia(0.00109731757378441, 0.00320050959020454, 0.00329195272135324, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 1.01200113242268, Vec3(0.1, 0.03, 0), Inertia(0.00109731757378441, 0.00320050959020454, 0.00329195272135324, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.175387925543209, Vec3(0.0346, 0.006, 0.0175), Inertia(8.09637483190601e-005, 0.000161918352325005, 0.000809600905938141, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.175387925543209, Vec3(0.0346, 0.006, -0.0175), Inertia(8.09637483190601e-005, 0.000161918352325005, 0.000809600905938141, 0, 0, 0));
    torso = new OpenSim::Body("torso", 20.7883858994007, Vec3(-0.0296672, 0.31645, 0), Inertia(1.08198530190117, 0.501693302382011, 1.08198530190117, 0, 0, 0));
    humerus_l = new OpenSim::Body("humerus_l", 2.17645574940628, Vec3(0, -0.173391, 0), Inertia(0.0142119584950294, 0.00490268549790861, 0.0159524653825422, 0, 0, 0));
    humerus_r = new OpenSim::Body("humerus_r", 2.17645574940628, Vec3(0, -0.173391, 0), Inertia(0.0142119584950294, 0.00490268549790861, 0.0159524653825422, 0, 0, 0));
    ulna_l = new OpenSim::Body("ulna_l", 0.763119972691866, Vec3(0, -0.133981, 0), Inertia(0.00459792509391667, 0.000959324006765869, 0.00498755345265168, 0, 0, 0));
    ulna_r = new OpenSim::Body("ulna_r", 0.763119972691866, Vec3(0, -0.133981, 0), Inertia(0.00459792509391667, 0.000959324006765869, 0.00498755345265168, 0, 0, 0));
    radius_l = new OpenSim::Body("radius_l", 0.763119972691866, Vec3(0, -0.133981, 0), Inertia(0.00459792509391667, 0.000959324006765869, 0.00498755345265168, 0, 0, 0));
    radius_r = new OpenSim::Body("radius_r", 0.763119972691866, Vec3(0, -0.133981, 0), Inertia(0.00459792509391667, 0.000959324006765869, 0.00498755345265168, 0, 0, 0));
    hand_l = new OpenSim::Body("hand_l", 0.57469528807659, Vec3(0, -0.0756973, 0), Inertia(0.00138465536251643, 0.000849110407283059, 0.00208008765221078, 0, 0, 0));
    hand_r = new OpenSim::Body("hand_r", 0.57469528807659, Vec3(0, -0.0756973, 0), Inertia(0.00138465536251643, 0.000849110407283059, 0.00208008765221078, 0, 0, 0));

    /// Joints
    /// Ground-Pelvis transform
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
    st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
    st_ground_pelvis[3].setFunction(new LinearFunction());
    st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
    st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
    st_ground_pelvis[4].setFunction(new LinearFunction());
    st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
    st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
    st_ground_pelvis[5].setFunction(new LinearFunction());
    st_ground_pelvis[5].setAxis(Vec3(0, 0, 1));
    /// Hip_l transform
    SpatialTransform st_hip_l;
    st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
    st_hip_l[0].setFunction(new LinearFunction());
    st_hip_l[0].setAxis(Vec3(0, 0, 1));
    st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
    st_hip_l[1].setFunction(new LinearFunction());
    st_hip_l[1].setAxis(Vec3(-1, 0, 0));
    st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
    st_hip_l[2].setFunction(new LinearFunction());
    st_hip_l[2].setAxis(Vec3(0, -1, 0));
    /// Hip_r transform
    SpatialTransform st_hip_r;
    st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
    st_hip_r[0].setFunction(new LinearFunction());
    st_hip_r[0].setAxis(Vec3(0, 0, 1));
    st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
    st_hip_r[1].setFunction(new LinearFunction());
    st_hip_r[1].setAxis(Vec3(1, 0, 0));
    st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
    st_hip_r[2].setFunction(new LinearFunction());
    st_hip_r[2].setAxis(Vec3(0, 1, 0));
    /// Knee_l transform
    SpatialTransform st_knee_l;
    st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    st_knee_l[2].setFunction(new LinearFunction());
    st_knee_l[2].setAxis(Vec3(0, 0, -1));
    /// Knee_r transform
    SpatialTransform st_knee_r; // should it be flexion, adduction, rotation?
    st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[0].setFunction(new LinearFunction());
    st_knee_r[0].setAxis(Vec3(0, 0, 1));
    st_knee_r[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_r", 1, 1));
    st_knee_r[1].setFunction(new LinearFunction());
    st_knee_r[1].setAxis(Vec3(1, 0, 0));
    st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_rotation_r", 1, 1));
    st_knee_r[2].setFunction(new LinearFunction());
    st_knee_r[2].setAxis(Vec3(0, -1, 0));
    st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_tx_r", 1, 1));
    st_knee_r[3].setFunction(new LinearFunction());
    st_knee_r[3].setAxis(Vec3(1, 0, 0));
    st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_ty_r", 1, 1));
    st_knee_r[4].setFunction(new LinearFunction());
    st_knee_r[4].setAxis(Vec3(0, 1, 0));
    st_knee_r[5].setCoordinateNames(OpenSim::Array<std::string>("knee_tz_r", 1, 1));
    st_knee_r[5].setFunction(new LinearFunction());
    st_knee_r[5].setAxis(Vec3(0, 0, 1));

    /// Ankle_l transform
    SpatialTransform st_ankle_l;
    st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
    st_ankle_l[0].setFunction(new LinearFunction());
    st_ankle_l[0].setAxis(Vec3(0.16149, 0.027235, 0.986));
    /// Ankle_r transform
    SpatialTransform st_ankle_r;
    st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
    st_ankle_r[0].setFunction(new LinearFunction());
    st_ankle_r[0].setAxis(Vec3(-0.16149, -0.027235, 0.986));
    /// Subtalar_l transform
    SpatialTransform st_subtalar_l;
    st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
    st_subtalar_l[0].setFunction(new LinearFunction());
    st_subtalar_l[0].setAxis(Vec3(-0.72089, -0.68935, -0.073236));
    /// Subtalar_r transform
    SpatialTransform st_subtalar_r;
    st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
    st_subtalar_r[0].setFunction(new LinearFunction());
    st_subtalar_r[0].setAxis(Vec3(0.72089, 0.68935, -0.073236));
    /// Back transform
    SpatialTransform st_back;
    st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
    st_back[0].setFunction(new LinearFunction());
    st_back[0].setAxis(Vec3(0, 0, 1));
    st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
    st_back[1].setFunction(new LinearFunction());
    st_back[1].setAxis(Vec3(1, 0, 0));
    st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
    st_back[2].setFunction(new LinearFunction());
    st_back[2].setAxis(Vec3(0, 1, 0));
    /// Shoulder_l transform
    SpatialTransform st_sho_l;
    st_sho_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
    st_sho_l[0].setFunction(new LinearFunction());
    st_sho_l[0].setAxis(Vec3(0, 0, 1));
    st_sho_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
    st_sho_l[1].setFunction(new LinearFunction());
    st_sho_l[1].setAxis(Vec3(-1, 0, 0));
    st_sho_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
    st_sho_l[2].setFunction(new LinearFunction());
    st_sho_l[2].setAxis(Vec3(0, -1, 0));
    /// Shoulder_r transform
    SpatialTransform st_sho_r;
    st_sho_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
    st_sho_r[0].setFunction(new LinearFunction());
    st_sho_r[0].setAxis(Vec3(0, 0, 1));
    st_sho_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
    st_sho_r[1].setFunction(new LinearFunction());
    st_sho_r[1].setAxis(Vec3(1, 0, 0));
    st_sho_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
    st_sho_r[2].setFunction(new LinearFunction());
    st_sho_r[2].setAxis(Vec3(0, 1, 0));
    /// Elbow_l transform
    SpatialTransform st_elb_l;
    st_elb_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
    st_elb_l[0].setFunction(new LinearFunction());
    st_elb_l[0].setAxis(Vec3(-0.22604696, -0.022269, 0.97386183));
    /// Elbow_r transform
    SpatialTransform st_elb_r;
    st_elb_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
    st_elb_r[0].setFunction(new LinearFunction());
    st_elb_r[0].setAxis(Vec3(0.22604696, 0.022269, 0.97386183));
    /// Radioulnar_l transform
    SpatialTransform st_radioulnar_l;
    st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
    st_radioulnar_l[0].setFunction(new LinearFunction());
    st_radioulnar_l[0].setAxis(Vec3(-0.05639803, -0.99840646, 0.001952));
    /// Radioulnar_r transform
    SpatialTransform st_radioulnar_r;
    st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
    st_radioulnar_r[0].setFunction(new LinearFunction());
    st_radioulnar_r[0].setAxis(Vec3(0.05639803, 0.99840646, 0.001952));
    /// Joint specifications
    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.041283, -0.098, -0.084432), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.041283, -0.098, 0.084432), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *tibia_l, Vec3(0), Vec3(0), *femur_l, Vec3(0, -0.39335, 0), Vec3(0), st_knee_l);
    femoral_component_weld = new WeldJoint("femoral_component_weld", *femur_r, Vec3(0, -0.39335, 0), Vec3(3.083, 0, -3.083), *femoral_component, Vec3(0), Vec3(0));
    knee_r = new CustomJoint("knee_r", *tibial_tray, Vec3(0, 0, 0), Vec3(0), *femoral_component, Vec3(0), Vec3(0), st_knee_r);
    tibial_tray_weld = new WeldJoint("tibial_tray_weld", *tibial_tray, Vec3(0, 0.044254, 0), Vec3(-3.14159, -7.34641e-06, -3.14159), *tibia_r, Vec3(0, 0, 0), Vec3(0, 0, 0));
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.44751, 0), Vec3(-0.041214, 0.0031538, -0.050218), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.44751, 0), Vec3(-0.041214, 0.0031538, -0.050218), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.043062, -0.04869, -0.018), Vec3(0.11332, -0.12598, 0.019146), *calcn_l, Vec3(0), Vec3(0), st_subtalar_l);
    subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.043062, -0.04869, 0.018), Vec3(-0.11332, 0.12598, 0.019146), *calcn_r, Vec3(0), Vec3(0), st_subtalar_r);
    mtp_l = new WeldJoint("mtp_l", *calcn_l, Vec3(0.17664, 0.0043199, -0.0021898), Vec3(-0.20251, -0.096222, 0.23426), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new WeldJoint("mtp_r", *calcn_r, Vec3(0.17664, 0.0043199, 0.0021898), Vec3(0.20251, 0.096222, 0.23426), *toes_r, Vec3(0), Vec3(0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.1007, 0.0815, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.00312000085692857, 0.367378864769877, -0.147313405645013), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
    shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.00312000085692857, 0.367378864769877, 0.147313405645013), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
    elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.0138542690293625, -0.301742480054983, 0.0101134899069334), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
    elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.0138542690293625, -0.301742480054983, -0.0101134899069334), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
    radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00747801922778799, -0.0144591342494185, -0.0289949718326734), Vec3(0), *radius_l, Vec3(0), Vec3(0), st_radioulnar_l);
    radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00747801922778799, -0.0144591342494185, 0.0289949718326734), Vec3(0), *radius_r, Vec3(0), Vec3(0), st_radioulnar_r);
    radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.00977911924287958, -0.262170883410249, -0.0151294546885974), Vec3(0), *hand_l, Vec3(0), Vec3(0));
    radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.00977911924287958, -0.262170883410249, 0.0151294546885974), Vec3(0), *hand_r, Vec3(0), Vec3(0));

    /// Add bodies and joints to model
    model->addBody(pelvis);		        model->addJoint(ground_pelvis);
    model->addBody(femur_l);		    model->addJoint(hip_l);
    model->addBody(femur_r);		    model->addJoint(hip_r);
    model->addBody(tibia_l);		    model->addJoint(knee_l);
    model->addBody(femoral_component);  model->addJoint(femoral_component_weld);
    model->addBody(tibial_tray);        model->addJoint(knee_r);
    model->addBody(tibia_r);		    model->addJoint(tibial_tray_weld);
    model->addBody(talus_l);		    model->addJoint(ankle_l);
    model->addBody(talus_r);		    model->addJoint(ankle_r);
    model->addBody(calcn_l);		    model->addJoint(subtalar_l);
    model->addBody(calcn_r);		    model->addJoint(subtalar_r);
    model->addBody(toes_l);		        model->addJoint(mtp_l);
    model->addBody(toes_r);		        model->addJoint(mtp_r);
    model->addBody(torso);              model->addJoint(back);
    model->addBody(humerus_l);          model->addJoint(shoulder_l);
    model->addBody(humerus_r);          model->addJoint(shoulder_r);
    model->addBody(ulna_l);             model->addJoint(elbow_l);
    model->addBody(ulna_r);             model->addJoint(elbow_r);
    model->addBody(radius_l);           model->addJoint(radioulnar_l);
    model->addBody(radius_r);           model->addJoint(radioulnar_r);
    model->addBody(hand_l);             model->addJoint(radius_hand_l);
    model->addBody(hand_r);             model->addJoint(radius_hand_r);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    std::cout << model->getStateVariableNames().size() << std::endl;
    for (int i = 0; i < model->getStateVariableNames().size(); i++) {
        std::cout << i << " " << model->getStateVariableNames().get(i) << std::endl;
    }
    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);

    // States and controls
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// pro_sup dofs are locked so Qs and Qdots are hard coded
    QsUs[NX] = 1.51;
    QsUs[NX+1] = 0;
    QsUs[NX+2] = 1.51;
    QsUs[NX+3] = 0;

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Get position, velocity and transform of calcn and toes
    SpatialVec vel_calcn_l =  calcn_l->getVelocityInGround(*state);
    SpatialVec vel_calcn_r =  calcn_r->getVelocityInGround(*state);
    Vec3 pos_calcn_l = calcn_l->getPositionInGround(*state);
    Vec3 pos_calcn_r = calcn_r->getPositionInGround(*state);
    Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
    Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);

    SpatialVec vel_toes_l =  toes_l->getVelocityInGround(*state);
    SpatialVec vel_toes_r =  toes_r->getVelocityInGround(*state);
    Vec3 pos_toes_l = toes_l->getPositionInGround(*state);
    Vec3 pos_toes_r = toes_r->getPositionInGround(*state);
    Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);
    Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);

    // Extract results
    int nc = 3;
    for (int i = 0; i < nc; ++i) {
        res[0][i] = value<T>(vel_calcn_l[0][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i+nc] = value<T>(vel_calcn_l[1][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i+nc+nc] = value<T>(vel_calcn_r[0][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i+nc+nc+nc] = value<T>(vel_calcn_r[1][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i+nc+nc+nc+nc] = value<T>(pos_calcn_l[i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i+nc+nc+nc+nc+nc] = value<T>(pos_calcn_r[i]);
    }
    int count = 0;
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            res[0][nc+nc+nc+nc+nc+nc+count] = TR_GB_calcn_l.R().get(i,j);
            ++ count;
        }
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+i] = TR_GB_calcn_l.T().get(i);
    }
    int count2 = 0;
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+count2] = TR_GB_calcn_r.R().get(i,j);
            ++ count2;
        }
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+i] = TR_GB_calcn_r.T().get(i);
    }

    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i] = value<T>(vel_toes_l[0][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i+nc] = value<T>(vel_toes_l[1][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i+nc+nc] = value<T>(vel_toes_r[0][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i+nc+nc+nc] = value<T>(vel_toes_r[1][i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i+nc+nc+nc+nc] = value<T>(pos_toes_l[i]);
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+i+nc+nc+nc+nc+nc] = value<T>(pos_toes_r[i]);
    }
    int count3 = 0;
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+nc+nc+nc+nc+nc+nc+count3] = TR_GB_toes_l.R().get(i,j);
            ++ count3;
        }
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+nc+nc+nc+nc+nc+nc+nc*nc+i] = TR_GB_toes_l.T().get(i);
    }
    int count4 = 0;
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+nc+nc+nc+nc+nc+nc+nc*nc+nc+count4] = TR_GB_toes_r.R().get(i,j);
            ++ count4;
        }
    }
    for (int i = 0; i < nc; ++i) {
        res[0][nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+nc+nc+nc+nc+nc+nc+nc+nc*nc+nc+nc*nc+i] = TR_GB_toes_r.T().get(i);
    }

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

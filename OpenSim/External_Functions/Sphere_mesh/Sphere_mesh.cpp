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
constexpr int n_in = 3;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 6;        // # degrees of freedom (excluding locked)
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NP = 2;          // # parameters, 1) E, 2) poisson, 3) method to compute pressures
constexpr int NR = ndof+3+3;    // # residual torques + # GRFs + # GRMs

constexpr int numpairs = 800; // all faces of the sphere
constexpr int nfacesSphere = 800;
//constexpr int nfacesFem = 188;
constexpr const char radForPairs[] = "05"; // 1 is 1 cm, 05 is 0.5 cm
constexpr char* multiplier_method = "cylinders"; // multiplier method: "cylinders" or "spheres"
constexpr char* method_pen = "mesh"; // "mesh", "centerdist"
constexpr char* mod = "HC"; //HC or EF

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

std::unordered_map<std::string, int>
createSystemYIndexMap(const Model& model) {
    std::unordered_map<std::string, int> sysYIndices;
    auto s = model.getWorkingState();
    const auto svNames = model.getStateVariableNames();
    s.updY() = 0;
    for (int iy = 0; iy < s.getNY(); ++iy) {
        s.updY()[iy] = SimTK::NaN;
        const auto svValues = model.getStateVariableValues(s);
        for (int isv = 0; isv < svNames.size(); ++isv) {
            if (SimTK::isNaN(svValues[isv])) {
                sysYIndices[svNames[isv]] = iy;
                s.updY()[iy] = 0;
                break;
            }
        }
    }
    SimTK_ASSERT2_ALWAYS(svNames.size() == (int)sysYIndices.size(),
        "Expected to find %i state indices but found %i.", svNames.size(),
        sysYIndices.size());
    return sysYIndices;
}

std::vector<Vec4> ReadDataDoublex4columns(std::string filename) {


    std::ifstream  file(filename);

    //
    double v1, v2, v3, v4;
    std::vector<Vec4> out_csv;

    if (file.is_open()) {

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            std::vector<double> values;

            while (std::getline(iss, token, ',')) {
                values.push_back(std::stod(token)); // Convert string to double and add to values vector
            }

            if (values.size() >= 4) {
                out_csv.emplace_back(values[0], values[1], values[2], values[3]); // Create Vec4 object and add to out_csv
            }
        }

        file.close();
    }
    else {
        std::cerr << "Error opening file." << std::endl;
    }
    return out_csv;
}

std::vector<std::vector<int>> ReadDataIntx3columns(std::string filename) {


    std::ifstream  file(filename);

    int v1, v2, v3;
    std::vector<std::vector<int>> out_csv;

    if (file.is_open()) {

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            std::vector<int> row;

            while (std::getline(iss, token, ',')) {
                row.push_back(std::stoi(token)); // Convert string to integer and add to the row vector
            }

            out_csv.push_back(row); // Add the row to the out_csv vector
        }

        file.close();

    }
    else {
        std::cerr << "Error opening file." << std::endl;
    }
    return out_csv;
}

std::vector<std::vector<int>> ReadDataIntx2columns(std::string filename) {


    std::ifstream  file(filename);

    int v1, v2;
    std::vector<std::vector<int>> out_csv;
    if (file.is_open()) {
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            std::vector<int> row;

            while (std::getline(iss, token, ',')) {
                row.push_back(std::stoi(token)); // Convert string to integer and add to the row vector
            }

            out_csv.push_back(row); // Add the row to the out_csv vector
        }
        file.close();
        }else {
        std::cerr << "Error opening file." << std::endl;
    }
    return out_csv;
}

std::vector<Vec3> ReadDataDoublex3columns(std::string filename) {

    
    std::ifstream  file(filename);

    //
    double v1, v2, v3;
    std::vector<Vec3> out_csv;
    
    int nrows = 0;
    if (file.is_open()) {

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);

            std::string token;
            std::vector<std::string> tokens;

            // Split the line by comma and store tokens in a vector
            while (std::getline(iss, token, ',')) {
                tokens.push_back(token);
            }

            // Check if the line has at least three tokens (assuming three values in each row)
            if (tokens.size() >= 3) {
                double v1 = std::stod(tokens[0]);
                double v2 = std::stod(tokens[1]);
                double v3 = std::stod(tokens[2]);

                out_csv.emplace_back(v1, v2, v3); // Add Vec3 object to the vector
            }

            nrows = nrows + 1;

        }

        file.close();
    }else {
        std::cerr << "Error opening file." << std::endl;
    }
  
    return out_csv;
}

Mat33 R_aux(Vec3 sphere_trans, Vec3 sphere_rot) {
    Recorder psi = sphere_rot[0];
    Recorder theta = sphere_rot[1];
    Recorder phi = sphere_rot[2];

    Mat33 R1(0);
    R1.set(0, 0, cos(psi));
    R1.set(0, 1, -sin(psi));
    R1.set(0, 2, 0.0);
    R1.set(1, 0, sin(psi));
    R1.set(1, 1, cos(psi));
    R1.set(1, 2, 0.0);
    R1.set(2, 0, 0.0);
    R1.set(2, 1, 0.0);
    R1.set(2, 2, 1.0);
    Mat33 R2(0);
    R2.set(0, 0, 1.0);
    R2.set(0, 1, 0.0);
    R2.set(0, 2, 0.0);
    R2.set(1, 0, 0.0);
    R2.set(1, 1, cos(theta));
    R2.set(1, 2, -sin(theta));
    R2.set(2, 0, 0.0);
    R2.set(2, 1, sin(theta));
    R2.set(2, 2, cos(theta));
    Mat33 R3(0);
    R3.set(0, 0, cos(phi));
    R3.set(0, 1, 0.0);
    R3.set(0, 2, sin(phi));
    R3.set(1, 0, 0.0);
    R3.set(1, 1, 1.0);
    R3.set(1, 2, 0.0);
    R3.set(2, 0, -sin(phi));
    R3.set(2, 1, 0.0);
    R3.set(2, 2, cos(phi));
    Mat33 R = R1*R2*R3;

    return R;
}

Mat44 ftransf_function(Vec3 sphere_trans, Vec3 sphere_rot) {

    Mat33 R = R_aux(sphere_trans, sphere_rot);
    
    /*Vec3 aux = R*(-sphere_trans);*/
    
    //Mat44 Rtrans0042(1);
    //Rtrans0042.set(1, 3, 0.042);
        
    //Mat44 Rtranstib_4x4(1);
    //Rtranstib_4x4.set(0, 3, aux[0]);
    //Rtranstib_4x4.set(1, 3, aux[1]);
    //Rtranstib_4x4.set(2, 3, aux[2]);

    //std::cout << Rtranstib_4x4 << std::endl;

    Mat44 R_rot(1);
    R_rot.setToZero();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_rot.set(i, j, R[i][j]);
        }
    }
    R_rot.set(3, 3, 1.0);
    R_rot.set(0, 3, sphere_trans[0]);
    R_rot.set(1, 3, sphere_trans[1]);
    R_rot.set(2, 3, sphere_trans[2]);
    //Rtranstib_4x4.set(1, 3, aux[1]);
    //Rtranstib_4x4.set(2, 3, aux[2]);

    std::cout << R_rot << std::endl;
    //Mat44 R_tib = Rtranstib_4x4*R_rot;
    return R_rot;
}

Real CheckContact(Real overlap) {
    Real k = 1e3; // k subject to change
    Real multiplier = (tanh(k * overlap) + 1.0) / 2.0;
    /*std::cout << multiplier << std::endl;*/
    return multiplier;
}

//Vector_<Real> GenerateMultList_cylinders(std::vector<std::vector<int>> pairs_list, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf, Vector_<Vec4> tibplanes, Vector_<Vec4> femplanes, Vector_<Real> radfem) {
//    std::cout << pairs_list.size() << std::endl;
//    Vector_<Real> mlist(pairs_list.size(), 0.0);
//    for (int i = 0; i < pairs_list.size(); i++) {
//        Vec3 arg1 = Vec3(femplanes[pairs_list[i][1]-1][0], femplanes[pairs_list[i][1]-1][1], femplanes[pairs_list[i][1]-1][2]);
//        Vec3 arg2 = Vec3(cont_centers_tib_transf[pairs_list[i][0]-1][0], cont_centers_tib_transf[pairs_list[i][0]-1][1], cont_centers_tib_transf[pairs_list[i][0]-1][2]);
//        Vec3 arg3 = Vec3(tibplanes[pairs_list[i][0]-1][0], tibplanes[pairs_list[i][0]-1][1], tibplanes[pairs_list[i][0]-1][2]);
//        Real d = (dot(arg1, arg2) + femplanes[pairs_list[i][1]-1][3]) / (-dot(arg1, arg3));
//        std::cout << "arg1=" << arg1 << std::endl;
//        std::cout << "arg2=" << arg2 << std::endl;
//        std::cout << "arg3=" << arg3 << std::endl;
//        Vec3 prod_aux = Vec3(tibplanes[pairs_list[i][0]-1][0] * d, tibplanes[pairs_list[i][0]-1][1] * d, tibplanes[pairs_list[i][0]-1][2] * d);
//        Vec3 Pt = cont_centers_tib_transf[pairs_list[i][0]-1] + prod_aux;
//        Vec3 d_mult = cont_centers_Fem[pairs_list[i][1]-1] - Pt;
//        Real distaux = radfem[pairs_list[i][1]-1] -
//            sqrt(pow(d_mult[0], 2.0) + pow(d_mult[1], 2.0) + pow(d_mult[2], 2.0) + 1e-8);
//        
//        std::cout << "cont_centers_Fem.size=" << cont_centers_Fem.size() << std::endl;
//        std::cout << "femplanes.size=" << femplanes.size() << std::endl;
//        std::cout << "pairs_list[i] =" << pairs_list[i][0] << "," << pairs_list[i][1] << std::endl;
//        std::cout << "tibplanes[pairs_list[i][0]-1]=" << tibplanes[pairs_list[i][0] - 1] << std::endl;
//        std::cout << "femplanes[pairs_list[i][1]-1]=" << femplanes[pairs_list[i][1] - 1] << std::endl;
//        std::cout << "d=" << d << std::endl;
//        std::cout << "prod_aux" << prod_aux << std::endl;
//        std::cout << "Pt=" << Pt << std::endl;
//        std::cout << "d_mult=" << d_mult << std::endl;
//        std::cout << "distaux=" << distaux << std::endl;
//        mlist[i]= CheckContact(distaux);
//        std::cout << "mlist[i]=" << mlist[i] << std::endl;
//    }
//    std::cout << "mlist=" << mlist << std::endl;
//    return mlist;
//}

//Vector_<Real> GenerateMultList_spheres(std::vector<std::vector<int>> pairs_list, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf, Vector_<Real> conFem_r, Vector_<Real> conTib_r) {
//    Vector_<Real> mlist(pairs_list.size(),0.0);
//   
//    // method of spheres to generate multipliers
//    for (int i = 0; i < pairs_list.size(); i++) {
//        Vec3 d = cont_centers_Fem[pairs_list[i][1] - 1] - cont_centers_tib_transf[pairs_list[i][0] - 1];
//        /* std::cout << "d=" << d << std::endl;
//            std::cout << "operation=" << pow(pow(d[0], 2.0) + pow(d[1], 2.0) + pow(d[2], 2.0) + 1.0e-8, 0.5) << std::endl;*/
//        Real distaux = conFem_r[pairs_list[i][1] - 1] + conTib_r[pairs_list[i][0] - 1] - pow(pow(d[0], 2.0) + pow(d[1], 2.0) + pow(d[2], 2.0) + 1.0e-8, 0.5);
//        mlist[i] = CheckContact(distaux);
//    }
//
//    /*std::cout << "mlist" << std::endl; 
//    std::cout << mlist << std::endl;*/
//    return mlist;
//}

void CalculateIntersection_cylinders(Vec3  cont_centers_tib, Vec3 cont_centers_fem, Vec3& d) {
    d = cont_centers_fem - cont_centers_tib;
}

void CalculateIntersection_spheres(Vector_<Vec3> fem_Points, Vector_<Vec3> tib_Points, Vec3 &d) {
    Vec3 Ctib(0.0);
    Vec3 Cfem(0.0);
    // mean point of the tib face
    for (int i = 0; i < 3; i++) {
        Ctib = mean (Vec3(tib_Points[0][i], tib_Points[1][i], tib_Points[2][i]));
    }
    
    // mean point of the fem face
    for (int i = 0; i < 3; i++) {
        Cfem[i] = mean(Vec3(fem_Points[0][i], fem_Points[1][i], fem_Points[2][i]));
    }

    d = Cfem - Ctib;

}

void CalculateIntersection(Vector_<Vec3> sphere_Points, Vec3& d) {
    Vec3 Csphere(0.0);
    // mean point of the tib face
    for (int i = 0; i < 3; i++) {
        Csphere[i] = mean(Vec3(sphere_Points[0][i], sphere_Points[1][i], sphere_Points[2][i]));
        std::cout << Csphere << std::endl;
    }

    d = Csphere;

}

void CalculateMaximumPenetration(Vector_<Vec3> d_v, Vector_<Vec3> nt_v, Real &maxpen, Vec3 &nt_l) {
    Vector_<Real> proj(d_v.size());
    //Vector_<Vec3> dist_v(d_v.size());
    Vector_<Real> pen(d_v.size());
    /*std::cout << proj.size() << std::endl;
    std::cout << mults.size() << std::endl;
    std::cout << mults << std::endl;*/
    for (int i = 0; i < d_v.size(); i++) {
        //proj[i] = dot(d_v[i], nt_v[i]); // projection of distance between tibial and femoral faces to the normal of tibial face
        proj[i] = dot(d_v[i], Vec3(0, 1, 0));
        pen[i] = -proj[i]; // pen is for penetration
        //pen[i] = -proj[i]; // pen is for penetration
    }
    Real k = 1e4;
    std::cout << "proj[0]" << proj[0] << std::endl;
    std::cout << "pen[0]" << pen[0] << std::endl;
    std::cout << "nt_v[0]" << nt_v[0] << std::endl;
    std::cout << "d_v[0]=" << d_v[0] << std::endl;
    /*std::cout << "proj" << proj << std::endl;
    std::cout << "pen" << pen << std::endl;
    std::cout << "k*pen" << k*pen << std::endl;
    std::cout << "exp(k*pen))" << exp(k * pen) << std::endl;*/

    //maxpen = (log(sum(exp(k*pen))+1e-16) / k); //version logSum
    maxpen = (log((1.0/pen.size())*sum(exp(k * pen))+1e-16) / k); //version mellowmax
    //maxpen = max(pen); // version nosmooth
    
    //Real auxval(0.0);
    //for (int i = 0; i < pen.size(); i++) {
    //    auxval = auxval + pow(pen[i], 10.0);
    //}
    //Real auxval2 = pow(auxval, 1.0 / 10);
    //maxpen = sum(pen*(0.5 * tanh(k * (pen - auxval2 + 1e-5)) + 0.5)); // version p norm

    std::cout << "sum(exp(k * pen))" << sum(exp(k * pen)) << std::endl;
    std::cout << "maxpen=" << maxpen << std::endl;
    std::cout << "1.0 / pen.size()" << 1.0/pen.size() << std::endl;
    std::cout << "true max pen=" << max(pen) << std::endl;
    //std::cout << "nt_v=" << nt_v << std::endl;
    //maxpen = max(pen);
    nt_l = Vec3(0,1, 0);
}

Real CalculatePressure(Real poisson, Real E, Real d, Real h) {
    Real k = 5e5; // initially 5e5
    Real pen = d;
   
    Real p_init = ((1 - poisson)*E / ((1 + poisson)*(1 - 2 * poisson)))*pen / h;
    Real s = 0.5 + 0.5 * tanh(k * pen);


    //Real p = p_init*(1 + tanh(k*pen)) / 2;
    Real f = 1.7524e5 * (pow(pen, 3.0)) + 5.8571e5 *pow(pen,2.0)+ 5.3548e5 * pen;
    Real p = s * p_init + (1 - s) * f;

    std::cout << "pen=" << pen << std::endl;
    std::cout << "p=" << p << std::endl;

    return p;
}

void CalculateForceCompartment(std::vector<Vec3> spherePoints, std::vector<std::vector<int>> facesSphere, Vector_<Vec3> spherePoints_transf, std::vector<std::vector<int>> pairs_list, Vec3& SumForces, Vec3& SumMoments, Real poisson, Real E, Real h, Vec3 originSphere_G, Vec3 sphere_trans, Vec3 sphere_vel, Vec3 sphere_rot, Vector_<Real> At, Vector_<Vec4> sphereplanes, Vector_<Vec3> cont_centers_Sphere, Vector_<Vec3> cont_centers_sphere_transf) {
    Vector_<Vec3> d(pairs_list.size());
    Vector_<Vec3> nt(pairs_list.size());
    Vector_<Vec3> force_s_l(0);
    Vector_<Real> face_s(0);
    Vector_<Vec3> mom_O(0);

    int k = 1; //number of contacting elements in femur
    int l = 1; //count number of elements contacting element k of femur

    //std::cout << "facesFem.size()" << facesFem.size() << std::endl;
    //for (int i = 0; i < facesFem.size(); i++) {
    //    std::cout << "faces_fem[" << i << "]=" << facesFem.at(i).at(0) << " " << facesFem.at(i).at(1) << facesFem.at(i).at(2) << std::endl;
    //}
    std::cout << "spherePoints" << spherePoints[0] << std::endl;
    for (int i = 0; i < sphereplanes.size(); i++) {
        std::cout << "sphereplanes[pairs_list[i][0]]=" << sphereplanes[pairs_list[i][1]-1] << std::endl;
    }
    if (strcmp(method_pen, "mesh") == 0) {
        for (int i = 0; i < pairs_list.size(); i++) {
            Vec3 d_aux;
            ///////// check Csphere_aux

     /*       if (strcmp(multiplier_method, "cylinders") == 0) {
                std::cout << multiplier_method << std::endl;
                CalculateIntersection_cylinders(cont_centers_tib_transf[pairs_list[i][0] - 1], cont_centers_Fem[pairs_list[i][1] - 1], d_aux);
            }
            else if (strcmp(multiplier_method, "spheres") == 0) {
                std::vector<int> facesFem_ind = facesFem[pairs_list[i][1] - 1];
                Vector_<Vec3> fem_Points_ind(3);
                for (int j = 0; j < 3; j++) {
                    fem_Points_ind[j] = femPoints[facesFem_ind[j] - 1];
                }*/

            std::vector<int> facesSphere_ind = facesSphere[pairs_list[i][1] - 1];
            Vector_<Vec3> sphere_Points_ind(3);
            for (int j = 0; j < 3; j++) {
                sphere_Points_ind[j] = spherePoints_transf[facesSphere_ind[j] - 1];
            }
            std::cout << "facesSphere_ind=" << facesSphere_ind[0] << "-" << facesSphere_ind[1] << "-" << facesSphere_ind[2] << std::endl;
            std::cout << "sphere_Points_ind[0]" << sphere_Points_ind[0] << std::endl;
            std::cout << "spherePoints_transf=" << spherePoints_transf[0] << std::endl;
            CalculateIntersection(sphere_Points_ind, d_aux); // first element in pairs is tibia, second femur
        //}
            std::cout << "sphere_Points_ind" << sphere_Points_ind << std::endl;
            std::cout << "d_aux=" << d_aux << std::endl;
            d[i] = d_aux;
            nt[i] = Vec3(sphereplanes[pairs_list[i][1] - 1][0], sphereplanes[pairs_list[i][1] - 1][1], sphereplanes[pairs_list[i][1] - 1][2]);
            /*std::cout << multipliers << std::endl;*/
            /*std::cout << d << std::endl;*/

            if (i > 0) {
                if ((pairs_list[i][0] == pairs_list[i - 1][0]) && (i < pairs_list.size() - 1)) {
                    l = l + 1;
                }
                else {
                    Vector_<Vec3> d_aux_list(l - 1);
                    Vector_<Vec3> nt_aux_list(l - 1);
                    Real maxpen;
                    Vec3 nt_l;
                    for (int j = 0; j < l - 1; j++) {
                        d_aux_list[j] = d[i - l + j + 1];
                        nt_aux_list[j] = nt[i - l + j + 1];
                    }
                    //std::cout << "d_aux_list=" << d_aux_list << std::endl;
                    std::cout << "i=" << i << std::endl;
                    //std::cout << "d=" << d << std::endl;
                    CalculateMaximumPenetration(d_aux_list, nt_aux_list, maxpen, nt_l);
                    std::cout << "maxpen=" << maxpen << std::endl;
                    std::cout << At[pairs_list[i][1] - 1] << std::endl;
                    std::cout << "nt_l=" << nt_l << std::endl;

                    /*E = 1e9;
                    poisson = 0.45;
                    maxpen = -0.01;*/

                    force_s_l.resizeKeep(k);
                    std::cout << mod << std::endl;

                    if (strcmp(mod, "HC")==0){
                        const Real vt = 0.2;
                        const Real us = 0.8;
                        const Real ud = 0.8;
                        const Real uv = 0.5;
                        const Real bd = 300;
                        const Real bv = 50;
                        // Calculate the Hertz force.
                        const Real k2 = (1. / 2.) * pow(E, 2. / 3.); //E should be "stiffness" here
                        const Real fh_pos = (4. / 3.) * k2 * sqrt(0.041 * k2) *
                            pow(sqrt(maxpen * maxpen + 1e-16), 3. / 2.);
                        const Real fh_smooth = fh_pos * (1. / 2. + (1. / 2.) * tanh(bd * maxpen));

                        Real vpen = -sphere_vel[1];
                        std::cout << sphere_vel << std::endl;

                        const Real c = 1.5; //dissipation
                        const Real fhc_pos = fh_smooth * (1. + (3. / 2.) * c * vpen);
                        const Real fhc_smooth = fhc_pos * (1. / 2. + (1. / 2.) * tanh(bv * (vpen + (2. / (3. * c)))));

                        force_s_l[k - 1] = fhc_smooth * nt_l;

                    }
                    else if (strcmp(mod, "EF") == 0) {
                            Real p = CalculatePressure(poisson, E, maxpen, h);
                        std::cout << "p=" << p << std::endl;

                        
                        force_s_l[k - 1] = p * At[pairs_list[i][1] - 1] * nt_l;
                    }
                    face_s.resizeKeep(k);
                    face_s[k - 1] = pairs_list[i - 2][1];
                    mom_O.resizeKeep(k);
                    mom_O[k - 1] = cross(cont_centers_sphere_transf[pairs_list[i][1] - 1] - originSphere_G, force_s_l[k - 1]);
                    std::cout << "force_s_l[k - 1])" << force_s_l[k - 1] << std::endl;
                    k = k + 1;
                    l = 2;
                    while (k < pairs_list[i][1]) {
                        k = k + 1;
                    }
                    ///////////////

                }
            }
            else {
                l = l + 1;
            }

        }
       /* for (int i = 0; i < nt.size(); i++) {
            std::cout << "nt=" << nt[i] << std::endl;
        }*/
        std::cout << "force_s_l=" << force_s_l << std::endl;
        Vec3 Sum_Force_G;
        Vec3 Sum_Moments_G;
        Sum_Force_G.setToZero();
        Sum_Moments_G.setToZero();
        for (int i = 0; i < force_s_l.size(); i++) {
            Sum_Force_G = Sum_Force_G + force_s_l[i];
            Sum_Moments_G = Sum_Moments_G + mom_O[i];
        }
        SumForces = Sum_Force_G;
        SumMoments = Sum_Moments_G;
    }
    else if (strcmp(method_pen, "centerdist") == 0) {
        Real maxpen(0.0);
        maxpen = -sphere_trans[1];
        Vec3 SumForces(0.0);
        Vec3 SumMoments(0.0);
        //Real f = CalculateContactForce_centerdist(poisson, E, maxpen, h);

    }
    
 /*   Mat33 Mrottib = R_aux(sphere_trans, sphere_rot);
    
    Mat33 MrottibTrans = Mrottib.transpose();*/
    //SumForces = MrottibTrans*Sum_Force_G;
    //SumMoments = MrottibTrans*Sum_Moments_G;

    
    std::cout << SumForces << std::endl;

}

void ComputeKneeContactForces(Vec3 sphere_trans, Vec3 sphere_vel, Vec3 sphere_rot, Vec3 &SumForces, Vec3 &SumMoments, Real E, Real poisson) {
    //// Read Geometry Information
    std::string root_folder = "C:/Gil/Docencia_UPC/TFGs_TFMs/AlbertMataro/DadesExpEsfera/3D_solid/";
    // Read Points
    std::string filename_spherePoints = root_folder + "spherePoints_" + std::to_string(nfacesSphere) + ".csv";
    std::vector<Vec3> spherePoints = ReadDataDoublex3columns(filename_spherePoints);
    for (int i = 0; i < spherePoints.size(); i++) {
        spherePoints[i] = spherePoints[i] / 1000 - 0.042; // to subtract the radius since the origin was not at the center
    }
   
    std::cout << filename_spherePoints << std::endl;

    // Read Faces
    std::string filename_facesSphere = root_folder + "facesSphere_" + std::to_string(nfacesSphere) + ".csv";
    std::vector<std::vector<int>> facesSphere = ReadDataIntx3columns(filename_facesSphere);
    
    // Read centers
    std::string filename_conSphere = root_folder + "conSphere_" + std::to_string(nfacesSphere) + ".csv";
    std::vector<Vec4> conSphere = ReadDataDoublex4columns(filename_conSphere);
    for (int i = 0; i < conSphere.size(); i++) {
        conSphere[i] = conSphere[i] / 1000 ;
        for (int j = 0; j < 3; j++) {
            conSphere[j] = conSphere[j] - 0.042; // to subtract the radius since the origin was not at the center
        }
    }
  
    //Read pairs
    std::string filename_pairs = root_folder + "pairs_" + std::to_string(nfacesSphere) + ".csv";

    std::vector<std::vector<int>> pairs_list = ReadDataIntx2columns(filename_pairs);
    
    std::cout << filename_pairs << std::endl;
    std::cout << pairs_list[4][0] << " " << pairs_list[4][1] << std::endl;  

    Vector_<Vec3> cont_centers_Sphere(conSphere.size(), Vec3(0));
    Vector_<Real> conSphere_r(conSphere.size());
    for (int i = 0; i < conSphere.size(); i++) {
        cont_centers_Sphere[i]=Vec3(conSphere[i][0],conSphere[i][1], conSphere[i][2]);
        conSphere_r[i] = conSphere[i][3];
    }

    std::cout << facesSphere[0][0] << facesSphere[0][1] << facesSphere[0][2] << std::endl;
    std::cout << facesSphere[1][0] << facesSphere[1][1] << facesSphere[1][2] << std::endl;

    //test with numerical values
    std::cout << "sphere_trans=" << sphere_trans << std::endl;
    std::cout << "sphere_rot=" << sphere_rot << std::endl;
    //// Apply transformations
    // Get matrix transformation of sphere with respect to ground
    Mat44 Mtransf_sphere = ftransf_function(sphere_trans, sphere_rot);
    std::cout << "Mtransf_sphere=" << Mtransf_sphere << std::endl;

    // Apply the transformation to all points of the sphere
    Vec4 spherePoints_transf_aux(1);
    Vector_<Vec3> spherePoints_transf(size(spherePoints));
    Vec4 cont_centers_sphere_aux(1);
    Vector_<Vec3> cont_centers_sphere_transf(conSphere.size());
    for (int i = 0; i < size(spherePoints); i++) {
        spherePoints_transf_aux = Mtransf_sphere*Vec4(spherePoints[i][0], spherePoints[i][1], spherePoints[i][2], 1.0);

        for (int j = 0; j < 3; j++) {
            spherePoints_transf[i][j] = spherePoints_transf_aux[j];
        }
    }
    
    std::cout << "spherePoints[0]=" << spherePoints[0] << std::endl;
    std::cout << "spherePoints_transf[0]=" << spherePoints_transf[0] << std::endl;
    std::cout << conSphere.size() << std::endl;
    for (int i = 0; i < conSphere.size(); i++) { // is the size of conTib correct?
        cont_centers_sphere_aux = Mtransf_sphere * Vec4(conSphere[i][0], conSphere[i][1], conSphere[i][2], 1.0);

        for (int j = 0; j < 3; j++) {
            cont_centers_sphere_transf[i][j] = cont_centers_sphere_aux[j];
        }
    }

    // Calculate origin of the sphere...
    Vec4 originSphere_G4 = Mtransf_sphere*Vec4(0, 0, 0, 1);
    Vec3 originSphere_G = Vec3(originSphere_G4[0], originSphere_G4[1], originSphere_G4[2]);
    // 
    std::cout << "originSphere_G=" << originSphere_G << std::endl;
    // Calculate multipliers for all pairs at this instant, and areas of tibia faces
    //Vector_<Real> multipliers(pairs_list.size(), 0.0);
    Vector_<Real> At1(conSphere_r.size(), 0.0);
    Vector_<Vec4> sphereplanes(conSphere_r.size(), Vec4(0.0));
    for (int i = 0; i < conSphere_r.size(); i++) {
        Vec3 edge1_t = spherePoints_transf[facesSphere[i][1] - 1] - spherePoints_transf[facesSphere[i][0] - 1];
        Vec3 edge2_t = spherePoints_transf[facesSphere[i][2] - 1] - spherePoints_transf[facesSphere[i][1] - 1];
        //// unit normal vector to triangle of tibia1
        Vec3 nt_i = cross(edge1_t, edge2_t);
        nt_i = nt_i.normalize();
        Real Dplanesphere = -dot(nt_i, cont_centers_sphere_transf[i]);
        sphereplanes[i][0] = nt_i[0];
        sphereplanes[i][1] = nt_i[1];
        sphereplanes[i][2] = nt_i[2];
        sphereplanes[i][3] = Dplanesphere;
        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
        At1[i] = (1.0 / 2.0) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
    }
    std::cout << "At1=" << At1 << std::endl;
    //if (strcmp(multiplier_method, "cylinders") == 0) {
    //    // method of cylinders to generate multipliers
    //    
    //    // Calculate multipliers for all pairs at this instant
    //    multipliers = GenerateMultList_cylinders(pairs_list, cont_centers_Sphere, cont_centers_sphere_transf, sphereplanes, conSphere_r);
    //   
    //}

    //else if (strcmp(multiplier_method, "spheres") == 0) {
    //   /* for (int i = 0; i < conTib1_r.size(); i++) {
    //        Vec3 edge1_t = tibPoints_transf[facesTib1[i][1] - 1] - tibPoints_transf[facesTib1[i][0] - 1];
    //        Vec3 edge2_t = tibPoints_transf[facesTib1[i][2] - 1] - tibPoints_transf[facesTib1[i][1] - 1];
    //        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
    //        At1[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
    //    }
    //    for (int i = 0; i < conTib2_r.size(); i++) {
    //        Vec3 edge1_t = tibPoints_transf[facesTib2[i][1] - 1] - tibPoints_transf[facesTib2[i][0] - 1];
    //        Vec3 edge2_t = tibPoints_transf[facesTib2[i][2] - 1] - tibPoints_transf[facesTib2[i][1] - 1];
    //        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
    //        At2[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
    //    }*/
    //    multipliers = GenerateMultList_spheres(pairs_list, cont_centers_Sphere, cont_centers_sphere_transf, conSphere_r);
    //    
    //}

    //std::cout << "multipliers1" << multipliers << std::endl;


   /* Real poisson =0.46; // this is input
    Real E = 400 * 1e6;*/ // this is input
    Real h = 0.006;
    SumForces.setToZero();
    SumMoments.setToZero();


    CalculateForceCompartment(spherePoints, facesSphere, spherePoints_transf, pairs_list, SumForces, SumMoments, poisson, E, h, originSphere_G, sphere_trans, sphere_vel, sphere_rot, At1, sphereplanes, cont_centers_Sphere, cont_centers_sphere_transf);
    std::cout << "forces=" << SumForces << " moments=" << SumMoments << std::endl;
 
    //SumForces_vert_Med = SumForces1[1]; // tibial part 1 is medial
    //SumForces_vert_Lat = SumForces2[1]; // tibial part 2 is lateral
    std::cout << "sumforces=" << SumForces << std::endl;

}
// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    // OpenSim model: create components
    /// Model
    OpenSim::Model* model;
    /// Bodies
    OpenSim::Body* sphere;
   
    /// Joints
    OpenSim::CustomJoint* ground_sphere;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    sphere = new OpenSim::Body("sphere", 0.79, Vec3(0, 0, 0), Inertia(5.5742e-04, 5.5742e-04, 5.5742e-04, 0, 0, 0));
    
    /// Joints
    /// Ground-Sphere transform
    SpatialTransform st_ground_sphere;
    st_ground_sphere[0].setCoordinateNames(OpenSim::Array<std::string>("sphere_tilt", 1, 1));
    st_ground_sphere[0].setFunction(new LinearFunction());
    st_ground_sphere[0].setAxis(Vec3(0, 0, 1));
    st_ground_sphere[1].setCoordinateNames(OpenSim::Array<std::string>("sphere_list", 1, 1));
    st_ground_sphere[1].setFunction(new LinearFunction());
    st_ground_sphere[1].setAxis(Vec3(1, 0, 0));
    st_ground_sphere[2].setCoordinateNames(OpenSim::Array<std::string>("sphere_rotation", 1, 1));
    st_ground_sphere[2].setFunction(new LinearFunction());
    st_ground_sphere[2].setAxis(Vec3(0, 1, 0));
    st_ground_sphere[3].setCoordinateNames(OpenSim::Array<std::string>("sphere_tx", 1, 1));
    st_ground_sphere[3].setFunction(new LinearFunction());
    st_ground_sphere[3].setAxis(Vec3(1, 0, 0));
    st_ground_sphere[4].setCoordinateNames(OpenSim::Array<std::string>("sphere_ty", 1, 1));
    st_ground_sphere[4].setFunction(new LinearFunction());
    st_ground_sphere[4].setAxis(Vec3(0, 1, 0));
    st_ground_sphere[5].setCoordinateNames(OpenSim::Array<std::string>("sphere_tz", 1, 1));
    st_ground_sphere[5].setFunction(new LinearFunction());
    st_ground_sphere[5].setAxis(Vec3(0, 0, 1));

    /// Joint specifications
    ground_sphere = new CustomJoint("ground_sphere", model->getGround(), Vec3(0), Vec3(0), *sphere, Vec3(0), Vec3(0), st_ground_sphere);
    
    /// Add bodies and joints to model
    model->addBody(sphere);		        model->addJoint(ground_sphere);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    //CALL STATE INDEX MAPPING FUNCTION TO ACCOUNT FOR OPENSIM VS SIMBODY STATE ORDERS
    Array<std::string> stateVars = model->getStateVariableNames();					//Assign string array with the state variable names
    std::cout << stateVars << std::endl;
    std::unordered_map<std::string, int> mapping = createSystemYIndexMap(*model);    //Call function
    for (int i = 0; i < mapping.size(); ++i) std::cout << mapping[stateVars[i]] << " " << stateVars[i] << " " << i << " OpenSim" << std::endl; //Loop through each state name and print to the cmd window the corresponding Simbody index and the name

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);
    std::vector<T> p(arg[2], arg[2] + NP);

    // States and controls
    T ua[NU]; /// joint accelerations (Qdotdots) - controls
    T up[NP]; /// choose model - parameters
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
    
    

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
   
    /// Controls
    for (int i = 0; i < NU; ++i) ua[i] = u[i];
   
    /// Parameters
    for (int i = 0; i < NP; ++i) up[i] = p[i];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    
    std::cout << "QsUs=" << QsUs[9] << std::endl;

    model->realizeVelocity(*state);

    std::cout << "ty=" << model->getStateVariableValue(*state, "ground_sphere/sphere_ty/value") << std::endl;

    Vec3 sphere_trans = Vec3(model->getStateVariableValue(*state, "ground_sphere/sphere_tx/value"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_ty/value"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_tz/value"));
    Vec3 sphere_rot = Vec3(model->getStateVariableValue(*state, "ground_sphere/sphere_tilt/value"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_list/value"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_rotation/value"));
    
    Vec3 sphere_vel = Vec3(model->getStateVariableValue(*state, "ground_sphere/sphere_tx/speed"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_ty/speed"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_tz/speed"));
    Vec3 sphere_rot_vel = Vec3(model->getStateVariableValue(*state, "ground_sphere/sphere_tilt/speed"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_list/speed"),
        model->getStateVariableValue(*state, "ground_sphere/sphere_rotation/speed"));

    std::cout << "sphere_vel " << sphere_vel << std::endl;
    std::cout << "sphere_rot_vel " << sphere_rot_vel << std::endl;
    std::cout << "QsUs " << QsUs << std::endl;

    // Compute Resulting contact wrench at the center of the sphere
    Vec3 Cont_SumForces;
    Vec3 Cont_SumMoments;
    Cont_SumForces.setToZero();
    Cont_SumMoments.setToZero();
    Real SumForces_vert_Lat;
    Real SumForces_vert_Med;
    
    Real E = up[0];
    Real poisson = up[1];

    ComputeKneeContactForces(sphere_trans, sphere_vel, sphere_rot, Cont_SumForces, Cont_SumMoments, E, poisson);
    Vec3 SphereCont_SumForces_onSphere_inG = Cont_SumForces;
    Vec3 SphereCont_SumMoments_onSphere_inG = Cont_SumMoments;

   /* Vec3 SphereCont_SumForces_onSphere_inG=sphere->expressVectorInGround(*state, KneeCont_SumForces_onTibialTray_inTibialTrayFrame);
    Vec3 KneeCont_SumMoments_onTibialTray_inG = tibial_tray->expressVectorInGround(*state, KneeCont_SumMoments_onTibialTray_inTibialTrayFrame);*/

   /* Vec3 KneeCont_SumForces_onFemoralComp_inG = -KneeCont_SumForces_onTibialTray_inG;
    Vec3 tibial_tray_Center_inG = tibial_tray->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 femoral_comp_Center_inG = femoral_component->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 KneeCont_SumMoments_onFemoralComp_inG = -KneeCont_SumMoments_onTibialTray_inG + cross(tibial_tray_Center_inG - femoral_comp_Center_inG, KneeCont_SumForces_onFemoralComp_inG);*/

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
    /// Add weights to appliedBodyForces
    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
        model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
    }
    std::cout << "sphere_trans= " << sphere_trans << std::endl;
    std::cout << "SphereCont_SumForces_onSphere_inG= " << SphereCont_SumForces_onSphere_inG << std::endl;
    std::cout << "weight=" << model->getBodySet().get("sphere").getMass() * gravity << std::endl;
    std::cout << "appliedBodyForces=" << appliedBodyForces << std::endl;
    /// Add sphere contact forces to appliedBodyForces
    model->getMatterSubsystem().addInStationForce(*state, sphere->getMobilizedBodyIndex(), Vec3(0, 0, 0), SphereCont_SumForces_onSphere_inG, appliedBodyForces);
    model->getMatterSubsystem().addInBodyTorque(*state, sphere->getMobilizedBodyIndex(), SphereCont_SumMoments_onSphere_inG, appliedBodyForces);
    std::cout << "appliedBodyForces=" << appliedBodyForces << std::endl;

   /* Vector_<SpatialVec> aaa;
    aaa.resize(nbodies);
    aaa.setToZero();
    model->getMatterSubsystem().addInStationForce(*state, femoral_component->getMobilizedBodyIndex(), Vec3(0, 0, 0), KneeCont_SumForces_G_onFemoralComp, aaa);
    model->getMatterSubsystem().addInBodyTorque(*state, femoral_component->getMobilizedBodyIndex(), KneeCont_SumMoments_G_onFemoralComp, aaa);
    Vector_<SpatialVec> bbb;
    bbb.resize(nbodies);
    bbb.setToZero();
    Vec3 TibialTrayC_onFemComp = tibial_tray->findStationLocationInAnotherFrame(*state, Vec3(0, 0, 0), *femoral_component);
    model->getMatterSubsystem().addInStationForce(*state, femoral_component->getMobilizedBodyIndex(), TibialTrayC_onFemComp, KneeCont_SumForces_G_onFemoralComp, bbb);
    Vec3 mKneeCont_SumMoments_G_onTibialTray = -KneeCont_SumMoments_G_onTibialTray;
    model->getMatterSubsystem().addInBodyTorque(*state, femoral_component->getMobilizedBodyIndex(), mKneeCont_SumMoments_G_onTibialTray, bbb);
    
    std::cout << "aaa=" << aaa << std::endl;
    std::cout << "bbb=" << bbb << std::endl;*/

    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    ///  Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Compute contact torques about the ground frame origin
 
    // Extract results
    /// Residual forces
    for (int i = 0; i < ndof; ++i) {
        res[0][i] = value<T>(residualMobilityForces[i]);
    }

    /// Contact forces
    for (int i = 0; i < 3; ++i) {
        res[0][i + ndof] = value<T>(SphereCont_SumForces_onSphere_inG[i]);      /// Sphere contact forces
    }
    for (int i = 0; i < 3; ++i) {
        res[0][i + ndof + 3] = value<T>(SphereCont_SumMoments_onSphere_inG[i]); /// Sphere contact moments
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
    Recorder u[NU];
    Recorder p[NP];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;
    for (int i = 0; i < NP; ++i) p[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u,p };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}

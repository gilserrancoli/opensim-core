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
constexpr int ndof = 6;        // # degrees of freedom (excluding locked)
constexpr int NX = ndof;      // # states
//constexpr int NU = ndof;        // # controls

//constexpr int numpairs = 932; //499 is right cycle with 5 mm radius sphere threshold, 932 is with 10 mm threshold
constexpr int nfacesTib = 49; //before 49
constexpr int nfacesFem = 188;
constexpr const char radForPairs[] = "05"; // 1 is 1 cm, 05 is 0.5 cm
constexpr char* multiplier_method = "cylinders"; // multiplier method: "cylinders" or "spheres"

constexpr int NR = 2 + nfacesTib;    // # residual torques + # GRFs + # GRMs

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

Mat33 R_aux(Vec3 knee_trans, Vec3 knee_rot) {
    Recorder psi = knee_rot[0];
    Recorder theta = knee_rot[1];
    Recorder phi = knee_rot[2];

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

Mat44 ftransf_function(Vec3 knee_trans, Vec3 knee_rot) {

    Mat33 R = R_aux(knee_trans, knee_rot);
    
    Vec3 aux = R*(-knee_trans);
    
    Mat44 Rtrans0042(1);
    Rtrans0042.set(1, 3, 0.042);
        
    Mat44 Rtranstib_4x4(1);
    Rtranstib_4x4.set(0, 3, aux[0]);
    Rtranstib_4x4.set(1, 3, aux[1]);
    Rtranstib_4x4.set(2, 3, aux[2]);

    std::cout << Rtranstib_4x4 << std::endl;

    Mat44 R_rot(1);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_rot.set(i, j, R[i][j]);
        }
    }

    Mat44 R_tib = Rtranstib_4x4*Rtrans0042*R_rot;
    return R_tib;
}

Real CheckContact(Real overlap) {
    Real k = 1e3; // k subject to change
    Real multiplier = (tanh(k * overlap) + 1.0) / 2.0;
    /*std::cout << multiplier << std::endl;*/
    return multiplier;
}

Vector_<Real> GenerateMultList_cylinders(std::vector<std::vector<int>> pairs_list, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf, Vector_<Vec4> tibplanes, Vector_<Vec4> femplanes, Vector_<Real> radfem) {
    std::cout << pairs_list.size() << std::endl;
    Vector_<Real> mlist(pairs_list.size(), 0.0);
    for (int i = 0; i < pairs_list.size(); i++) {
        Vec3 arg1 = Vec3(femplanes[pairs_list[i][1]-1][0], femplanes[pairs_list[i][1]-1][1], femplanes[pairs_list[i][1]-1][2]);
        Vec3 arg2 = Vec3(cont_centers_tib_transf[pairs_list[i][0]-1][0], cont_centers_tib_transf[pairs_list[i][0]-1][1], cont_centers_tib_transf[pairs_list[i][0]-1][2]);
        Vec3 arg3 = Vec3(tibplanes[pairs_list[i][0]-1][0], tibplanes[pairs_list[i][0]-1][1], tibplanes[pairs_list[i][0]-1][2]);
        Real d = (dot(arg1, arg2) + femplanes[pairs_list[i][1]-1][3]) / (-dot(arg1, arg3));
        std::cout << "arg1=" << arg1 << std::endl;
        std::cout << "arg2=" << arg2 << std::endl;
        std::cout << "arg3=" << arg3 << std::endl;
        Vec3 prod_aux = Vec3(tibplanes[pairs_list[i][0]-1][0] * d, tibplanes[pairs_list[i][0]-1][1] * d, tibplanes[pairs_list[i][0]-1][2] * d);
        Vec3 Pt = cont_centers_tib_transf[pairs_list[i][0]-1] + prod_aux;
        Vec3 d_mult = cont_centers_Fem[pairs_list[i][1]-1] - Pt;
        Real distaux = radfem[pairs_list[i][1]-1] -
            sqrt(pow(d_mult[0], 2.0) + pow(d_mult[1], 2.0) + pow(d_mult[2], 2.0) + 1e-8);
        
        std::cout << "cont_centers_Fem.size=" << cont_centers_Fem.size() << std::endl;
        std::cout << "femplanes.size=" << femplanes.size() << std::endl;
        std::cout << "pairs_list[i] =" << pairs_list[i][0] << "," << pairs_list[i][1] << std::endl;
        std::cout << "tibplanes[pairs_list[i][0]-1]=" << tibplanes[pairs_list[i][0] - 1] << std::endl;
        std::cout << "femplanes[pairs_list[i][1]-1]=" << femplanes[pairs_list[i][1] - 1] << std::endl;
        std::cout << "d=" << d << std::endl;
        std::cout << "prod_aux" << prod_aux << std::endl;
        std::cout << "Pt=" << Pt << std::endl;
        std::cout << "d_mult=" << d_mult << std::endl;
        std::cout << "distaux=" << distaux << std::endl;
        mlist[i]= CheckContact(distaux);
        std::cout << "mlist[i]=" << mlist[i] << std::endl;
    }
    std::cout << "mlist=" << mlist << std::endl;
    return mlist;
}

Vector_<Real> GenerateMultList_spheres(std::vector<std::vector<int>> pairs_list, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf, Vector_<Real> conFem_r, Vector_<Real> conTib_r) {
    Vector_<Real> mlist(pairs_list.size(),0.0);
   
    // method of spheres to generate multipliers
    for (int i = 0; i < pairs_list.size(); i++) {
        Vec3 d = cont_centers_Fem[pairs_list[i][1] - 1] - cont_centers_tib_transf[pairs_list[i][0] - 1];
        /* std::cout << "d=" << d << std::endl;
            std::cout << "operation=" << pow(pow(d[0], 2.0) + pow(d[1], 2.0) + pow(d[2], 2.0) + 1.0e-8, 0.5) << std::endl;*/
        Real distaux = conFem_r[pairs_list[i][1] - 1] + conTib_r[pairs_list[i][0] - 1] - pow(pow(d[0], 2.0) + pow(d[1], 2.0) + pow(d[2], 2.0) + 1.0e-8, 0.5);
        mlist[i] = CheckContact(distaux);
    }

    /*std::cout << "mlist" << std::endl; 
    std::cout << mlist << std::endl;*/
    return mlist;
}

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

void CalculateMaximumPenetration(Vector_<Vec3> d_v, Vector_<Vec3> nt_v, Real &maxpen, Vec3 &nt_l, Vector_<Real> mults) {
    Vector_<Real> proj(d_v.size());
    //Vector_<Vec3> dist_v(d_v.size());
    Vector_<Real> pen(d_v.size());
    /*std::cout << proj.size() << std::endl;
    std::cout << mults.size() << std::endl;
    std::cout << mults << std::endl;*/
    for (int i = 0; i < d_v.size(); i++) {
        proj[i] = dot(d_v[i], nt_v[i]); // projection of distance between tibial and femoral faces to the normal of tibial face
        pen[i] = -proj[i]*mults[i]; // pen is for penetration
        //pen[i] = -proj[i]; // pen is for penetration
    }
    Real k = 1e4;
    std::cout << "proj" << proj << std::endl;
    std::cout << "pen" << pen << std::endl;
    std::cout << "mults=" << mults << std::endl;
    std::cout << "sum mults=" << sum(mults) << std::endl;
    std::cout << "k*pen" << k*pen << std::endl;
    std::cout << "exp(k*pen))" << exp(k * pen) << std::endl;

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
    std::cout << "nt_v=" << nt_v << std::endl;
    //maxpen = max(pen);
    nt_l = nt_v[0];
}

Real CalculatePressure(Real poisson, Real E, Real d, Real h) {
    //Real k = 5e5;
    Real pen = d;

    Real p_init = ((1 - poisson)*E / ((1 + poisson)*(1 - 2 * poisson)))*pen / h;

    // Real p = p_init*(1 + tanh(k*pen)) / 2; //smoothing using a tanh curve
    Real fpen_negative = 10000 * pen;
    Real k = 8e7;
    Real stiffness = ((1 - poisson) * E / ((1 + poisson) * (1 - 2 * poisson))) / h;;
    Real trans = -1e-7; // transition point from constant slope to start to change the slope 
    Real midpoint = trans / 2;
    Real alfa = (1000 + stiffness) / 2;
    Real beta = ((stiffness - 1000) / 2) / tanh(-k * midpoint);

    Real p;
    if (pen <= trans) {
        p = 1000 * pen;
    }
    else if ((pen > trans)& (pen < 0.0)) {
        p = pen * (alfa + beta) - (beta * log(tanh(k * pen + 4) + 1)) / k + (beta * log(tanh(4) + 1)) / k;
    }
    else if (pen >= 0.0) {
        p = p_init;
    }

    std::cout << "pen" << pen << std::endl;
    std::cout << "p" << p << std::endl;

    return p;
}

void CalculateForceCompartment(std::vector<Vec3> femPoints, std::vector<std::vector<int>> facesFem, Vector_<Vec3> tibPoints_transf, std::vector<std::vector<int>> facesTib, std::vector<std::vector<int>> pairs_list, Vec3& SumForces, Vec3& SumMoments, Real poisson, Real E, Real h, Vec3 originTib_G, Vector_<Real> multipliers, Vec3 knee_trans, Vec3 knee_rot, Vector_<Real> At, Vector_<Vec4> tibplanes, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf, Vector& p_vec) {
    Vector_<Vec3> d(pairs_list.size());
    Vector_<Vec3> nt(pairs_list.size());
    Vector_<Vec3> force_s_l(0);
    Vector_<Real> face_s(0);
    Vector_<Vec3> mom_O(0);

    Vector maxpen_v(size(facesTib));
    maxpen_v.setToZero();

    int k = 1; //number of contacting elements in femur
    int l = 1; //count number of elements contacting element k of femur

    //std::cout << "facesFem.size()" << facesFem.size() << std::endl;
    //for (int i = 0; i < facesFem.size(); i++) {
    //    std::cout << "faces_fem[" << i << "]=" << facesFem.at(i).at(0) << " " << facesFem.at(i).at(1) << facesFem.at(i).at(2) << std::endl;
    //}
    std::cout << "femPoints" << std::endl;
    for (int i = 0; i < tibplanes.size(); i++) {
        std::cout << "tibplanes[pairs_list[i][0]]=" << tibplanes[pairs_list[i][0]] << std::endl;
    }

    for (int i = 0; i < size(femPoints); i++) {
        std::cout << femPoints[i] << std::endl;
    }

    for (int i = 0; i < pairs_list.size(); i++) {
        Vec3 d_aux;
        ///////// check Ctib_aux

        if (strcmp(multiplier_method, "cylinders") == 0) {
            std::cout << multiplier_method << std::endl;
            CalculateIntersection_cylinders(cont_centers_tib_transf[pairs_list[i][0] - 1], cont_centers_Fem[pairs_list[i][1] - 1], d_aux);
        }
        else if (strcmp(multiplier_method, "spheres") == 0) {
            std::vector<int> facesFem_ind = facesFem[pairs_list[i][1] - 1];
            Vector_<Vec3> fem_Points_ind(3);
            for (int j = 0; j < 3; j++) {
                fem_Points_ind[j] = femPoints[facesFem_ind[j] - 1];
            }

            std::vector<int> facesTib_ind = facesTib[pairs_list[i][0] - 1];
            Vector_<Vec3> tib_Points_ind(3);
            for (int j = 0; j < 3; j++) {
                tib_Points_ind[j] = tibPoints_transf[facesTib_ind[j] - 1];
            }
            CalculateIntersection_spheres(fem_Points_ind, tib_Points_ind, d_aux); // first element in pairs is tibia, second femur
        }
        d[i] = d_aux;
        nt[i] = Vec3(tibplanes[pairs_list[i][0]-1][0], tibplanes[pairs_list[i][0]-1][1], tibplanes[pairs_list[i][0]-1][2]);
        


        if (i > 0) {
            if ((pairs_list[i][0] == pairs_list[i - 1][0]) && (i < pairs_list.size() - 1)) {
                l = l + 1;
            }
            else {
                Vector_<Vec3> d_aux_list(l - 1);
                Vector_<Vec3> nt_aux_list(l - 1);
                Vector_<Real> multipliers_list(l - 1);
                Real maxpen;
                Vec3 nt_l;
                for (int j = 0; j < l - 1; j++) {
                    d_aux_list[j] = d[i - l + j + 1];
                    nt_aux_list[j] = nt[i - l + j + 1];
                    multipliers_list[j] = multipliers[i - l + j + 1];
                }
                CalculateMaximumPenetration(d_aux_list, nt_aux_list, maxpen, nt_l, multipliers_list);
                std::cout << maxpen << std::endl;
                maxpen_v[k - 1] = maxpen;
                //std::cout << At[pairs_list[i][0] - 1] << std::endl;
                //std::cout << "nt_l=" << nt_l << std::endl;
                Real p = CalculatePressure(poisson, E, maxpen, h);
                p_vec[k - 1] = p;

                force_s_l.resizeKeep(k);
                force_s_l[k - 1] = p * At[pairs_list[i][0] - 1] * nt_l;
                face_s.resizeKeep(k);
                face_s[k - 1] = pairs_list[i - 2][0];
                mom_O.resizeKeep(k);
                mom_O[k - 1] = cross(cont_centers_tib_transf[pairs_list[i][0] - 1] - originTib_G, force_s_l[k - 1]);

                k = k + 1;
                l = 2;
                while (k < pairs_list[i][0]) {
                    p_vec[k - 1] = p * 0.0;
                    k = k + 1;
                }
                ///////////////

            }
        }
        else {
            l = l + 1;
        }

    }
    for (int i = 0; i < nt.size(); i++){
        std::cout << "nt=" << nt[i] << std::endl;
    }

    std::cout << "multipliers" << multipliers << std::endl;
    std::cout << "pvec= " << p_vec << std::endl;
    std::cout << "maxpen_v=" << maxpen_v << std::endl;
    
    Vec3 Sum_Force_G;
    Vec3 Sum_Moments_G;
    Sum_Force_G.setToZero();
    Sum_Moments_G.setToZero();
    for (int i = 0; i < force_s_l.size(); i++) {
        Sum_Force_G = Sum_Force_G + force_s_l[i];
        Sum_Moments_G = Sum_Moments_G + mom_O[i];
    }
    
    Mat33 Mrottib = R_aux(knee_trans, knee_rot);
    
    Mat33 MrottibTrans = Mrottib.transpose();
    SumForces = MrottibTrans*Sum_Force_G;
    SumMoments = MrottibTrans*Sum_Moments_G;
    
    std::cout << SumForces << std::endl;

}

void ComputeKneeContactForces(Vec3 knee_trans, Vec3 knee_rot, Vec3 &SumForces, Vec3 &SumMoments, Real &SumForces_vert_Lat, Real &SumForces_vert_Med, Vector& pvec1, Vector& pvec2) {
    //// Read Geometry Information
    std::string root_folder = "C:/Gil/MeshesInAD/WorkshopBarcelona2024/Seminar3-MeshBasedContact/ExperimentalData/";
    // Read Points
    std::string filename_femPoints = root_folder + "femPoints_" + std::to_string(nfacesFem) + ".csv";
    std::vector<Vec3> femPoints = ReadDataDoublex3columns(filename_femPoints);
    std::string filename_tibPoints=root_folder+"tibPoints_" + std::to_string(nfacesTib) + ".csv";
    std::vector<Vec3> tibPoints = ReadDataDoublex3columns(filename_tibPoints);
    // Read Faces
    std::string filename_facesFem = root_folder + "facesFem_" + std::to_string(nfacesFem) + ".csv";
    std::vector<std::vector<int>> facesFem = ReadDataIntx3columns(filename_facesFem);
    std::string filename_facesTib1 = root_folder + "facesTib1_" +std::to_string(nfacesTib) + ".csv";
    std::vector<std::vector<int>> facesTib1 = ReadDataIntx3columns(filename_facesTib1);
    std::string filename_facesTib2 = root_folder + "facesTib2_" + std::to_string(nfacesTib) + ".csv";
    std::vector<std::vector<int>> facesTib2 = ReadDataIntx3columns(filename_facesTib2);
    // Read centers
    std::string filename_conFem = root_folder + "ConFem_" + std::to_string(nfacesFem) + ".csv";
    std::vector<Vec4> conFem = ReadDataDoublex4columns(filename_conFem);
    std::string filename_conTibia1 = root_folder + "ConTib1_" +std::to_string(nfacesTib) + ".csv";
    std::vector<Vec4> conTib1 = ReadDataDoublex4columns(filename_conTibia1);
    std::string filename_conTibia2 = root_folder + "ConTib2_" + std::to_string(nfacesTib) + ".csv";
    std::vector<Vec4> conTib2 = ReadDataDoublex4columns(filename_conTibia2);


    //Read pairs
    std::string filename_pairs1 = root_folder + "pairs1_" + std::to_string(nfacesTib) + "x" + std::to_string(nfacesFem) + "_at" + radForPairs + "cm.csv";
    std::string filename_pairs2 = root_folder + "pairs2_" + std::to_string(nfacesTib) + "x" + std::to_string(nfacesFem) + "_at" + radForPairs + "cm.csv";

    std::cout << filename_pairs1 << std::endl;

    std::vector<std::vector<int>> pairs1_list = ReadDataIntx2columns(filename_pairs1);
    std::vector<std::vector<int>> pairs2_list = ReadDataIntx2columns(filename_pairs2);
    
    std::cout << filename_pairs1 << std::endl;
    std::cout << pairs1_list[4][0] << " " << pairs1_list[4][1] << std::endl;  

    // Sum shift translation to the femur
    for (int i = 0; i < size(femPoints); i++) {
        femPoints[i] = femPoints[i] + Vec3(0.0, 0.042, 0.0);
    }

    Vector_<Vec3> cont_centers_Fem(conFem.size(), Vec3(0));
    Vector_<Real> conFem_r(conFem.size());
    for (int i = 0; i < conFem.size(); i++) {
        cont_centers_Fem[i]=Vec3(conFem[i][0],conFem[i][1], conFem[i][2]) + Vec3(0.0, 0.042, 0.0);
        conFem_r[i] = conFem[i][3];
    }

    std::cout << facesFem[0][0] << facesFem[0][1] << facesFem[0][2] << std::endl;
    std::cout << facesFem[1][0] << facesFem[1][1] << facesFem[1][2] << std::endl;

 
    //test with numerical values
    //knee_trans = Vec3(0.01, 0.03, 0.008);
    //knee_rot = Vec3(0.1, 0.03, 0.025);
    std::cout << "knee_trans=" << knee_trans << std::endl;

    //// Apply transformations
    // Get matrix transformation of tibia with respect to the femur
    Mat44 Mtransf_tib = ftransf_function(knee_trans, knee_rot);
    std::cout << "Mtransf_tib=" << Mtransf_tib << std::endl;

    // Apply the transformation to all points of the tibia
    Vec4 tibPoints_transf_aux(1);
    Vector_<Vec3> tibPoints_transf(size(tibPoints));
    Vec4 cont_centers_tib_aux1(1);
    Vector_<Vec3> cont_centers_tib_transf1(conTib1.size());
    Vector_<Real> conTib1_r(conTib1.size());
    Vec4 cont_centers_tib_aux2(1);
    Vector_<Vec3> cont_centers_tib_transf2(conTib2.size());
    Vector_<Real> conTib2_r(conTib2.size());
    std::cout << "tibPoints=" << tibPoints[0] << std::endl;
    for (int i = 0; i < size(tibPoints); i++) {
        tibPoints_transf_aux = Mtransf_tib*Vec4(tibPoints[i][0], tibPoints[i][1], tibPoints[i][2], 1.0);

        for (int j = 0; j < 3; j++) {
            tibPoints_transf[i][j] = tibPoints_transf_aux[j];
        }
    }
    std::cout << conTib1.size() << std::endl;
    for (int i = 0; i < conTib1.size(); i++) { // is the size of conTib correct?
        cont_centers_tib_aux1 = Mtransf_tib * Vec4(conTib1[i][0], conTib1[i][1], conTib1[i][2], 1.0);
        conTib1_r[i] = conTib1[i][3];

        for (int j = 0; j < 3; j++) {
            cont_centers_tib_transf1[i][j] = cont_centers_tib_aux1[j];
        }
    }
    std::cout << conTib2.size() << std::endl;
    for (int i = 0; i < conTib2.size(); i++) { // is the size of conTib correct?
        cont_centers_tib_aux2 = Mtransf_tib * Vec4(conTib2[i][0], conTib2[i][1], conTib2[i][2], 1.0);
        conTib2_r[i] = conTib2[i][3];

        for (int j = 0; j < 3; j++) {
            cont_centers_tib_transf2[i][j] = cont_centers_tib_aux2[j];
        }
    }

    // Calculate origin of the tibia...
    Vec4 originTib_G4 = Mtransf_tib*Vec4(0, 0, 0, 1);
    Vec3 originTib_G = Vec3(originTib_G4[0], originTib_G4[1], originTib_G4[2]);
    // 

    // Calculate multipliers for all pairs at this instant, and areas of tibia faces
    Vector_<Real> multipliers1(pairs1_list.size(), 0.0);
    Vector_<Real> multipliers2(pairs1_list.size(), 0.0);
    Vector_<Real> At1(conTib1_r.size(), 0.0);
    Vector_<Real> At2(conTib1_r.size(), 0.0);
    Vector_<Vec4> tibplanes1(conTib1_r.size(), Vec4(0.0));
    Vector_<Vec4> tibplanes2(conTib2_r.size(), Vec4(0.0));
    std::cout << "tibPoints_transf=" << tibPoints_transf << std::endl;

    for (int i = 0; i < conTib1_r.size(); i++) {
        std::cout << "facesTib1[i][0]-1=" << facesTib1[i][0] << std::endl;
        std::cout << "facesTib1[i][1]-1=" << facesTib1[i][1] << std::endl;
        std::cout << "facesTib1[i][2]-1=" << facesTib1[i][2] << std::endl;
        std::cout << "tibPoints_transf[facesTib1[i][0] - 1]=" << tibPoints_transf[facesTib1[i][0] - 1] << std::endl;
        std::cout << "tibPoints_transf[facesTib1[i][1] - 1]=" << tibPoints_transf[facesTib1[i][1] - 1] << std::endl;
        std::cout << "tibPoints_transf[facesTib1[i][2] - 1]=" << tibPoints_transf[facesTib1[i][2] - 1] << std::endl;
        Vec3 edge1_t = tibPoints_transf[facesTib1[i][1] - 1] - tibPoints_transf[facesTib1[i][0] - 1];
        Vec3 edge2_t = tibPoints_transf[facesTib1[i][2] - 1] - tibPoints_transf[facesTib1[i][1] - 1];
        //// unit normal vector to triangle of tibia1
        std::cout << "edge1_t=" << edge1_t << std::endl;
        std::cout << "edge2_t=" << edge2_t << std::endl;
        Vec3 nt_i = cross(edge1_t, edge2_t);
        //std::cout << "nt_i=" << nt_i << std::endl;
        nt_i = nt_i.normalize();
        /*std::cout << "nt_i normalized=" << nt_i << std::endl;*/
        Real Dplanetib = -dot(nt_i, cont_centers_tib_transf1[i]);
        tibplanes1[i][0] = nt_i[0];
        tibplanes1[i][1] = nt_i[1];
        tibplanes1[i][2] = nt_i[2];
        tibplanes1[i][3] = Dplanetib;
        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
        At1[i] = (1.0 / 2.0) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
    }
    for (int i = 0; i < conTib2_r.size(); i++) {
        Vec3 edge1_t = tibPoints_transf[facesTib2[i][1] - 1] - tibPoints_transf[facesTib2[i][0] - 1];
        Vec3 edge2_t = tibPoints_transf[facesTib2[i][2] - 1] - tibPoints_transf[facesTib2[i][1] - 1];
        //// unit normal vector to triangle of tibia1
        Vec3 nt_i = cross(edge1_t, edge2_t);
        nt_i = nt_i.normalize();
        Real Dplanetib = -dot(nt_i, cont_centers_tib_transf1[i]);
        tibplanes2[i][0] = nt_i[0];
        tibplanes2[i][1] = nt_i[1];
        tibplanes2[i][2] = nt_i[2];
        tibplanes2[i][3] = Dplanetib;
        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
        At2[i] = (1.0 / 2.0) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
    }
    std::cout << "At1=" << At1 << std::endl;
    std::cout << "At2=" << At2 << std::endl;
    if (strcmp(multiplier_method, "cylinders") == 0) {
        // method of cylinders to generate multipliers

        /// first compute planes
        Vector_<Vec4> femplanes(conFem_r.size(), Vec4(0.0));
        for (int i = 0; i < conFem_r.size(); i++) {
            Vec3 edge1_s = femPoints[facesFem[i][1]-1] - femPoints[facesFem[i][0]-1];
            Vec3 edge2_s = femPoints[facesFem[i][2]-1] - femPoints[facesFem[i][1]-1];
            //// unit normal vector to triangle of femur
            Vec3 ns_i = cross(edge1_s, edge2_s);
            ns_i = ns_i.normalize();
            Real Dplanefem = -dot(ns_i, cont_centers_Fem[i]);
            femplanes[i][0] = ns_i[0];
            femplanes[i][1] = ns_i[1];
            femplanes[i][2] = ns_i[2];
            femplanes[i][3] = Dplanefem;
        }
        
        // Calculate multipliers for all pairs at this instant
        std::cout << "cont_centers_Fem=" << cont_centers_Fem << "/n" << std::endl;
        std::cout << "cont_centers_tib_transf1= " << cont_centers_tib_transf1 << std::endl;
        std::cout << "tibplanes1= " << tibplanes1 << std::endl;
        std::cout << "femplanes= " << femplanes << std::endl;
        multipliers1 = GenerateMultList_cylinders(pairs1_list, cont_centers_Fem, cont_centers_tib_transf1, tibplanes1, femplanes, conFem_r);
        multipliers2 = GenerateMultList_cylinders(pairs2_list, cont_centers_Fem, cont_centers_tib_transf2, tibplanes2, femplanes, conFem_r);
       
    }

    else if (strcmp(multiplier_method, "spheres") == 0) {
       /* for (int i = 0; i < conTib1_r.size(); i++) {
            Vec3 edge1_t = tibPoints_transf[facesTib1[i][1] - 1] - tibPoints_transf[facesTib1[i][0] - 1];
            Vec3 edge2_t = tibPoints_transf[facesTib1[i][2] - 1] - tibPoints_transf[facesTib1[i][1] - 1];
            Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
            At1[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
        }
        for (int i = 0; i < conTib2_r.size(); i++) {
            Vec3 edge1_t = tibPoints_transf[facesTib2[i][1] - 1] - tibPoints_transf[facesTib2[i][0] - 1];
            Vec3 edge2_t = tibPoints_transf[facesTib2[i][2] - 1] - tibPoints_transf[facesTib2[i][1] - 1];
            Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
            At2[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
        }*/
        multipliers1 = GenerateMultList_spheres(pairs1_list, cont_centers_Fem, cont_centers_tib_transf1, conFem_r, conTib1_r);
        multipliers2 = GenerateMultList_spheres(pairs2_list, cont_centers_Fem, cont_centers_tib_transf2, conFem_r, conTib2_r);
    }

    std::cout << "multipliers1" << multipliers1 << std::endl;
    std::cout << "multipliers2" << multipliers2 << std::endl;


    Real poisson =0.46;
    Real E = 400 * 1e6;
    Real h = 0.006;
    Vec3 SumForces1;
    Vec3 SumMoments1;
    SumForces1.setToZero();
    SumMoments1.setToZero();
    Vec3 SumForces2;
    Vec3 SumMoments2;
    SumForces2.setToZero();
    SumMoments2.setToZero();

    CalculateForceCompartment(femPoints, facesFem, tibPoints_transf, facesTib1, pairs1_list, SumForces1, SumMoments1, poisson, E, h, originTib_G, multipliers1, knee_trans, knee_rot, At1, tibplanes1, cont_centers_Fem, cont_centers_tib_transf1,pvec1);
    std::cout << "forces=" << SumForces1 << " moments=" << SumMoments1 << std::endl;
    CalculateForceCompartment(femPoints, facesFem, tibPoints_transf, facesTib2, pairs2_list, SumForces2, SumMoments2, poisson, E, h, originTib_G, multipliers2, knee_trans, knee_rot, At2, tibplanes2, cont_centers_Fem, cont_centers_tib_transf2,pvec2);
    SumForces = SumForces1 + SumForces2;
    SumMoments = SumMoments1 + SumMoments2;
    SumForces_vert_Med = SumForces1[1]; // tibial part 1 is medial
    SumForces_vert_Lat = SumForces2[1]; // tibial part 2 is lateral
    std::cout << "sumforces=" << SumForces << std::endl;

}
// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    //// OpenSim model: create components
    ///// Model
    //OpenSim::Model* model;
    ///// Bodies
    //OpenSim::Body* femoral_component;
    //OpenSim::Body* tibial_tray;
    ///// Joints
    //OpenSim::CustomJoint* ground_femoralcomponent;
    //OpenSim::CustomJoint* knee;
    //
    //// OpenSim model: initialize components
    ///// Model
    //model = new OpenSim::Model();
    ///// Body specifications
    //femoral_component = new OpenSim::Body("femoral_component", 0.00914431311487012, Vec3(0, 0, 0), Inertia(9.14431311487012e-005, 9.14431311487012e-005, 9.14431311487012e-005, 0, 0, 0));
    //tibial_tray = new OpenSim::Body("tibial_tray", 0.00914431311487012, Vec3(0, 0, 0), Inertia(9.14431311487012e-005, 9.14431311487012e-005, 9.14431311487012e-005, 0, 0, 0));

    ///// Joints
    ///// Knee_r transform
    //SpatialTransform st_knee_r; // should it be flexion, adduction, rotation?
    //st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    //st_knee_r[0].setFunction(new LinearFunction());
    //st_knee_r[0].setAxis(Vec3(0, 0, 1));
    //st_knee_r[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_r", 1, 1));
    //st_knee_r[1].setFunction(new LinearFunction());
    //st_knee_r[1].setAxis(Vec3(1, 0, 0));
    //st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_rotation_r", 1, 1));
    //st_knee_r[2].setFunction(new LinearFunction());
    //st_knee_r[2].setAxis(Vec3(0, -1, 0));
    //st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_tx_r", 1, 1));
    //st_knee_r[3].setFunction(new LinearFunction());
    //st_knee_r[3].setAxis(Vec3(1, 0, 0));
    //st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_ty_r", 1, 1));
    //st_knee_r[4].setFunction(new LinearFunction());
    //st_knee_r[4].setAxis(Vec3(0, 1, 0));
    //st_knee_r[5].setCoordinateNames(OpenSim::Array<std::string>("knee_tz_r", 1, 1));
    //st_knee_r[5].setFunction(new LinearFunction());
    //st_knee_r[5].setAxis(Vec3(0, 0, 1));

    ///// Joint specifications
    //ground_femoralcomponent = new WeldJoint("ground_femoralcomponent", model->getGround(), Vec3(0), Vec3(0), *femoral_component, Vec3(0), Vec3(0));
    //knee = new CustomJoint("knee_r", *tibial_tray, Vec3(0, 0, 0), Vec3(0), *femoral_component, Vec3(0), Vec3(0), st_knee_r);
   
    ///// Add bodies and joints to model
    //model->addBody(femoral_component);	model->addJoint(ground_femoral_component);
    //model->addBody(tibial_tray);	    model->addJoint(knee);

    //// Initialize system and state
    //SimTK::State* state;
    //state = new State(model->initSystem());

    ////CALL STATE INDEX MAPPING FUNCTION TO ACCOUNT FOR OPENSIM VS SIMBODY STATE ORDERS
    //Array<std::string> stateVars = model->getStateVariableNames();					//Assign string array with the state variable names
    //std::unordered_map<std::string, int> mapping = createSystemYIndexMap(*model);    //Call function
    //for (int i = 0; i < mapping.size(); ++i) std::cout << mapping[stateVars[i]] << " " << stateVars[i] << " " << i << " OpenSim" << std::endl; //Loop through each state name and print to the cmd window the corresponding Simbody index and the name

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);

    //// States and controls
    //Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states
    //
    //// Assign inputs to model variables
    ///// States
    //for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    //
    ///// Controls
    //for (int i = 0; i < 12; ++i) ua[i] = u[i];
    ///// OpenSim and Simbody have different state orders so we adjust manually
    ////...

    //// Set state variables and realize
    //model->setStateVariableValues(*state, QsUs);
    //model->realizeVelocity(*state);

    //Vec3 knee_trans = Vec3(model->getStateVariableValue(*state, "knee_r/knee_tx_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_ty_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_tz_r/value"));
    //Vec3 knee_rot_0 = Vec3(model->getStateVariableValue(*state, "knee_r/knee_angle_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_adduction_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_rotation_r/value"));
    //
    Vec3 knee_rot_0(0);
    Vec3 knee_rot(0);
    Vec3 knee_trans(0);
    for (int i = 0; i < 3; i++) {
        knee_rot_0[i] = x[i];
        knee_trans[i] = x[i + 3];
    }
    
    knee_rot = Vec3(-knee_rot_0[0], -knee_rot_0[1], knee_rot_0[2]); // change sign to be consistent with joint definition (reverse in .osim)

    // Compute Resulting contact wrench at the center of the tibia
    Vec3 KneeCont_SumForces;
    Vec3 KneeCont_SumMoments;
    KneeCont_SumForces.setToZero();
    KneeCont_SumMoments.setToZero();
    Real SumForces_vert_Lat;
    Real SumForces_vert_Med;
    
    int nfacestib1 = 0;
    int nfacestib2 = 0;
    if (nfacesTib == 49) {
        nfacestib1 = 26;
        nfacestib2 = 23;
    }
    else if (nfacesTib == 100) {
        nfacestib1 = 51;
        nfacestib2 = 49;
    }
    Vector pvec1(nfacestib1);
    pvec1.setToZero();
    Vector pvec2(nfacestib2);
    pvec2.setToZero();

    ComputeKneeContactForces(knee_trans, knee_rot, KneeCont_SumForces, KneeCont_SumMoments, SumForces_vert_Lat, SumForces_vert_Med, pvec1,pvec2);
    Vec3 KneeCont_SumForces_onTibialTray_inTibialTrayFrame = -KneeCont_SumForces;
    Vec3 KneeCont_SumMoments_onTibialTray_inTibialTrayFrame = -KneeCont_SumMoments;

    /*Vec3 KneeCont_SumForces_onTibialTray_inG=tibial_tray->
    VectorInGround(*state, KneeCont_SumForces_onTibialTray_inTibialTrayFrame);
    Vec3 KneeCont_SumMoments_onTibialTray_inG = tibial_tray->expressVectorInGround(*state, KneeCont_SumMoments_onTibialTray_inTibialTrayFrame);

    Vec3 KneeCont_SumForces_onFemoralComp_inG = -KneeCont_SumForces_onTibialTray_inG;
    Vec3 tibial_tray_Center_inG = tibial_tray->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 femoral_comp_Center_inG = femoral_component->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 KneeCont_SumMoments_onFemoralComp_inG = -KneeCont_SumMoments_onTibialTray_inG + cross(tibial_tray_Center_inG - femoral_comp_Center_inG, KneeCont_SumForces_onFemoralComp_inG);*/

  
    // Extract results
    /// Knee contact forces
    res[0][0] = value<T>(SumForces_vert_Med);
    res[0][1] = value<T>(SumForces_vert_Lat);
    std::cout << "SumForces_vert_Lat" << SumForces_vert_Lat << std::endl;
    std::cout << "SumForces_vert_Med" << SumForces_vert_Med << std::endl;

    for (int i = 0; i < nfacestib1; ++i) {
        res[0][i + 2] = value<T>(pvec1[i]);
    }
    std::cout << "nfacestib2=" << nfacestib2 << std::endl;
    for (int i = 0; i < nfacestib2; ++i) {
        res[0][i + 2 + nfacestib1] = value<T>(pvec2[i]);
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
    //Recorder u[NU];
    //Recorder p[NP];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    //for (int i = 0; i < NU; ++i) u[i] <<= 0;
    //for (int i = 0; i < NP; ++i) p[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}

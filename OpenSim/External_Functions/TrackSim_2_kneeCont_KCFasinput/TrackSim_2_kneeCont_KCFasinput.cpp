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

#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
    #include <limits.h>
#endif


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
constexpr int n_in = 4;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 34;        // # degrees of freedom (excluding locked)
constexpr int ndofr = 36;       // # degrees of freedom (including locked)
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NP = 54;          // # parameters
constexpr int NKCFM = 6;         // # KCF
constexpr int NR = ndof+6+6;    // # residual torques + # GRFs + # GRMs


//constexpr int nfacesTib = 49; //before 49
//constexpr int nfacesFem = 258;
//constexpr const char radForPairs[] = "05"; // 1 is 1 cm, 05 is 0.5 cm
//constexpr char* multiplier_method = "cylinders"; // multiplier method: "cylinders" or "spheres"

std::string getHostname() {
    char hostname[256];
#ifdef _WIN32
    DWORD size = sizeof(hostname);
    GetComputerNameA(hostname, &size);
#else
    gethostname(hostname, sizeof(hostname));
#endif
    return std::string(hostname);
}

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

//std::vector<Vec4> ReadDataDoublex4columns(std::string filename) {
//
//
//    std::ifstream  file(filename);
//
//    //
//    double v1, v2, v3, v4;
//    std::vector<Vec4> out_csv;
//
//    if (file.is_open()) {
//
//        std::string line;
//
//        while (std::getline(file, line)) {
//            std::istringstream iss(line);
//            std::string token;
//            std::vector<double> values;
//
//            while (std::getline(iss, token, ',')) {
//                values.push_back(std::stod(token)); // Convert string to double and add to values vector
//            }
//
//            if (values.size() >= 4) {
//                out_csv.emplace_back(values[0], values[1], values[2], values[3]); // Create Vec4 object and add to out_csv
//            }
//        }
//
//        file.close();
//    }
//    else {
//        std::cerr << "Error opening file." << std::endl;
//    }
//    return out_csv;
//}
//
//std::vector<std::vector<int>> ReadDataIntx3columns(std::string filename) {
//
//
//    std::ifstream  file(filename);
//
//    int v1, v2, v3;
//    std::vector<std::vector<int>> out_csv;
//
//    if (file.is_open()) {
//
//        std::string line;
//
//        while (std::getline(file, line)) {
//            std::istringstream iss(line);
//            std::string token;
//            std::vector<int> row;
//
//            while (std::getline(iss, token, ',')) {
//                row.push_back(std::stoi(token)); // Convert string to integer and add to the row vector
//            }
//
//            out_csv.push_back(row); // Add the row to the out_csv vector
//        }
//
//        file.close();
//
//    }
//    else {
//        std::cerr << "Error opening file." << std::endl;
//    }
//    return out_csv;
//}
//
//std::vector<std::vector<int>> ReadDataIntx2columns(std::string filename) {
//
//
//    std::ifstream  file(filename);
//
//    int v1, v2;
//    std::vector<std::vector<int>> out_csv;
//    if (file.is_open()) {
//        std::string line;
//
//        while (std::getline(file, line)) {
//            std::istringstream iss(line);
//            std::string token;
//            std::vector<int> row;
//
//            while (std::getline(iss, token, ',')) {
//                row.push_back(std::stoi(token)); // Convert string to integer and add to the row vector
//            }
//
//            out_csv.push_back(row); // Add the row to the out_csv vector
//        }
//        file.close();
//        }else {
//        std::cerr << "Error opening file." << std::endl;
//    }
//    return out_csv;
//}
//
//std::vector<Vec3> ReadDataDoublex3columns(std::string filename) {
//
//    
//    std::ifstream  file(filename);
//
//    //
//    double v1, v2, v3;
//    std::vector<Vec3> out_csv;
//    
//    int nrows = 0;
//    if (file.is_open()) {
//
//        std::string line;
//
//        while (std::getline(file, line)) {
//            std::istringstream iss(line);
//
//            std::string token;
//            std::vector<std::string> tokens;
//
//            // Split the line by comma and store tokens in a vector
//            while (std::getline(iss, token, ',')) {
//                tokens.push_back(token);
//            }
//
//            // Check if the line has at least three tokens (assuming three values in each row)
//            if (tokens.size() >= 3) {
//                double v1 = std::stod(tokens[0]);
//                double v2 = std::stod(tokens[1]);
//                double v3 = std::stod(tokens[2]);
//
//                out_csv.emplace_back(v1, v2, v3); // Add Vec3 object to the vector
//            }
//
//            nrows = nrows + 1;
//
//        }
//
//        file.close();
//    }else {
//        std::cerr << "Error opening file." << std::endl;
//    }
//  
//    return out_csv;
//}
//
//Mat33 R_aux(Vec3 knee_trans, Vec3 knee_rot) {
//    Recorder psi = knee_rot[0];
//    Recorder theta = knee_rot[1];
//    Recorder phi = knee_rot[2];
//
//    Mat33 R1(0);
//    R1.set(0, 0, cos(psi));
//    R1.set(0, 1, -sin(psi));
//    R1.set(0, 2, 0.0);
//    R1.set(1, 0, sin(psi));
//    R1.set(1, 1, cos(psi));
//    R1.set(1, 2, 0.0);
//    R1.set(2, 0, 0.0);
//    R1.set(2, 1, 0.0);
//    R1.set(2, 2, 1.0);
//    Mat33 R2(0);
//    R2.set(0, 0, 1.0);
//    R2.set(0, 1, 0.0);
//    R2.set(0, 2, 0.0);
//    R2.set(1, 0, 0.0);
//    R2.set(1, 1, cos(theta));
//    R2.set(1, 2, -sin(theta));
//    R2.set(2, 0, 0.0);
//    R2.set(2, 1, sin(theta));
//    R2.set(2, 2, cos(theta));
//    Mat33 R3(0);
//    R3.set(0, 0, cos(phi));
//    R3.set(0, 1, 0.0);
//    R3.set(0, 2, sin(phi));
//    R3.set(1, 0, 0.0);
//    R3.set(1, 1, 1.0);
//    R3.set(1, 2, 0.0);
//    R3.set(2, 0, -sin(phi));
//    R3.set(2, 1, 0.0);
//    R3.set(2, 2, cos(phi));
//    Mat33 R = R1*R2*R3;
//
//    return R;
//}
//
//Mat44 ftransf_function(Vec3 knee_trans, Vec3 knee_rot) {
//
//    Mat33 R = R_aux(knee_trans, knee_rot);
//    
//    Vec3 aux = R*(-knee_trans);
//    
//    Mat44 Rtrans0042(1);
//    Rtrans0042.set(1, 3, 0.042);
//        
//    Mat44 Rtranstib_4x4(1);
//    Rtranstib_4x4.set(0, 3, aux[0]);
//    Rtranstib_4x4.set(1, 3, aux[1]);
//    Rtranstib_4x4.set(2, 3, aux[2]);
//
//    std::cout << Rtranstib_4x4 << std::endl;
//
//    Mat44 R_rot(1);
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++) {
//            R_rot.set(i, j, R[i][j]);
//        }
//    }
//
//    Mat44 R_tib = Rtranstib_4x4*Rtrans0042*R_rot;
//    return R_tib;
//}
//
//Real CheckContact(Real overlap) {
//    Real k = 1e3; // k subject to change
//    Real multiplier = (tanh(k * overlap) + 1.0) / 2.0;
//    /*std::cout << multiplier << std::endl;*/
//    return multiplier;
//}
//
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
//
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
//
//void CalculateIntersection_cylinders(Vec3  cont_centers_tib, Vec3 cont_centers_fem, Vec3& d) {
//    d = cont_centers_fem - cont_centers_tib;
//}
//
//void CalculateIntersection_spheres(Vector_<Vec3> fem_Points, Vector_<Vec3> tib_Points, Vec3 &d) {
//    Vec3 Ctib(0.0);
//    Vec3 Cfem(0.0);
//    // mean point of the tib face
//    for (int i = 0; i < 3; i++) {
//        Ctib = mean (Vec3(tib_Points[0][i], tib_Points[1][i], tib_Points[2][i]));
//    }
//    
//    // mean point of the fem face
//    for (int i = 0; i < 3; i++) {
//        Cfem[i] = mean(Vec3(fem_Points[0][i], fem_Points[1][i], fem_Points[2][i]));
//    }
//
//    d = Cfem - Ctib;
//
//}
//
//void CalculateMaximumPenetration(Vector_<Vec3> d_v, Vector_<Vec3> nt_v, Real &maxpen, Vec3 &nt_l, Vector_<Real> mults) {
//    Vector_<Real> proj(d_v.size());
//    //Vector_<Vec3> dist_v(d_v.size());
//    Vector_<Real> pen(d_v.size());
//    /*std::cout << proj.size() << std::endl;
//    std::cout << mults.size() << std::endl;
//    std::cout << mults << std::endl;*/
//    for (int i = 0; i < d_v.size(); i++) {
//        proj[i] = dot(d_v[i], nt_v[i]); // projection of distance between tibial and femoral faces to the normal of tibial face
//        pen[i] = -proj[i]*mults[i]; // pen is for penetration
//        //pen[i] = -proj[i]; // pen is for penetration
//    }
//    Real k = 1e4;
//    std::cout << "proj" << proj << std::endl;
//    std::cout << "pen" << pen << std::endl;
//    std::cout << "mults=" << mults << std::endl;
//    std::cout << "sum mults=" << sum(mults) << std::endl;
//    std::cout << "k*pen" << k*pen << std::endl;
//    std::cout << "exp(k*pen))" << exp(k * pen) << std::endl;
//
//    //maxpen = (log(sum(exp(k*pen))+1e-16) / k); //version logSum
//    maxpen = (log((1.0/pen.size())*sum(exp(k * pen))+1e-16) / k); //version mellowmax
//    //maxpen = max(pen); // version nosmooth
//    
//    //Real auxval(0.0);
//    //for (int i = 0; i < pen.size(); i++) {
//    //    auxval = auxval + pow(pen[i], 10.0);
//    //}
//    //Real auxval2 = pow(auxval, 1.0 / 10);
//    //maxpen = sum(pen*(0.5 * tanh(k * (pen - auxval2 + 1e-5)) + 0.5)); // version p norm
//
//    std::cout << "sum(exp(k * pen))" << sum(exp(k * pen)) << std::endl;
//    std::cout << "maxpen=" << maxpen << std::endl;
//    std::cout << "1.0 / pen.size()" << 1.0/pen.size() << std::endl;
//    std::cout << "true max pen=" << max(pen) << std::endl;
//    std::cout << "nt_v=" << nt_v << std::endl;
//    //maxpen = max(pen);
//    nt_l = nt_v[0];
//}
//
//Real CalculatePressure(Real poisson, Real E, Real d, Real h) {
//    Real k = 5e5;
//    Real pen = d;
//
//    Real p_init = ((1 - poisson)*E / ((1 + poisson)*(1 - 2 * poisson)))*pen / h;
//
//    Real p = p_init*(1 + tanh(k*pen)) / 2;
//
//    std::cout << "pen" << pen << std::endl;
//    std::cout << "p" << p << std::endl;
//
//    return p;
//}
//
//void CalculateForceCompartment(std::vector<Vec3> femPoints, std::vector<std::vector<int>> facesFem, Vector_<Vec3> tibPoints_transf, std::vector<std::vector<int>> facesTib, std::vector<std::vector<int>> pairs_list, Vec3& SumForces, Vec3& SumMoments, Real poisson, Real E, Real h, Vec3 originTib_G, Vector_<Real> multipliers, Vec3 knee_trans, Vec3 knee_rot, Vector_<Real> At, Vector_<Vec4> tibplanes, Vector_<Vec3> cont_centers_Fem, Vector_<Vec3> cont_centers_tib_transf) {
//    Vector_<Vec3> d(pairs_list.size());
//    Vector_<Vec3> nt(pairs_list.size());
//    Vector_<Vec3> force_s_l(0);
//    Vector_<Real> face_s(0);
//    Vector_<Vec3> mom_O(0);
//
//    int k = 1; //number of contacting elements in femur
//    int l = 1; //count number of elements contacting element k of femur
//
//    //std::cout << "facesFem.size()" << facesFem.size() << std::endl;
//    //for (int i = 0; i < facesFem.size(); i++) {
//    //    std::cout << "faces_fem[" << i << "]=" << facesFem.at(i).at(0) << " " << facesFem.at(i).at(1) << facesFem.at(i).at(2) << std::endl;
//    //}
//    std::cout << "femPoints" << std::endl;
//    for (int i = 0; i < tibplanes.size(); i++) {
//        std::cout << "tibplanes[pairs_list[i][0]]=" << tibplanes[pairs_list[i][0]] << std::endl;
//    }
//
//    for (int i = 0; i < size(femPoints); i++) {
//        std::cout << femPoints[i] << std::endl;
//    }
//
//    for (int i = 0; i < pairs_list.size(); i++) {
//        Vec3 d_aux;
//        ///////// check Ctib_aux
//
//        if (strcmp(multiplier_method, "cylinders") == 0) {
//            std::cout << multiplier_method << std::endl;
//            CalculateIntersection_cylinders(cont_centers_tib_transf[pairs_list[i][0] - 1], cont_centers_Fem[pairs_list[i][1] - 1], d_aux);
//        }
//        else if (strcmp(multiplier_method, "spheres") == 0) {
//            std::vector<int> facesFem_ind = facesFem[pairs_list[i][1] - 1];
//            Vector_<Vec3> fem_Points_ind(3);
//            for (int j = 0; j < 3; j++) {
//                fem_Points_ind[j] = femPoints[facesFem_ind[j] - 1];
//            }
//
//            std::vector<int> facesTib_ind = facesTib[pairs_list[i][0] - 1];
//            Vector_<Vec3> tib_Points_ind(3);
//            for (int j = 0; j < 3; j++) {
//                tib_Points_ind[j] = tibPoints_transf[facesTib_ind[j] - 1];
//            }
//            CalculateIntersection_spheres(fem_Points_ind, tib_Points_ind, d_aux); // first element in pairs is tibia, second femur
//        }
//        d[i] = d_aux;
//        nt[i] = Vec3(tibplanes[pairs_list[i][0]-1][0], tibplanes[pairs_list[i][0]-1][1], tibplanes[pairs_list[i][0]-1][2]);
//        std::cout << multipliers << std::endl;
//
//
//        if (i > 0) {
//            if ((pairs_list[i][0] == pairs_list[i - 1][0]) && (i < pairs_list.size() - 1)) {
//                l = l + 1;
//            }
//            else {
//                Vector_<Vec3> d_aux_list(l - 1);
//                Vector_<Vec3> nt_aux_list(l - 1);
//                Vector_<Real> multipliers_list(l - 1);
//                Real maxpen;
//                Vec3 nt_l;
//                for (int j = 0; j < l - 1; j++) {
//                    d_aux_list[j] = d[i - l + j + 1];
//                    nt_aux_list[j] = nt[i - l + j + 1];
//                    multipliers_list[j] = multipliers[i - l + j + 1];
//                }
//                CalculateMaximumPenetration(d_aux_list, nt_aux_list, maxpen, nt_l, multipliers_list);
//                std::cout << maxpen << std::endl;
//                std::cout << At[pairs_list[i][0] - 1] << std::endl;
//                std::cout << "nt_l=" << nt_l << std::endl;
//                Real p = CalculatePressure(poisson, E, maxpen, h);
//
//                force_s_l.resizeKeep(k);
//                force_s_l[k - 1] = p * At[pairs_list[i][0] - 1] * nt_l;
//                face_s.resizeKeep(k);
//                face_s[k - 1] = pairs_list[i - 2][0];
//                mom_O.resizeKeep(k);
//                mom_O[k - 1] = cross(cont_centers_tib_transf[pairs_list[i][0] - 1] - originTib_G, force_s_l[k - 1]);
//
//                k = k + 1;
//                l = 2;
//                while (k < pairs_list[i][0]) {
//                    k = k + 1;
//                }
//                ///////////////
//
//            }
//        }
//        else {
//            l = l + 1;
//        }
//
//    }
//    for (int i = 0; i < nt.size(); i++){
//        std::cout << "nt=" << nt[i] << std::endl;
//    }
//    Vec3 Sum_Force_G;
//    Vec3 Sum_Moments_G;
//    Sum_Force_G.setToZero();
//    Sum_Moments_G.setToZero();
//    for (int i = 0; i < force_s_l.size(); i++) {
//        Sum_Force_G = Sum_Force_G + force_s_l[i];
//        Sum_Moments_G = Sum_Moments_G + mom_O[i];
//    }
//    
//    Mat33 Mrottib = R_aux(knee_trans, knee_rot);
//    
//    Mat33 MrottibTrans = Mrottib.transpose();
//    SumForces = MrottibTrans*Sum_Force_G;
//    SumMoments = MrottibTrans*Sum_Moments_G;
//    
//    std::cout << SumForces << std::endl;
//
//}
//
//void ComputeKneeContactForces(Vec3 knee_trans, Vec3 knee_rot, Vec3 &SumForces, Vec3 &SumMoments, Real &SumForces_vert_Lat, Real &SumForces_vert_Med) {
//    //// Read Geometry Information
//    std::string hostname = getHostname();
//    std::string root_folder;
//    if (hostname == "DESKTOP-U8CF7T5") {
//        root_folder = "C:/Gil/MeshesInAD/contactsKneeProsthesis/";
//    }
//    else
//    {
//        root_folder = "";
//    };
//
//    // Read Points
//    std::string filename_femPoints = root_folder + "femPoints_" + std::to_string(nfacesFem) + ".csv";
//    std::vector<Vec3> femPoints = ReadDataDoublex3columns(filename_femPoints);
//    std::string filename_tibPoints=root_folder+"tibPoints_" + std::to_string(nfacesTib) + ".csv";
//    std::vector<Vec3> tibPoints = ReadDataDoublex3columns(filename_tibPoints);
//    // Read Faces
//    std::string filename_facesFem = root_folder + "facesFem_" + std::to_string(nfacesFem) + ".csv";
//    std::vector<std::vector<int>> facesFem = ReadDataIntx3columns(filename_facesFem);
//    std::string filename_facesTib1 = root_folder + "facesTib1_" +std::to_string(nfacesTib) + ".csv";
//    std::vector<std::vector<int>> facesTib1 = ReadDataIntx3columns(filename_facesTib1);
//    std::string filename_facesTib2 = root_folder + "facesTib2_" + std::to_string(nfacesTib) + ".csv";
//    std::vector<std::vector<int>> facesTib2 = ReadDataIntx3columns(filename_facesTib2);
//    // Read centers
//    std::string filename_conFem = root_folder + "ConFem_" + std::to_string(nfacesFem) + ".csv";
//    std::vector<Vec4> conFem = ReadDataDoublex4columns(filename_conFem);
//    std::string filename_conTibia1 = root_folder + "ConTib1_" +std::to_string(nfacesTib) + ".csv";
//    std::vector<Vec4> conTib1 = ReadDataDoublex4columns(filename_conTibia1);
//    std::string filename_conTibia2 = root_folder + "ConTib2_" + std::to_string(nfacesTib) + ".csv";
//    std::vector<Vec4> conTib2 = ReadDataDoublex4columns(filename_conTibia2);
//
//
//    //Read pairs
//    std::string filename_pairs1 = root_folder + "pairs1_" + std::to_string(nfacesTib) + "x" + std::to_string(nfacesFem) + "_at" + radForPairs + "cm.csv";
//    std::string filename_pairs2 = root_folder + "pairs2_" + std::to_string(nfacesTib) + "x" + std::to_string(nfacesFem) + "_at" + radForPairs + "cm.csv";
//
//    std::cout << filename_pairs1 << std::endl;
//
//    std::vector<std::vector<int>> pairs1_list = ReadDataIntx2columns(filename_pairs1);
//    std::vector<std::vector<int>> pairs2_list = ReadDataIntx2columns(filename_pairs2);
//    
//    std::cout << filename_pairs1 << std::endl;
//    std::cout << pairs1_list[4][0] << " " << pairs1_list[4][1] << std::endl;  
//
//    // Sum shift translation to the femur
//    for (int i = 0; i < size(femPoints); i++) {
//        femPoints[i] = femPoints[i] + Vec3(0.0, 0.042, 0.0);
//    }
//
//    Vector_<Vec3> cont_centers_Fem(conFem.size(), Vec3(0));
//    Vector_<Real> conFem_r(conFem.size());
//    for (int i = 0; i < conFem.size(); i++) {
//        cont_centers_Fem[i]=Vec3(conFem[i][0],conFem[i][1], conFem[i][2]) + Vec3(0.0, 0.042, 0.0);
//        conFem_r[i] = conFem[i][3];
//    }
//
//    std::cout << facesFem[0][0] << facesFem[0][1] << facesFem[0][2] << std::endl;
//    std::cout << facesFem[1][0] << facesFem[1][1] << facesFem[1][2] << std::endl;
//
// 
//    //test with numerical values
//    //knee_trans = Vec3(0.01, 0.03, 0.008);
//    //knee_rot = Vec3(0.1, 0.03, 0.025);
//    std::cout << "knee_trans=" << knee_trans << std::endl;
//
//    //// Apply transformations
//    // Get matrix transformation of tibia with respect to the femur
//    Mat44 Mtransf_tib = ftransf_function(knee_trans, knee_rot);
//    std::cout << "Mtransf_tib=" << Mtransf_tib << std::endl;
//
//    // Apply the transformation to all points of the tibia
//    Vec4 tibPoints_transf_aux(1);
//    Vector_<Vec3> tibPoints_transf(size(tibPoints));
//    Vec4 cont_centers_tib_aux1(1);
//    Vector_<Vec3> cont_centers_tib_transf1(conTib1.size());
//    Vector_<Real> conTib1_r(conTib1.size());
//    Vec4 cont_centers_tib_aux2(1);
//    Vector_<Vec3> cont_centers_tib_transf2(conTib2.size());
//    Vector_<Real> conTib2_r(conTib2.size());
//    for (int i = 0; i < size(tibPoints); i++) {
//        tibPoints_transf_aux = Mtransf_tib*Vec4(tibPoints[i][0], tibPoints[i][1], tibPoints[i][2], 1.0);
//
//        for (int j = 0; j < 3; j++) {
//            tibPoints_transf[i][j] = tibPoints_transf_aux[j];
//        }
//    }
//    std::cout << conTib1.size() << std::endl;
//    for (int i = 0; i < conTib1.size(); i++) { // is the size of conTib correct?
//        cont_centers_tib_aux1 = Mtransf_tib * Vec4(conTib1[i][0], conTib1[i][1], conTib1[i][2], 1.0);
//        conTib1_r[i] = conTib1[i][3];
//
//        for (int j = 0; j < 3; j++) {
//            cont_centers_tib_transf1[i][j] = cont_centers_tib_aux1[j];
//        }
//    }
//    std::cout << conTib2.size() << std::endl;
//    for (int i = 0; i < conTib2.size(); i++) { // is the size of conTib correct?
//        cont_centers_tib_aux2 = Mtransf_tib * Vec4(conTib2[i][0], conTib2[i][1], conTib2[i][2], 1.0);
//        conTib2_r[i] = conTib2[i][3];
//
//        for (int j = 0; j < 3; j++) {
//            cont_centers_tib_transf2[i][j] = cont_centers_tib_aux2[j];
//        }
//    }
//
//    // Calculate origin of the tibia...
//    Vec4 originTib_G4 = Mtransf_tib*Vec4(0, 0, 0, 1);
//    Vec3 originTib_G = Vec3(originTib_G4[0], originTib_G4[1], originTib_G4[2]);
//    // 
//
//    // Calculate multipliers for all pairs at this instant, and areas of tibia faces
//    Vector_<Real> multipliers1(pairs1_list.size(), 0.0);
//    Vector_<Real> multipliers2(pairs1_list.size(), 0.0);
//    Vector_<Real> At1(conTib1_r.size(), 0.0);
//    Vector_<Real> At2(conTib1_r.size(), 0.0);
//    Vector_<Vec4> tibplanes1(conTib1_r.size(), Vec4(0.0));
//    Vector_<Vec4> tibplanes2(conTib2_r.size(), Vec4(0.0));
//    for (int i = 0; i < conTib1_r.size(); i++) {
//        Vec3 edge1_t = tibPoints_transf[facesTib1[i][1] - 1] - tibPoints_transf[facesTib1[i][0] - 1];
//        Vec3 edge2_t = tibPoints_transf[facesTib1[i][2] - 1] - tibPoints_transf[facesTib1[i][1] - 1];
//        //// unit normal vector to triangle of tibia1
//        Vec3 nt_i = cross(edge1_t, edge2_t);
//        nt_i = nt_i.normalize();
//        Real Dplanetib = -dot(nt_i, cont_centers_tib_transf1[i]);
//        tibplanes1[i][0] = nt_i[0];
//        tibplanes1[i][1] = nt_i[1];
//        tibplanes1[i][2] = nt_i[2];
//        tibplanes1[i][3] = Dplanetib;
//        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
//        At1[i] = (1.0 / 2.0) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
//    }
//    for (int i = 0; i < conTib2_r.size(); i++) {
//        Vec3 edge1_t = tibPoints_transf[facesTib2[i][1] - 1] - tibPoints_transf[facesTib2[i][0] - 1];
//        Vec3 edge2_t = tibPoints_transf[facesTib2[i][2] - 1] - tibPoints_transf[facesTib2[i][1] - 1];
//        //// unit normal vector to triangle of tibia1
//        Vec3 nt_i = cross(edge1_t, edge2_t);
//        nt_i = nt_i.normalize();
//        Real Dplanetib = -dot(nt_i, cont_centers_tib_transf1[i]);
//        tibplanes2[i][0] = nt_i[0];
//        tibplanes2[i][1] = nt_i[1];
//        tibplanes2[i][2] = nt_i[2];
//        tibplanes2[i][3] = Dplanetib;
//        Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
//        At2[i] = (1.0 / 2.0) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
//    }
//    std::cout << "At1=" << At1 << std::endl;
//    std::cout << "At2=" << At2 << std::endl;
//    if (strcmp(multiplier_method, "cylinders") == 0) {
//        // method of cylinders to generate multipliers
//
//        /// first compute planes
//        Vector_<Vec4> femplanes(conFem_r.size(), Vec4(0.0));
//        for (int i = 0; i < conFem_r.size(); i++) {
//            Vec3 edge1_s = femPoints[facesFem[i][1]-1] - femPoints[facesFem[i][0]-1];
//            Vec3 edge2_s = femPoints[facesFem[i][2]-1] - femPoints[facesFem[i][1]-1];
//            //// unit normal vector to triangle of femur
//            Vec3 ns_i = cross(edge1_s, edge2_s);
//            ns_i = ns_i.normalize();
//            Real Dplanefem = -dot(ns_i, cont_centers_Fem[i]);
//            femplanes[i][0] = ns_i[0];
//            femplanes[i][1] = ns_i[1];
//            femplanes[i][2] = ns_i[2];
//            femplanes[i][3] = Dplanefem;
//        }
//        
//        // Calculate multipliers for all pairs at this instant
//        multipliers1 = GenerateMultList_cylinders(pairs1_list, cont_centers_Fem, cont_centers_tib_transf1, tibplanes1, femplanes, conFem_r);
//        multipliers2 = GenerateMultList_cylinders(pairs2_list, cont_centers_Fem, cont_centers_tib_transf2, tibplanes2, femplanes, conFem_r);
//       
//    }
//
//    else if (strcmp(multiplier_method, "spheres") == 0) {
//       /* for (int i = 0; i < conTib1_r.size(); i++) {
//            Vec3 edge1_t = tibPoints_transf[facesTib1[i][1] - 1] - tibPoints_transf[facesTib1[i][0] - 1];
//            Vec3 edge2_t = tibPoints_transf[facesTib1[i][2] - 1] - tibPoints_transf[facesTib1[i][1] - 1];
//            Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
//            At1[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
//        }
//        for (int i = 0; i < conTib2_r.size(); i++) {
//            Vec3 edge1_t = tibPoints_transf[facesTib2[i][1] - 1] - tibPoints_transf[facesTib2[i][0] - 1];
//            Vec3 edge2_t = tibPoints_transf[facesTib2[i][2] - 1] - tibPoints_transf[facesTib2[i][1] - 1];
//            Real anglec = acos(dot(edge1_t, edge2_t) / (edge1_t.norm() * edge2_t.norm()));
//            At2[i] = (1 / 2) * edge1_t.norm() * edge2_t.norm() * sin(anglec);
//        }*/
//        multipliers1 = GenerateMultList_spheres(pairs1_list, cont_centers_Fem, cont_centers_tib_transf1, conFem_r, conTib1_r);
//        multipliers2 = GenerateMultList_spheres(pairs2_list, cont_centers_Fem, cont_centers_tib_transf2, conFem_r, conTib2_r);
//    }
//
//    std::cout << "multipliers1" << multipliers1 << std::endl;
//    std::cout << "multipliers2" << multipliers2 << std::endl;
//
//
//    Real poisson =0.46;
//    Real E = 400 * 1e6;
//    Real h = 0.006;
//    Vec3 SumForces1;
//    Vec3 SumMoments1;
//    SumForces1.setToZero();
//    SumMoments1.setToZero();
//    Vec3 SumForces2;
//    Vec3 SumMoments2;
//    SumForces2.setToZero();
//    SumMoments2.setToZero();
//
//    CalculateForceCompartment(femPoints, facesFem, tibPoints_transf, facesTib1, pairs1_list, SumForces1, SumMoments1, poisson, E, h, originTib_G, multipliers1, knee_trans, knee_rot, At1, tibplanes1, cont_centers_Fem, cont_centers_tib_transf1);
//    std::cout << "forces=" << SumForces1 << " moments=" << SumMoments1 << std::endl;
//    CalculateForceCompartment(femPoints, facesFem, tibPoints_transf, facesTib2, pairs2_list, SumForces2, SumMoments2, poisson, E, h, originTib_G, multipliers2, knee_trans, knee_rot, At2, tibplanes2, cont_centers_Fem, cont_centers_tib_transf2);
//    SumForces = SumForces1 + SumForces2;
//    SumMoments = SumMoments1 + SumMoments2;
//    SumForces_vert_Med = SumForces1[1]; // tibial part 1 is medial
//    SumForces_vert_Lat = SumForces2[1]; // tibial part 2 is lateral
//    std::cout << "sumforces=" << SumForces << std::endl;
//
//}
//// Function F
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
    knee_l = new CustomJoint("knee_l", *tibia_l, Vec3(0), Vec3(0), *femur_l, Vec3(-0.00451221232146798, -0.396907245921447, 0), Vec3(0), st_knee_l);
    femoral_component_weld = new WeldJoint("femoral_component_weld", *femur_r, Vec3(0, -0.39335, 0), Vec3(3.083, 0, -3.083), *femoral_component, Vec3(0), Vec3(0));
    knee_r = new CustomJoint("knee_r", *tibial_tray, Vec3(0, 0, 0), Vec3(0), *femoral_component, Vec3(0), Vec3(0), st_knee_r);
    tibial_tray_weld = new WeldJoint("tibial_tray_weld", *tibial_tray, Vec3(0, 0.044254, 0), Vec3(0, 3.1416, 0), *tibia_r, Vec3(0, 0, 0), Vec3(0, 0, 0));
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.44751, 0), Vec3(-0.041214, 0.0031538, -0.050218), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.44751, 0), Vec3(0.041214, 0.0031538, -0.050218), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
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

    //CALL STATE INDEX MAPPING FUNCTION TO ACCOUNT FOR OPENSIM VS SIMBODY STATE ORDERS
    Array<std::string> stateVars = model->getStateVariableNames();					//Assign string array with the state variable names
    std::unordered_map<std::string, int> mapping = createSystemYIndexMap(*model);    //Call function
    for (int i = 0; i < mapping.size(); ++i) std::cout << mapping[stateVars[i]] << " " << stateVars[i] << " " << i << " OpenSim" << std::endl; //Loop through each state name and print to the cmd window the corresponding Simbody index and the name

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);
    std::vector<T> p(arg[2], arg[2] + NP);
    std::vector<T> KCF_in(arg[3], arg[3] + NKCFM);

    // States and controls
    T ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    T up[NP]; /// contact model parameters - parameters
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states
    T KCF[NKCFM];

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// pro_sup dofs are locked so Qs and Qdots are hard coded
    QsUs[NX] = 1.51;
    QsUs[NX+1] = 0;
    QsUs[NX+2] = 1.51;
    QsUs[NX+3] = 0;
    /// Controls
    for (int i = 0; i < 12; ++i) ua[i] = u[i];
    /// OpenSim and Simbody have different state orders so we adjust manually
    ua[12] = u[23]; /// 12 Simbody is 23 OpenSim
    ua[13] = u[24]; /// 13 Simbody is 24 OpenSim
    ua[14] = u[25]; /// 14 Simbody is 25 OpenSim
    ua[15] = u[12]; /// 15 Simbody is 12 OpenSim
    ua[16] = u[26]; /// 16 Simbody is 26 OpenSim
    ua[17] = u[27]; /// 17 Simbody is 27 OpenSim
    ua[18] = u[28]; /// 18 Simbody is 28 OpenSim
    ua[19] = u[29]; /// 19 Simbody is 29 OpenSim
    ua[20] = u[30]; /// 20 Simbody is 30 OpenSim
    ua[21] = u[31]; /// 21 Simbody is 31 OpenSim
    ua[22] = u[13]; /// 22 Simbody is 13 OpenSim
    ua[23] = u[14]; /// 23 Simbody is 14 OpenSim
    ua[24] = u[15]; /// 24 Simbody is 15 OpenSim
    ua[25] = u[16]; /// 25 Simbody is 16 OpenSim
    ua[26] = u[17]; /// 26 Simbody is 17 OpenSim
    ua[27] = u[18]; /// 27 Simbody is 18 OpenSim
    ua[28] = u[19]; /// 28 Simbody is 19 OpenSim
    ua[29] = u[32]; /// 29 Simbody is 32 OpenSim
    ua[30] = u[33]; /// 30 Simbody is 33 OpenSim
    ua[31] = u[21]; /// 31 Simbody is 21 OpenSim

    /// pro_sup dofs are locked so Qs and Qdots are hard coded
    ua[32] = 0;
    ua[33] = 0;
    ua[34] = u[20]; /// 34 Simbody is 20 OpenSim
    ua[35] = u[22]; /// 34 Simbody is 22 OpenSim

    /// Parameters
    for (int i = 0; i < NP; ++i) up[i] = p[i];

    /// KCF
    for (int i = 0; i < NKCFM; ++i) KCF[i] = KCF_in[i];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    //Vec3 knee_trans = Vec3(model->getStateVariableValue(*state, "knee_r/knee_tx_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_ty_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_tz_r/value"));
    //Vec3 knee_rot_0 = Vec3(model->getStateVariableValue(*state, "knee_r/knee_angle_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_adduction_r/value"),
    //    model->getStateVariableValue(*state, "knee_r/knee_rotation_r/value"));
    //
    //Vec3 knee_rot(0);
    //knee_rot = Vec3(-knee_rot_0[0], -knee_rot_0[1], knee_rot_0[2]); // change sign to be consistent with joint definition (reverse in .osim)

    //std::cout << "pelvis_ty=" << model->getStateVariableValue(*state, "knee_r/pelvis_ty/value") << std::endl;

    //// Compute Resulting contact wrench at the center of the tibia
    //Vec3 KneeCont_SumForces;
    //Vec3 KneeCont_SumMoments;
    //KneeCont_SumForces.setToZero();
    //KneeCont_SumMoments.setToZero();
    //Real SumForces_vert_Lat;
    //Real SumForces_vert_Med;
    //

    //ComputeKneeContactForces(knee_trans, knee_rot, KneeCont_SumForces, KneeCont_SumMoments, SumForces_vert_Lat, SumForces_vert_Med);
    
    Vec3 KneeCont_SumForces;
    Vec3 KneeCont_SumMoments;
    for (int i = 0; i < 3; ++i) {
        KneeCont_SumForces[i] = KCF[i];
        KneeCont_SumMoments[i] = KCF[i + 3];
    }
    
    Vec3 KneeCont_SumForces_onTibialTray_inTibialTrayFrame = -KneeCont_SumForces;
    Vec3 KneeCont_SumMoments_onTibialTray_inTibialTrayFrame = -KneeCont_SumMoments;

    Vec3 KneeCont_SumForces_onTibialTray_inG=tibial_tray->expressVectorInGround(*state, KneeCont_SumForces_onTibialTray_inTibialTrayFrame);
    Vec3 KneeCont_SumMoments_onTibialTray_inG = tibial_tray->expressVectorInGround(*state, KneeCont_SumMoments_onTibialTray_inTibialTrayFrame);

    Vec3 KneeCont_SumForces_onFemoralComp_inG = -KneeCont_SumForces_onTibialTray_inG;
    Vec3 tibial_tray_Center_inG = tibial_tray->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 femoral_comp_Center_inG = femoral_component->findStationLocationInGround(*state, Vec3(0, 0, 0));
    Vec3 KneeCont_SumMoments_onFemoralComp_inG = -KneeCont_SumMoments_onTibialTray_inG + cross(tibial_tray_Center_inG - femoral_comp_Center_inG, KneeCont_SumForces_onFemoralComp_inG);

    // Compute residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndofr);
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
    /// Extract contact forces
    Vec3 AppliedPointForce_s1_l, AppliedPointForce_s2_l;
    Vec3 AppliedPointForce_s3_l, AppliedPointForce_s4_l;
    Vec3 AppliedPointForce_s5_l, AppliedPointForce_s6_l;
    Vec3 AppliedPointForce_s1_r, AppliedPointForce_s2_r;
    Vec3 AppliedPointForce_s3_r, AppliedPointForce_s4_r;
    Vec3 AppliedPointForce_s5_r, AppliedPointForce_s6_r;
    int nc = 3;
    for (int i = 0; i < nc; ++i) {
        AppliedPointForce_s1_l[i]   = up[i];
        AppliedPointForce_s2_l[i]   = up[i + nc];
        AppliedPointForce_s3_l[i]   = up[i + nc + nc];
        AppliedPointForce_s4_l[i]   = up[i + nc + nc + nc];
        AppliedPointForce_s5_l[i]   = up[i + nc + nc + nc + nc];
        AppliedPointForce_s6_l[i]   = up[i + nc + nc + nc + nc + nc];
        AppliedPointForce_s1_r[i]   = up[i + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s2_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s3_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s4_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s5_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s6_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
    }
    /// Extract contact sphere locations
    Vec3 locSphere_s1_r, locSphere_s2_r;
    Vec3 locSphere_s3_r, locSphere_s4_r;
    Vec3 locSphere_s5_r, locSphere_s6_r;
    /// Vertical positions are fixed
    locSphere_s1_r[1] = -0.021859;
    locSphere_s2_r[1] = -0.021859;
    locSphere_s3_r[1] = -0.021859;
    locSphere_s4_r[1] = -0.0214476;
    locSphere_s5_r[1] = -0.021859;
    locSphere_s6_r[1] = -0.0214476;
    int count = 0;
    for (int i = 0; i < nc; i+=2) {
        locSphere_s1_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        locSphere_s2_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1];
        locSphere_s3_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1];
        locSphere_s4_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1];
        locSphere_s5_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1];
        locSphere_s6_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1];
        ++count;
    }
    Vec3 locSphere_s1_l(locSphere_s1_r[0],locSphere_s1_r[1],-locSphere_s1_r[2]);
    Vec3 locSphere_s2_l(locSphere_s2_r[0],locSphere_s2_r[1],-locSphere_s2_r[2]);
    Vec3 locSphere_s3_l(locSphere_s3_r[0],locSphere_s3_r[1],-locSphere_s3_r[2]);
    Vec3 locSphere_s4_l(locSphere_s4_r[0],locSphere_s4_r[1],-locSphere_s4_r[2]);
    Vec3 locSphere_s5_l(locSphere_s5_r[0],locSphere_s5_r[1],-locSphere_s5_r[2]);
    Vec3 locSphere_s6_l(locSphere_s6_r[0],locSphere_s6_r[1],-locSphere_s6_r[2]);
    /// Extract radii
    osim_double_adouble radius_s1, radius_s2, radius_s3, radius_s4, radius_s5, radius_s6;
    radius_s1 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1];
    radius_s2 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 1];
    radius_s3 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 2];
    radius_s4 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 3];
    radius_s5 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 4];
    radius_s6 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 5];
    /// Compute contact point positions in body frames
    Vec3 normal = Vec3(0, 1, 0);
    /// sphere 1 left
    Vec3 pos_InGround_HC_s1_l = calcn_l->findStationLocationInGround(*state, locSphere_s1_l);
    Vec3 contactPointpos_InGround_HC_s1_l = pos_InGround_HC_s1_l - radius_s1*normal;
    Vec3 contactPointpos_InGround_HC_s1_l_adj = contactPointpos_InGround_HC_s1_l - 0.5*contactPointpos_InGround_HC_s1_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_l_adj, *calcn_l);
    /// sphere 2 left
    Vec3 pos_InGround_HC_s2_l = calcn_l->findStationLocationInGround(*state, locSphere_s2_l);
    Vec3 contactPointpos_InGround_HC_s2_l = pos_InGround_HC_s2_l - radius_s2*normal;
    Vec3 contactPointpos_InGround_HC_s2_l_adj = contactPointpos_InGround_HC_s2_l - 0.5*contactPointpos_InGround_HC_s2_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_l_adj, *calcn_l);
    /// sphere 3 left
    Vec3 pos_InGround_HC_s3_l = calcn_l->findStationLocationInGround(*state, locSphere_s3_l);
    Vec3 contactPointpos_InGround_HC_s3_l = pos_InGround_HC_s3_l - radius_s3*normal;
    Vec3 contactPointpos_InGround_HC_s3_l_adj = contactPointpos_InGround_HC_s3_l - 0.5*contactPointpos_InGround_HC_s3_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_l_adj, *calcn_l);
    /// sphere 4 left
    Vec3 pos_InGround_HC_s4_l = toes_l->findStationLocationInGround(*state, locSphere_s4_l);
    Vec3 contactPointpos_InGround_HC_s4_l = pos_InGround_HC_s4_l - radius_s4*normal;
    Vec3 contactPointpos_InGround_HC_s4_l_adj = contactPointpos_InGround_HC_s4_l - 0.5*contactPointpos_InGround_HC_s4_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_l_adj, *toes_l);
    /// sphere 5 left
    Vec3 pos_InGround_HC_s5_l = calcn_l->findStationLocationInGround(*state, locSphere_s5_l);
    Vec3 contactPointpos_InGround_HC_s5_l = pos_InGround_HC_s5_l - radius_s5*normal;
    Vec3 contactPointpos_InGround_HC_s5_l_adj = contactPointpos_InGround_HC_s5_l - 0.5*contactPointpos_InGround_HC_s5_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s5_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_l_adj, *calcn_l);
    /// sphere 6 left
    Vec3 pos_InGround_HC_s6_l = toes_l->findStationLocationInGround(*state, locSphere_s6_l);
    Vec3 contactPointpos_InGround_HC_s6_l = pos_InGround_HC_s6_l - radius_s6*normal;
    Vec3 contactPointpos_InGround_HC_s6_l_adj = contactPointpos_InGround_HC_s6_l - 0.5*contactPointpos_InGround_HC_s6_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s6_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_l_adj, *toes_l);
    /// sphere 1 right
    Vec3 pos_InGround_HC_s1_r = calcn_r->findStationLocationInGround(*state, locSphere_s1_r);
    Vec3 contactPointpos_InGround_HC_s1_r = pos_InGround_HC_s1_r - radius_s1*normal;
    Vec3 contactPointpos_InGround_HC_s1_r_adj = contactPointpos_InGround_HC_s1_r - 0.5*contactPointpos_InGround_HC_s1_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_r_adj, *calcn_r);
    /// sphere 2 right
    Vec3 pos_InGround_HC_s2_r = calcn_r->findStationLocationInGround(*state, locSphere_s2_r);
    Vec3 contactPointpos_InGround_HC_s2_r = pos_InGround_HC_s2_r - radius_s2*normal;
    Vec3 contactPointpos_InGround_HC_s2_r_adj = contactPointpos_InGround_HC_s2_r - 0.5*contactPointpos_InGround_HC_s2_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_r_adj, *calcn_r);
    /// sphere 3 right
    Vec3 pos_InGround_HC_s3_r = calcn_r->findStationLocationInGround(*state, locSphere_s3_r);
    Vec3 contactPointpos_InGround_HC_s3_r = pos_InGround_HC_s3_r - radius_s3*normal;
    Vec3 contactPointpos_InGround_HC_s3_r_adj = contactPointpos_InGround_HC_s3_r - 0.5*contactPointpos_InGround_HC_s3_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_r_adj, *calcn_r);
    /// sphere 4 right
    Vec3 pos_InGround_HC_s4_r = toes_r->findStationLocationInGround(*state, locSphere_s4_r);
    Vec3 contactPointpos_InGround_HC_s4_r = pos_InGround_HC_s4_r - radius_s4*normal;
    Vec3 contactPointpos_InGround_HC_s4_r_adj = contactPointpos_InGround_HC_s4_r - 0.5*contactPointpos_InGround_HC_s4_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_r_adj, *toes_r);
    /// sphere 5 right
    Vec3 pos_InGround_HC_s5_r = calcn_r->findStationLocationInGround(*state, locSphere_s5_r);
    Vec3 contactPointpos_InGround_HC_s5_r = pos_InGround_HC_s5_r - radius_s5*normal;
    Vec3 contactPointpos_InGround_HC_s5_r_adj = contactPointpos_InGround_HC_s5_r - 0.5*contactPointpos_InGround_HC_s5_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s5_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_r_adj, *calcn_r);
    /// sphere 6 right
    Vec3 pos_InGround_HC_s6_r = toes_r->findStationLocationInGround(*state, locSphere_s6_r);
    Vec3 contactPointpos_InGround_HC_s6_r = pos_InGround_HC_s6_r - radius_s6*normal;
    Vec3 contactPointpos_InGround_HC_s6_r_adj = contactPointpos_InGround_HC_s6_r - 0.5*contactPointpos_InGround_HC_s6_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s6_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_r_adj, *toes_r);
    /// Add GRF contact forces to appliedBodyForces
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s1_l, AppliedPointForce_s1_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s2_l, AppliedPointForce_s2_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s3_l, AppliedPointForce_s3_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s4_l, AppliedPointForce_s4_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s5_l, AppliedPointForce_s5_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s6_l, AppliedPointForce_s6_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s1_r, AppliedPointForce_s1_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s2_r, AppliedPointForce_s2_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s3_r, AppliedPointForce_s3_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s4_r, AppliedPointForce_s4_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s5_r, AppliedPointForce_s5_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s6_r, AppliedPointForce_s6_r, appliedBodyForces);
    /// Add knee contact forces to appliedBodyForces
    model->getMatterSubsystem().addInStationForce(*state, tibial_tray->getMobilizedBodyIndex(), Vec3(0, 0, 0), KneeCont_SumForces_onTibialTray_inG, appliedBodyForces);
    model->getMatterSubsystem().addInBodyTorque(*state, tibial_tray->getMobilizedBodyIndex(), KneeCont_SumMoments_onTibialTray_inG, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, femoral_component->getMobilizedBodyIndex(), Vec3(0, 0, 0), KneeCont_SumForces_onFemoralComp_inG, appliedBodyForces);
    model->getMatterSubsystem().addInBodyTorque(*state, femoral_component->getMobilizedBodyIndex(), KneeCont_SumMoments_onFemoralComp_inG, appliedBodyForces);

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
    Vector knownUdot(ndofr);
    knownUdot.setToZero();
    for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];
    ///  Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Compute contact torques about the ground frame origin
    /// Get transforms
    SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);
    /// Calculate torques
    Vec3 AppliedPointTorque_s1_l, AppliedPointTorque_s2_l, AppliedPointTorque_s3_l, AppliedPointTorque_s4_l, AppliedPointTorque_s5_l, AppliedPointTorque_s6_l;
    Vec3 AppliedPointTorque_s1_r, AppliedPointTorque_s2_r, AppliedPointTorque_s3_r, AppliedPointTorque_s4_r, AppliedPointTorque_s5_r, AppliedPointTorque_s6_r;
    AppliedPointTorque_s1_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s1_l) % AppliedPointForce_s1_l;
    AppliedPointTorque_s2_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s2_l) % AppliedPointForce_s2_l;
    AppliedPointTorque_s3_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s3_l) % AppliedPointForce_s3_l;
    AppliedPointTorque_s4_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s4_l) % AppliedPointForce_s4_l;
    AppliedPointTorque_s5_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s5_l) % AppliedPointForce_s5_l;
    AppliedPointTorque_s6_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s6_l) % AppliedPointForce_s6_l;
    AppliedPointTorque_s1_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s1_r) % AppliedPointForce_s1_r;
    AppliedPointTorque_s2_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s2_r) % AppliedPointForce_s2_r;
    AppliedPointTorque_s3_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s3_r) % AppliedPointForce_s3_r;
    AppliedPointTorque_s4_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s4_r) % AppliedPointForce_s4_r;
    AppliedPointTorque_s5_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s5_r) % AppliedPointForce_s5_r;
    AppliedPointTorque_s6_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s6_r) % AppliedPointForce_s6_r;
    /// Contact torques
    Vec3 MOM_l, MOM_r;
    MOM_l = AppliedPointTorque_s1_l + AppliedPointTorque_s2_l + AppliedPointTorque_s3_l + AppliedPointTorque_s4_l + AppliedPointTorque_s5_l + AppliedPointTorque_s6_l;
    MOM_r = AppliedPointTorque_s1_r + AppliedPointTorque_s2_r + AppliedPointTorque_s3_r + AppliedPointTorque_s4_r + AppliedPointTorque_s5_r + AppliedPointTorque_s6_r;
    /// Contact forces
    Vec3 GRF_r = AppliedPointForce_s1_r + AppliedPointForce_s2_r + AppliedPointForce_s3_r + AppliedPointForce_s4_r + AppliedPointForce_s5_r + AppliedPointForce_s6_r;
    Vec3 GRF_l = AppliedPointForce_s1_l + AppliedPointForce_s2_l + AppliedPointForce_s3_l + AppliedPointForce_s4_l + AppliedPointForce_s5_l + AppliedPointForce_s6_l;

    // Extract results
    /// Residual forces
    for (int i = 0; i < 12; ++i) {
        res[0][i] = value<T>(residualMobilityForces[i]);
    }
    /// OpenSim and Simbody have different state orders so we adjust manually
    res[0][12] = value<T>(residualMobilityForces[15]);
    res[0][13] = value<T>(residualMobilityForces[22]);
    res[0][14] = value<T>(residualMobilityForces[23]);
    res[0][15] = value<T>(residualMobilityForces[24]);
    res[0][16] = value<T>(residualMobilityForces[25]);
    res[0][17] = value<T>(residualMobilityForces[26]);
    res[0][18] = value<T>(residualMobilityForces[27]);
    res[0][19] = value<T>(residualMobilityForces[28]);
    res[0][20] = value<T>(residualMobilityForces[34]);
    res[0][21] = value<T>(residualMobilityForces[31]);
    res[0][22] = value<T>(residualMobilityForces[35]);
    res[0][23] = value<T>(residualMobilityForces[12]);
    res[0][24] = value<T>(residualMobilityForces[13]);
    res[0][25] = value<T>(residualMobilityForces[14]);
    res[0][26] = value<T>(residualMobilityForces[16]);
    res[0][27] = value<T>(residualMobilityForces[17]);
    res[0][28] = value<T>(residualMobilityForces[18]);
    res[0][29] = value<T>(residualMobilityForces[19]);
    res[0][30] = value<T>(residualMobilityForces[20]);
    res[0][31] = value<T>(residualMobilityForces[21]);
    res[0][32] = value<T>(residualMobilityForces[29]);
    res[0][33] = value<T>(residualMobilityForces[30]);

    /// Contact forces
    for (int i = 0; i < nc; ++i) {
        res[0][i + ndof] = value<T>(GRF_r[i]);      /// GRF_r
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i + ndof + nc] = value<T>(GRF_l[i]); /// GRF_l
    }
    /// Contact torques
    for (int i = 0; i < nc; ++i) {
        res[0][i + ndof + nc + nc] = value<T>(MOM_r[i]);        /// GRM_r
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i + ndof + nc + nc + nc] = value<T>(MOM_l[i]);   /// GRM_l
    }
    
  /*  /// Knee contact forces
    res[0][ndof + nc + nc + nc + nc] = value<T>(SumForces_vert_Lat);
    res[0][ndof + nc + nc + nc + nc + 1] = value<T>(SumForces_vert_Med);
    std::cout << "SumForces_vert_Lat" << SumForces_vert_Lat << std::endl;
    std::cout << "SumForces_vert_Med" << SumForces_vert_Med << std::endl;*/
    std::cout << "computed until the end" << std::endl;

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
    Recorder KCF_in[NKCFM];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;
    for (int i = 0; i < NP; ++i) p[i] <<= 0;
    for (int i = 0; i < NKCFM; ++i) KCF_in[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u,p,KCF_in };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}

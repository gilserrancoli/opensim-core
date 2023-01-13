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
constexpr int ndof = 34;        // # degrees of freedom (excluding locked)
constexpr int ndofr = 36;       // # degrees of freedom (including locked)
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NP = 54;          // # parameters
constexpr int NR = ndof+6+6+2;    // # residual torques + # GRFs + # GRMs

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

std::vector<std::vector<int>> ReadDataIntx3columns(std::string filename, int nrows_in) {


    std::ifstream  file(filename);

    //
    int v1, v2, v3;
    std::vector<std::vector<int>> out_csv(nrows_in, std::vector<int>(3, 0));

    int nrows = 0;
    if (file.is_open()) {

        std::string line;

        while (!file.eof()) {

            std::getline(file, line);

            std::istringstream iss(line);

            std::cout << line << std::endl;

            std::string delimiter = ",";
            std::string token = line.substr(0, line.find(","));

            size_t pos0 = 0;
            pos0 = line.find(",");
            std::string token1 = line.substr(0, pos0);
            v1 = std::stoi(token1);
            line.erase(0, pos0 + delimiter.length());
            size_t pos1 = line.find(",");
            std::string token2 = line.substr(0, pos1);
            v2 = std::stoi(token2);
            line.erase(0, pos1 + delimiter.length());
            size_t pos2 = line.find(",");
            std::string token3 = line.substr(0, pos2);
            v3 = std::stoi(token3);

            out_csv[nrows] = { v1, v2, v3 };
            nrows = nrows + 1;

        }

        file.close();
    }
    return out_csv;
}

std::vector<std::vector<int>> ReadDataIntx2columns(std::string filename, int nrows_in) {


    std::ifstream  file(filename);

    //
    int v1, v2;
    std::vector<std::vector<int>> out_csv(nrows_in,std::vector<int>(2,0));
   
    int nrows = 0;
    if (file.is_open()) {

        std::string line;

        while (!file.eof()) {

            std::getline(file, line);

            std::istringstream iss(line);

            std::cout << line << std::endl;
            std::string delimiter = ",";
            std::string token = line.substr(0, line.find(","));

            
            size_t pos0 = 0;
            pos0 = line.find(",");
            std::string token1 = line.substr(0, pos0);
            v1 = std::stoi(token1);
            line.erase(0, pos0 + delimiter.length());
            size_t pos1 = line.find(",");
            std::string token2 = line.substr(0, pos1);
            v2 = std::stoi(token2);
            
            out_csv[nrows] = { v1, v2 };

            nrows = nrows + 1;           
        }
        
        file.close();
    }
    return out_csv;
}

Vector_<Vec3> ReadDataDoublex3columns(std::string filename, int nrows_in) {

    
    std::ifstream  file(filename);

    //
    double v1, v2, v3;
    Vector_<Vec3> out_csv(nrows_in);
    
    int nrows = 0;
    if (file.is_open()) {

        std::string line;
        
        while (!file.eof()) {

            std::getline(file, line);

            std::istringstream iss(line);

            std::string delimiter = ",";
            std::string token = line.substr(0, line.find(","));

            size_t pos0 = 0;
            pos0 = line.find(",");
            std::string token1 = line.substr(0, pos0);
            v1 = std::stod(token1);
            line.erase(0, pos0 + delimiter.length());
            size_t pos1 = line.find(",");
            std::string token2 = line.substr(0, pos1);
            v2 = std::stod(token2);
            line.erase(0, pos1 + delimiter.length());
            size_t pos2 = line.find(",");
            std::string token3 = line.substr(0, pos2);
            v3 = std::stod(token3);

            out_csv[nrows] = Vec3(v1, v2, v3);
            nrows = nrows + 1;

        }
        
        file.close();
    }
  

    Vector_<Vec3> tibPoints(36);
    tibPoints[0] = Vec3(0.0140, 0.0087, 0.0079);
    tibPoints[1] = Vec3(0.0037, 0.0105, 0.0065);
    tibPoints[2] = Vec3(0.0172, 0.0086, 0.0129);
    tibPoints[3] = Vec3(-0.0069, 0.0117, 0.0059);
    tibPoints[4] = Vec3(0.0052, 0.0060, 0.0181);
    tibPoints[5] = Vec3(-0.0055, 0.0110, -0.0073);
    tibPoints[6] = Vec3(0.0112, 0.0091, -0.0090);
    tibPoints[7] = Vec3(-0.0064, 0.0078, -0.0153);
    tibPoints[8] = Vec3(0.0173, 0.0067, 0.0199);
    tibPoints[9] = Vec3(-0.0076, 0.0072, 0.0174);
    tibPoints[10] = Vec3(0.0019, 0.0061, 0.0256);
    tibPoints[11] = Vec3(-0.0192, 0.0117, 0.0165);
    tibPoints[12] = Vec3(-0.0114, 0.0079, 0.0258);
    tibPoints[13] = Vec3(-0.0206, 0.0120, 0.0245);
    tibPoints[14] = Vec3(-0.0131, 0.0118, 0.0342);
    tibPoints[15] = Vec3(-0.0064, 0.0092, 0.0325);
    tibPoints[16] = Vec3(0.0161, 0.0069, 0.0273);
    tibPoints[17] = Vec3(0.0106, 0.0087, 0.0333);
    tibPoints[18] = Vec3(0.0010, 0.0118, 0.0373);
    tibPoints[19] = Vec3(0.0253, 0.0091, 0.0184);
    tibPoints[20] = Vec3(0.0212, 0.0089, 0.0300);
    tibPoints[21] = Vec3(0.0122, 0.0087, 0.0382);
    tibPoints[22] = Vec3(0.0198, 0.0089, -0.0095);
    tibPoints[23] = Vec3(0.0167, 0.0067, -0.0190);
    tibPoints[24] = Vec3(0.0231, 0.0089, -0.0299);
    tibPoints[25] = Vec3(0.0165, 0.0068, -0.0260);
    tibPoints[26] = Vec3(0.0122, 0.0088, -0.0327);
    tibPoints[27] = Vec3(0.0065, 0.0062, -0.0176);
    tibPoints[28] = Vec3(-0.0077, 0.0068, -0.0225);
    tibPoints[29] = Vec3(0.0058, 0.0064, -0.0279);
    tibPoints[30] = Vec3(-0.0115, 0.0119, -0.0351);
    tibPoints[31] = Vec3(-0.0198, 0.0113, -0.0262);
    tibPoints[32] = Vec3(-0.0071, 0.0084, -0.0305);
    tibPoints[33] = Vec3(-0.0174, 0.0113, -0.0146);
    tibPoints[34] = Vec3(0.0112, 0.0089, -0.0380);
    tibPoints[35] = Vec3(0.0018, 0.0116, -0.0371);

    Vector_<Vec3> facesTib(49);
    facesTib[0] = Vec3(1, 2, 3);
    facesTib[1] = Vec3(1, 3, 20);
    facesTib[2] = Vec3(2, 4, 5);
    facesTib[3] = Vec3(3, 2, 5);
    facesTib[4] = Vec3(3, 5, 9);
    facesTib[5] = Vec3(4, 10, 5);
    facesTib[6] = Vec3(5, 10, 11);
    facesTib[7] = Vec3(5, 11, 17);
    facesTib[8] = Vec3(6, 7, 8);
    facesTib[9] = Vec3(7, 24, 28);
    facesTib[10] = Vec3(8, 7, 28);
    facesTib[11] = Vec3(8, 34, 6);
    facesTib[12] = Vec3(9, 5, 17);
    facesTib[13] = Vec3(9, 17, 20);
    facesTib[14] = Vec3(10, 12, 13);
    facesTib[15] = Vec3(10, 13, 11);
    facesTib[16] = Vec3(11, 16, 18);
    facesTib[17] = Vec3(12, 10, 4);
    facesTib[18] = Vec3(13, 15, 16);
    facesTib[19] = Vec3(13, 14, 15);
    facesTib[20] = Vec3(13, 16, 11);
    facesTib[21] = Vec3(14, 13, 12);
    facesTib[22] = Vec3(16, 15, 19);
    facesTib[23] = Vec3(17, 11, 18);
    facesTib[24] = Vec3(17, 21, 20);
    facesTib[25] = Vec3(18, 16, 19);
    facesTib[26] = Vec3(18, 19, 22);
    facesTib[27] = Vec3(20, 3, 9);
    facesTib[28] = Vec3(21, 17, 18);
    facesTib[29] = Vec3(21, 18, 22);
    facesTib[30] = Vec3(23, 24, 7);
    facesTib[31] = Vec3(23, 25, 24);
    facesTib[32] = Vec3(25, 26, 24);
    facesTib[33] = Vec3(25, 27, 26);
    facesTib[34] = Vec3(26, 27, 30);
    facesTib[35] = Vec3(27, 25, 35);
    facesTib[36] = Vec3(27, 35, 36);
    facesTib[37] = Vec3(28, 29, 8);
    facesTib[38] = Vec3(28, 24, 30);
    facesTib[39] = Vec3(29, 28, 30);
    facesTib[40] = Vec3(30, 33, 29);
    facesTib[41] = Vec3(30, 24, 26);
    facesTib[42] = Vec3(31, 32, 33);
    facesTib[43] = Vec3(31, 33, 36);
    facesTib[44] = Vec3(32, 34, 29);
    facesTib[45] = Vec3(33, 32, 29);
    facesTib[46] = Vec3(34, 8, 29);
    facesTib[47] = Vec3(36, 33, 30);
    facesTib[48] = Vec3(36, 30, 27);

    Vector_<Vec3> femPoints(117);
    femPoints[0] = Vec3(-0.0027, -0.0321, 0.0111);
    femPoints[1] = Vec3(0.0144, -0.0332, 0.0130);
    femPoints[2] = Vec3(0.0055, -0.0364, 0.0176);
    femPoints[3] = Vec3(-0.0302, -0.0052, 0.0226);
    femPoints[4] = Vec3(-0.0334, -0.0053, 0.0172);
    femPoints[5] = Vec3(-0.0270, -0.0164, 0.0267);
    femPoints[6] = Vec3(0.0333, -0.0140, 0.0310);
    femPoints[7] = Vec3(0.0276, -0.0193, 0.0363);
    femPoints[8] = Vec3(0.0277, -0.0257, 0.0311);
    femPoints[9] = Vec3(0.0289, -0.0075, 0.0364);
    femPoints[10] = Vec3(0.0324, -0.0017, 0.0264);
    femPoints[11] = Vec3(0.0277, 0.0001, 0.0328);
    femPoints[12] = Vec3(0.0345, -0.0075, 0.0252);
    femPoints[13] = Vec3(0.0322, -0.0227, 0.0250);
    femPoints[14] = Vec3(0.0259, -0.0299, 0.0200);
    femPoints[15] = Vec3(0.0350, -0.0129, 0.0186);
    femPoints[16] = Vec3(0.0289, 0.0022, 0.0159);
    femPoints[17] = Vec3(0.0280, 0.0040, 0.0253);
    femPoints[18] = Vec3(0.0331, -0.0038, 0.0169);
    femPoints[19] = Vec3(0.0333, -0.0199, 0.0174);
    femPoints[20] = Vec3(0.0272, -0.0225, 0.0093);
    femPoints[21] = Vec3(0.0307, -0.0103, 0.0093);
    femPoints[22] = Vec3(0.0297, -0.0177, 0.0093);
    femPoints[23] = Vec3(0.0223, -0.0292, 0.0126);
    femPoints[24] = Vec3(0.0129, -0.0311, 0.0091);
    femPoints[25] = Vec3(0.0200, -0.0318, 0.0294);
    femPoints[26] = Vec3(0.0203, -0.0265, 0.0365);
    femPoints[27] = Vec3(-0.0070, -0.0306, 0.0331);
    femPoints[28] = Vec3(-0.0151, -0.0291, 0.0303);
    femPoints[29] = Vec3(-0.0039, -0.0341, 0.0311);
    femPoints[30] = Vec3(-0.0236, -0.0227, 0.0271);
    femPoints[31] = Vec3(-0.0231, -0.0260, 0.0217);
    femPoints[32] = Vec3(0.0075, -0.0312, 0.0364);
    femPoints[33] = Vec3(-0.0325, 0.0105, 0.0150);
    femPoints[34] = Vec3(-0.0315, 0.0188, 0.0079);
    femPoints[35] = Vec3(-0.0330, 0.0191, 0.0022);
    femPoints[36] = Vec3(-0.0351, 0.0033, 0.0073);
    femPoints[37] = Vec3(-0.0339, -0.0066, 0.0096);
    femPoints[38] = Vec3(-0.0346, 0.0080, 0.0108);
    femPoints[39] = Vec3(-0.0319, 0.0125, -0.0033);
    femPoints[40] = Vec3(-0.0325, 0.0031, 0.0012);
    femPoints[41] = Vec3(-0.0323, -0.0007, -0.0147);
    femPoints[42] = Vec3(-0.0338, 0.0159, -0.0200);
    femPoints[43] = Vec3(-0.0351, 0.0019, -0.0207);
    femPoints[44] = Vec3(-0.0308, 0.0216, -0.0146);
    femPoints[45] = Vec3(-0.0316, 0.0060, -0.0078);
    femPoints[46] = Vec3(-0.0306, 0.0203, -0.0071);
    femPoints[47] = Vec3(-0.0303, -0.0054, -0.0094);
    femPoints[48] = Vec3(-0.0307, -0.0047, 0.0001);
    femPoints[49] = Vec3(-0.0343, -0.0008, -0.0261);
    femPoints[50] = Vec3(-0.0313, -0.0100, -0.0295);
    femPoints[51] = Vec3(-0.0330, -0.0106, -0.0200);
    femPoints[52] = Vec3(-0.0341, 0.0122, -0.0245);
    femPoints[53] = Vec3(-0.0324, 0.0043, -0.0280);
    femPoints[54] = Vec3(-0.0294, -0.0075, -0.0306);
    femPoints[55] = Vec3(-0.0142, -0.0318, -0.0168);
    femPoints[56] = Vec3(-0.0037, -0.0355, -0.0171);
    femPoints[57] = Vec3(-0.0060, -0.0306, -0.0089);
    femPoints[58] = Vec3(-0.0033, -0.0352, -0.0284);
    femPoints[59] = Vec3(-0.0104, -0.0322, -0.0295);
    femPoints[60] = Vec3(-0.0057, -0.0310, -0.0346);
    femPoints[61] = Vec3(-0.0274, -0.0208, -0.0243);
    femPoints[62] = Vec3(-0.0237, -0.0214, -0.0308);
    femPoints[63] = Vec3(-0.0198, -0.0273, -0.0277);
    femPoints[64] = Vec3(-0.0109, -0.0340, -0.0225);
    femPoints[65] = Vec3(-0.0216, -0.0274, -0.0180);
    femPoints[66] = Vec3(-0.0301, -0.0156, -0.0172);
    femPoints[67] = Vec3(-0.0161, -0.0268, -0.0325);
    femPoints[68] = Vec3(0.0035, -0.0367, -0.0237);
    femPoints[69] = Vec3(-0.0273, -0.0133, -0.0088);
    femPoints[70] = Vec3(-0.0240, -0.0187, -0.0105);
    femPoints[71] = Vec3(-0.0250, -0.0158, -0.0021);
    femPoints[72] = Vec3(-0.0263, -0.0168, 0.0064);
    femPoints[73] = Vec3(-0.0303, -0.0159, 0.0122);
    femPoints[74] = Vec3(-0.0236, -0.0251, 0.0153);
    femPoints[75] = Vec3(-0.0143, -0.0277, -0.0103);
    femPoints[76] = Vec3(-0.0179, -0.0227, 0.0004);
    femPoints[77] = Vec3(-0.0164, -0.0264, 0.0085);
    femPoints[78] = Vec3(-0.0110, -0.0273, -0.0021);
    femPoints[79] = Vec3(-0.0103, -0.0287, 0.0064);
    femPoints[80] = Vec3(-0.0158, -0.0309, 0.0157);
    femPoints[81] = Vec3(-0.0082, -0.0334, 0.0147);
    femPoints[82] = Vec3(-0.0034, -0.0359, 0.0195);
    femPoints[83] = Vec3(-0.0100, -0.0341, 0.0246);
    femPoints[84] = Vec3(-0.0302, -0.0154, 0.0203);
    femPoints[85] = Vec3(0.0027, -0.0363, 0.0269);
    femPoints[86] = Vec3(0.0137, -0.0356, 0.0239);
    femPoints[87] = Vec3(0.0273, -0.0082, 0.0355);
    femPoints[88] = Vec3(-0.0220 - 0.0206, -0.0319);
    femPoints[89] = Vec3(0.0328, -0.0189, -0.0150);
    femPoints[90] = Vec3(0.0244, -0.0253, -0.0092);
    femPoints[91] = Vec3(0.0279, -0.0264, -0.0149);
    femPoints[92] = Vec3(0.0294, -0.0005, -0.0127);
    femPoints[93] = Vec3(0.0303, -0.0111, -0.0086);
    femPoints[94] = Vec3(0.0333, -0.0095, -0.0139);
    femPoints[95] = Vec3(0.0351, -0.0121, -0.0193);
    femPoints[96] = Vec3(0.0291, 0.0036, -0.0212);
    femPoints[97] = Vec3(0.0333, -0.0030, -0.0207);
    femPoints[98] = Vec3(0.0310, -0.0249, -0.0227);
    femPoints[99] = Vec3(0.0256, -0.0302, -0.0242);
    femPoints[100] = Vec3(0.0344, -0.0162, -0.0259);
    femPoints[101] = Vec3(0.0333, -0.0047, -0.0278);
    femPoints[102] = Vec3(0.0283, 0.0004, -0.0323);
    femPoints[103] = Vec3(0.0315, -0.0200, -0.0306);
    femPoints[104] = Vec3(0.0265, -0.0262, -0.0320);
    femPoints[105] = Vec3(0.0303, -0.0104, -0.0359);
    femPoints[106] = Vec3(0.0270, -0.0208, -0.0364);
    femPoints[107] = Vec3(0.0217, -0.0311, -0.0158);
    femPoints[108] = Vec3(0.0152, -0.0318, -0.0116);
    femPoints[109] = Vec3(0.0160, -0.0348, -0.0216);
    femPoints[110] = Vec3(0.0061, -0.0360, -0.0168);
    femPoints[111] = Vec3(0.0105, -0.0361, -0.0254);
    femPoints[112] = Vec3(0.0189, -0.0317, -0.0304);
    femPoints[113] = Vec3(0.0047, -0.0335, -0.0335);
    femPoints[114] = Vec3(0.0175, -0.0289, -0.0352);
    femPoints[115] = Vec3(-0.0005, -0.0331, -0.0123);
    femPoints[116] = Vec3(0.0271, -0.0061, -0.0091);
    femPoints = femPoints + Vec3(0, 0.042, 0);

    /*Vector_<Vec3> facesFem(188);
    faceFem[0] = Vec3(1     2     3
    faceFem[0] = 1    82    78
    faceFem[0] = 2    87     3
    faceFem[0] = 3    83    82
    faceFem[0] = 3    82     1
    faceFem[0] = 3    86    83
    faceFem[0] = 4     5     6
    faceFem[0] = 5     4    34
    faceFem[0] = 5    38    85
    faceFem[0] = 6     5    85
    faceFem[0] = 7     8     9
    faceFem[0] = 7     9    14
    faceFem[0] = 8     7    10
    faceFem[0] = 9     8    27
    faceFem[0] = 10    11    12
    faceFem[0] = 10    12    88
    faceFem[0] = 10    88     8
    faceFem[0] = 11    10     7
    faceFem[0] = 11    13    19
    faceFem[0] = 12    11    18
    faceFem[0] = 12    18    88
    faceFem[0] = 13    11     7
    faceFem[0] = 13     7    16
    faceFem[0] = 14     9    15
    faceFem[0] = 14    15    20
    faceFem[0] = 15     9    26
    faceFem[0] = 15    26    87
    faceFem[0] = 15    87     2
    16     7    14
    16    14    20
    16    20    23
    17    18    11
    17    11    19
    17    19    22
    19    13    16
    19    16    22
    20    15    21
    20    21    23
    21    15    24
    21    24    25
    22    16    23
    24    15     2
    24     2    25
    26    27    33
    26    33    87
    27    26     9
    28     6    29
    28    29    30
    29    31    32
    29     6    31
    29    32    84
    29    84    30
    30    33    28
    31     6    85
    32    74    75
    32    81    84
    34    39     5
    35    36    34
    36    39    34
    37    38     5
    37    36    40
    37    40    41
    38    37    41
    38    41    49
    38    49    74
    39    37     5
    39    36    37
    40    45    46
    40    46    41
    40    36    47
    41    46    49
    42    43    44
    42    52    67
    42    67    48
    44    50    52
    44    43    53
    45    43    46
    46    43    42
    46    42    48
    47    45    40
    48    67    70
    48    70    49
    48    49    46
    49    70    72
    49    72    73
    49    73    74
    50    51    52
    50    44    53
    50    53    54
    50    54    51
    51    62    52
    51    55    89
    52    42    44
    52    62    67
    55    51    54
    56    57    58
    56    65    57
    56    58    76
    57    59    69
    57   116    58
    59    60    61
    61    68    89
    61    60    68
    61   114    59
    62    63    64
    62    64    66
    63    68    64
    63    62    51
    64    60    65
    65    60    59
    65    59    57
    66    64    65
    66    65    56
    66    67    62
    66    56    76
    66    76    71
    67    66    71
    68    60    64
    68    63    89
    69    59   114
    70    67    71
    70    71    72
    71    76    72
    72    77    73
    73    77    78
    74    73    75
    75    81    32
    76    77    72
    76    58    79
    77    79    80
    77    76    79
    77    80    78
    78    75    73
    78    81    75
    78    82    81
    80     1    78
    81    82    83
    81    83    84
    84    83    86
    84    86    30
    85    38    74
    85    74    32
    85    32    31
    86    33    30
    87    33    86
    87    86     3
    89    63    51
    90    91    92
    90    92    99
    90    95    94
    91   109   108
    93    94    95
    93    95    98
    94    91    90
    95    90    96
    96    90   101
    97    93    98
    97    98   102
    98    95    96
    98    96   102
    99    92   100
    100    92   108
    100   110   113
    100   113   105
    101    90    99
    101    99   104
    102    96   101
    102   103    97
    102   101   106
    103   102   106
    104    99   105
    105    99   100
    105   113   115
    105   115   107
    106   101   104
    106   104   107
    107   104   105
    108    92    91
    108   109   110
    108   110   100
    110   109   111
    110   111   112
    111    57    69
    111   116    57
    112   111    69
    112    69   114
    113   110   112
    117    94    93);*/
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

void CalculateIntersection(Vector_<Vec3> fem_Points, Vector_<Vec3> tib_Points, Vec3 &d, Real &As, Real &Atib, Vec3 &ns, Vec3 &Cfem, Vec3 &Ctib, Vec3 &nt) {
    // mean point of the tib face
    for (int i = 0; i < 3; i++) {
        Ctib[i] = mean (Vec3(tib_Points[0][i], tib_Points[1][i], tib_Points[2][i]));
    }

    //edges of tib face
    Vec3 edge1_p = tib_Points[1] - tib_Points[0];
    Vec3 edge2_p = tib_Points[2] - tib_Points[0];

    // unit vector normal to tibia
    nt = cross(edge1_p, edge2_p);
    nt = nt.normalize();

    // Area of tib face
    Real anglec = acos(dot(edge1_p.normalize(), edge2_p.normalize()));
    Atib = 0.5*edge1_p.norm()*edge2_p.norm()*sin(anglec);
    
    // mean point of the fem face
    for (int i = 0; i < 3; i++) {
        Cfem[i] = mean(Vec3(fem_Points[0][i], fem_Points[1][i], fem_Points[2][i]));
    }

    //edges of fem face
    Vec3 edge1_s = fem_Points[1] - fem_Points[0];
    Vec3 edge2_s = fem_Points[2] - fem_Points[0];
    //// unit vector normal to femur
    //ns = cross(edge1_s, edge2_s);
    //ns = ns.normalize();
    Real angles = acos(dot(edge1_s.normalize(),edge2_s.normalize()));
    As = (1 / 2)*edge1_s.norm()*edge2_s.norm()*sin(angles);

    d = Cfem - Ctib;

}

void CalculateMaximumPenetration(Vector_<Vec3> d_v, Vector_<Vec3> nt_v, Real &mindist, Vec3 &nt_l, Vector_<Real> highligh_index) {
    Vector_<Real> proj(d_v.size());
    //Vector_<Vec3> dist_v(d_v.size());
    Vector_<Real> pen(d_v.size());

    for (int i = 0; i < d_v.size(); i++) {
        proj[i] = dot(d_v[i], nt_v[i]); // projection of distance between tibial and femoral faces to the normal of tibial face
        pen[i] = -proj[i]; // pen is for penetration
    }
    Real k = 1e4;

    mindist = (log(sum(exp(k*pen))) / k);
    nt_l = nt_v[0];
    for (int i = 0; i < d_v.size(); i++) {
        highligh_index[i] = 1.0 / (mindist - pen[i]);
    }
}

Real CalculatePressure(Real poisson, Real E, Real d, Real h) {
    Real k = 1e4;
    Real pen = d;

    Real p_init = ((1 - poisson)*E / ((1 + poisson)*(1 - 2 * poisson)))*pen / h;

    Real p = p_init*(1 + tanh(k*pen)) / 2;

    return p;
}

void CalculateForceCompartment(Vector_<Vec3> femPoints, Vector_<Vec3> v_femPoints_wrtTibia, std::vector<std::vector<int>> facesFem, Vector_<Vec3> tibPoints_transf, std::vector<std::vector<int>> facesTib, std::vector<std::vector<int>> pairs_list, Vec3 &SumForces, Vec3 &SumMoments, Real poisson, Real E, Real h, Vec3 originTib_G, Vec3 knee_trans, Vec3 knee_rot) {
    Vector_<Vec3> d(pairs_list.size());
    Vector_<Vec3> nt(pairs_list.size());
    Vector_<Vec3> force_s_l(0);
    Vector_<Vec3> Ctib(0);
    Vector_<Real> face_s(0);
    Vector_<Vec3> mom_O(0);
    Vector_<Vec3> v_rel_femur_wrtTibia(0);

    int k = 1; //number of contacting elements in femur
    int l = 1; //count number of elements contacting element k of femur
    
    //std::cout << "facesFem.size()" << facesFem.size() << std::endl;
    //for (int i = 0; i < facesFem.size(); i++) {
    //    std::cout << "faces_fem[" << i << "]=" << facesFem.at(i).at(0) << " " << facesFem.at(i).at(1) << facesFem.at(i).at(2) << std::endl;
    //}
    
    Vector_<Vec3> fem_v_ind(pairs_list.size(), Vec3(0));
    for (int i = 0; i < pairs_list.size(); i++) {
        std::cout << pairs_list[i][1] - 1 << std::endl;
        //get indices of points for each face. Faces with the same order as in pairs_list
        std::vector<int> facesFem_ind = facesFem[pairs_list[i][1]-1];
        // get the points that configure each face of facesFem_ind, 
        Vector_<Vec3> fem_Points_ind(3);
        Vec3 fem_v_ind_aux(0);
        for (int j = 0; j < 3; j++) {
            fem_Points_ind[j] = femPoints[facesFem_ind[j]-1];
            fem_v_ind_aux = fem_v_ind_aux + v_femPoints_wrtTibia[facesFem_ind[j] - 1];
        }
        //get the velocities of the middle point of each femoral face wrt tibia, faces in the same order as facesFem_ind
        for (int j = 0; j < 3; j++) {
            fem_v_ind[i][j] = fem_v_ind_aux[j] / 3.0; // Division by three because we average over the three nodes
        }

        // get indices of facesTib for each pair, with the same order as in pairs_list
        std::vector<int> facesTib_ind = facesTib[pairs_list[i][0]-1];
        // get the points that configure each face of facesTib_ind
        Vector_<Vec3> tib_Points_ind(3);
        for (int j = 0; j < 3; j++) {
            tib_Points_ind[j] = tibPoints_transf[facesTib_ind[j]-1];
        }
        Vec3 d_aux, ns_aux, Cfem, Ctib_aux, nt_aux;
        Real As, At;
        
        CalculateIntersection(fem_Points_ind, tib_Points_ind, d_aux, As, At, ns_aux, Cfem, Ctib_aux, nt_aux); // first element in pairs is tibia, second femur
        d[i] = d_aux;
        nt[i] = nt_aux;
        
        /*std::cout << "fem_Points_ind" << fem_Points_ind << std::endl;
        std::cout << "facesFem_ind" << facesFem_ind[0] << " " << facesFem_ind[1] << " " << facesFem_ind[2] << std::endl;
        std::cout << "tib_Points_ind=" << tib_Points_ind << std::endl;
        std::cout << "d= " << d << std::endl;
        std::cout << "pairs_list" << pairs_list[i][0] - 1 << " " << pairs_list[i][1] - 1 << std::endl;*/
        if (i > 0) {
            if ((pairs_list[i][0] == pairs_list[i - 1][0])&&(i<pairs_list.size()-1)) {
                l = l + 1;
            }
            else {
                Vector_<Vec3> d_aux_list(l - 1);
                Vector_<Vec3> nt_aux_list(l - 1);
                Vector_<Vec3> fem_v_aux_list(l - 1);
                Vector_<Real> highligh_index(l-1); //the index of penetration value highlighted with a high value
                Real mindist;
                Vec3 nt_l;
                for (int j = 0; j < l - 1; j++) {
                    d_aux_list[j] = d[i - l + j + 1];
                    nt_aux_list[j] = nt[i - l + j + 1];
                    fem_v_aux_list[j] = fem_v_ind[i - l + j + 1];
                }
                CalculateMaximumPenetration(d_aux_list, nt_aux_list,mindist,nt_l,highligh_index);
                
                Real p = CalculatePressure(poisson, E, mindist, h);

                // Compute velocity corresponding to the pair in contact (based on the face that it is in contact)
                //highlight index of maximum penetration
                Real b1 = 1e-7; // value to find the maximum
                Real max_tohighlight = log(sum(exp(b1 * highligh_index))) / b1;
                v_rel_femur_wrtTibia.resizeKeep(k);
                for (int l2 = 0; l2 < highligh_index.size(); l2++) {
                    v_rel_femur_wrtTibia[k - 1] = v_rel_femur_wrtTibia[k - 1] + (highligh_index[l2] / max_tohighlight) * fem_v_aux_list[l2];
                }
                // prepare computation tangential force
                Real vt = 0.2;
                Real us = 0.8;
                Real uv = 0.5;
                Real ud = 0.8;
                Real vnormal = dot(v_rel_femur_wrtTibia[k - 1], nt_aux_list[0]);
                Vec3 vtangent = v_rel_femur_wrtTibia[k - 1] - vnormal * nt_aux_list[0];
                Real vslip = pow(vtangent.normSqr() + 1e-10, 1. / 2.);
                Real vrel = vslip / vt;

                force_s_l.resizeKeep(k);
                // compute normal force
                Vec3 normal_force = p * At * nt_l;
                
                // compute tangential force
                Real tangential_force_scalar= p*At* (std::min(vrel, Real(1)) *
                    (ud + 2 * (us - ud) / (1 + vrel * vrel)) + uv * vslip);
                Vec3 tangential_force = tangential_force_scalar * vtangent / vslip;
                // compute total contact force at that face
                force_s_l[k - 1] = normal_force+tangential_force;
                
                Ctib.resizeKeep(k);
                Ctib[k - 1] = Ctib_aux;
                face_s.resizeKeep(k);
                face_s[k - 1] = pairs_list[i - 2][0];
                mom_O.resizeKeep(k);
                mom_O[k - 1] = cross(Ctib_aux - originTib_G, force_s_l[k - 1]);



                k = k + 1;
                l = 2;
                while (k < pairs_list[i][0]) {
                    k = k + 1;
                }
                ///////////////

            }
        }
        else {
            l = l + 1;
        }

    }
    
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
    

}

void ComputeKneeContactForces(Model *model, State *state,Vec3 knee_trans, Vec3 knee_rot, Vec3 &SumForces, Vec3 &SumMoments, Real &SumForces_vert_Lat, Real &SumForces_vert_Med) {
    // Read Geometry Information
    std::string filename_tibPoints("C:/Gil/MeshesInAD/contactsKneeProsthesis/tibPoints.csv");
    Vector_<Vec3> tibPoints = ReadDataDoublex3columns(filename_tibPoints, 36);
    std::string filename_femPoints("C:/Gil/MeshesInAD/contactsKneeProsthesis/femPoints.csv");
    Vector_<Vec3> femPoints = ReadDataDoublex3columns(filename_femPoints, 117);
    //std::string filename_tibconList1("C:/Gil/MeshesInAD/contactsKneeProsthesis/tib1_connectivityList.csv");
    //std::vector<std::vector<int>> tibConList1 = ReadDataIntx3columns(filename_tibconList1, 26);
    //std::string filename_tibconList2("C:/Gil/MeshesInAD/contactsKneeProsthesis/tib2_connectivityList.csv");
    //std::vector<std::vector<int>> tibConList2 = ReadDataIntx3columns(filename_tibconList2, 23);
    std::string filename_pairs1("C:/Gil/MeshesInAD/contactsKneeProsthesis/pairs1.csv");
    std::vector<std::vector<int>> pairs1_list = ReadDataIntx2columns(filename_pairs1, 287);
    std::string filename_pairs2("C:/Gil/MeshesInAD/contactsKneeProsthesis/pairs2.csv");
    std::vector<std::vector<int>> pairs2_list = ReadDataIntx2columns(filename_pairs2, 212);
    std::string filename_facesFem("C:/Gil/MeshesInAD/contactsKneeProsthesis/facesFem.csv");
    std::vector<std::vector<int>> facesFem = ReadDataIntx3columns(filename_facesFem, 185);
    std::string filename_facesTib1("C:/Gil/MeshesInAD/contactsKneeProsthesis/facesTib1.csv");
    std::vector<std::vector<int>> facesTib1 = ReadDataIntx3columns(filename_facesTib1, 26);
    std::string filename_facesTib2("C:/Gil/MeshesInAD/contactsKneeProsthesis/facesTib2.csv");
    std::vector<std::vector<int>> facesTib2 = ReadDataIntx3columns(filename_facesTib2, 23);

    // Compute relative velocities in tibia reference
    Vector_<Vec3> v_femPoints_wrtTibia(femPoints.size(), Vec3(0));
    MobilizedBody femoral_component = model->getBodySet().get("femoral_component").getMobilizedBody();
    MobilizedBody tibial_tray = model->getBodySet().get("tibial_tray").getMobilizedBody();
    for (int i = 0; i < femPoints.size(); i++) {
        v_femPoints_wrtTibia[i]= femoral_component.findStationVelocityInAnotherBody(*state, femPoints[i], tibial_tray);
    }

    // Sum shift translation to the femur
    femPoints = femPoints + Vec3(0.0, 0.042, 0.0);

    //test with numerical values
    //knee_trans = Vec3(0.01, 0.03, 0.008);
    //knee_rot = Vec3(0.1, 0.03, 0.025);

    //// Apply transformations
    // Get matrix transformation of tibia with respect to the femur
    Mat44 Mtransf_tib = ftransf_function(knee_trans, knee_rot);
    std::cout << "Mtransf_tib=" << Mtransf_tib << std::endl;

    // Apply the transformation to all points of the tibia
    Vec4 tibPoints_transf_aux(1);
    Vector_<Vec3> tibPoints_transf(36);
    for (int i = 0; i < 36; i++) {
        tibPoints_transf_aux = Mtransf_tib*Vec4(tibPoints[i][0], tibPoints[i][1], tibPoints[i][2], 1.0);
        for (int j = 0; j < 3; j++) {
            tibPoints_transf[i][j] = tibPoints_transf_aux[j];
        }
    }

    // Calculate origin of the tibia...
    Vec4 originTib_G4 = Mtransf_tib*Vec4(0, 0, 0, 1);
    Vec3 originTib_G = Vec3(originTib_G4[0], originTib_G4[1], originTib_G4[2]);
    // 

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

    CalculateForceCompartment(femPoints, v_femPoints_wrtTibia, facesFem, tibPoints_transf, facesTib1, pairs1_list, SumForces1, SumMoments1, poisson, E, h, originTib_G, knee_trans, knee_rot);
    std::cout << "forces=" << SumForces1 << " moments=" << SumMoments1 << std::endl;
    CalculateForceCompartment(femPoints, v_femPoints_wrtTibia, facesFem, tibPoints_transf, facesTib2, pairs2_list, SumForces2, SumMoments2, poisson, E, h, originTib_G, knee_trans, knee_rot);
    SumForces = SumForces1 + SumForces2;
    SumMoments = SumMoments1 + SumMoments2;
    SumForces_vert_Lat = SumForces1[1];
    SumForces_vert_Med = SumForces2[1];
}
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

    // States and controls
    T ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    T up[NP]; /// contact model parameters - parameters
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states
    
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

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    Vec3 knee_trans = Vec3(model->getStateVariableValue(*state, "knee_r/knee_tx_r/value"),
        model->getStateVariableValue(*state, "knee_r/knee_ty_r/value"),
        model->getStateVariableValue(*state, "knee_r/knee_tz_r/value"));
    Vec3 knee_rot_0 = Vec3(model->getStateVariableValue(*state, "knee_r/knee_angle_r/value"),
        model->getStateVariableValue(*state, "knee_r/knee_adduction_r/value"),
        model->getStateVariableValue(*state, "knee_r/knee_rotation_r/value"));
    
    Vec3 knee_rot(0);
    knee_rot = Vec3(-knee_rot_0[0], -knee_rot_0[1], knee_rot_0[2]); // change sign to be consistent with joint definition (reverse in .osim)

    // Compute Resulting contact wrench at the center of the tibia
    Vec3 KneeCont_SumForces;
    Vec3 KneeCont_SumMoments;
    KneeCont_SumForces.setToZero();
    KneeCont_SumMoments.setToZero();
    Real SumForces_vert_Lat;
    Real SumForces_vert_Med;
    
    ComputeKneeContactForces(model,state,knee_trans, knee_rot, KneeCont_SumForces, KneeCont_SumMoments, SumForces_vert_Lat, SumForces_vert_Med);
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
    
    /// Knee contact forces
    res[0][ndof + nc + nc + nc + nc] = value<T>(SumForces_vert_Lat);
    res[0][ndof + nc + nc + nc + nc + 1] = value<T>(SumForces_vert_Med);
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

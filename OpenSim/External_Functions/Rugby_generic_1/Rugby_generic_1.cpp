#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/PolynomialFunction.h>
#include <OpenSim/Common/MultiplierFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
#include <OpenSim/Simulation/Model/BushingForce.h>
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

constexpr int n_in = 2;
constexpr int n_out = 1;
constexpr int nCoordinates = 91+6; //91 MSK model + 6 PB
constexpr int NX = nCoordinates * 2;
constexpr int NU = nCoordinates; // joint accelerations, one for each dof
constexpr int NR = nCoordinates +52; // joint moments/forces, one for each dof + 6 resultant GRF, 10*3 contact forces for each of the 10 spheres, 6 resultant GRM, 10 foot-ground contact deformation power

template<typename T>
T value(const Recorder& e) { return e; };
template<>
double value(const Recorder& e) { return e.getValue(); };

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

Vec3 find_point_projection_along_a_line(Vec3 P0, Vec3 P1, Vec3 P2) {
	Vec3 r = P2 - P1;
	Vec3 r_hat = r.normalize();
	Vec3 q = P0 - P1;
	Real t = dot(q, r);
	Real cos_theta = t / (q.norm() * r.norm());
	Real projection = q.norm() * cos_theta;
	Vec3 Q = projection * r_hat;
	Vec3 Q_inG = P1 + Q;

	return Q_inG;
}

void calc_ShoulderContact(State* state, OpenSim::Body* clavicle, Vec3 sphere_loc, OpenSim::Body* cylinder, Vec3& Contact_Forces) {
	Real cylinder_radius = 0.1702085;
	Real cylinder_height = 0.5386;
	Real sphere_radius = 0.025;

	Vec3 cylinder_bottom = Vec3(0, -cylinder_height/2, 0); //in local punching_bag frame //TO BE MODIFIED
	Vec3 cylinder_top = Vec3(0, cylinder_height/2, 0); // in local punching_bag frame // TO BE MODIFIED

	Vec3 sphere_loc_inG= clavicle->findStationLocationInGround(*state, sphere_loc);
	Vec3 cylinder_bottom_inG = cylinder->findStationLocationInGround(*state, cylinder_bottom);
	Vec3 cylinder_top_inG = cylinder->findStationLocationInGround(*state, cylinder_top);

	// projection of sphere_loc_inG on the longitudinal axis of the punching bag

	Vec3 Q_inG = find_point_projection_along_a_line(sphere_loc_inG, cylinder_bottom_inG, cylinder_top_inG);
	
	Vec3 d, cylinder_edge, sphere_edge, n;
	get_distance_between_edges(Q_inG, sphere_loc_inG, sphere_radius,cylinder_radius,)


}

void calc_BushingForce(State* state, const MobilizedBody* body1, const MobilizedBody* body2, Vec3 loc_body1, Vec3 loc_body2, Vec3 Kt, Vec3 Dt, Vec3 Kr, Vec3 Dr, Vector_<SpatialVec>& BF_Forces, Vec3& qt, Vec3& qdott, Vec3& qr, Vec3& qdotr) {
	// This function calculates the forces between subject and exoskeleton, 
	// assuming bushing forces between bodies

	osim_double_adouble dist = body1->calcStationToStationDistance(*state, loc_body1, *body2, loc_body2);
	Transform X_GB1 = body1->getBodyTransform(*state);
	Transform X_GB2 = body2->getBodyTransform(*state);

	Rotation NulRot;
	NulRot.setRotationToIdentityMatrix();
	Transform T_B1(NulRot, loc_body1);
	Transform T_B2(NulRot, loc_body2);

	Transform X_GF;
	Transform X_GM;
	Transform X_FM;
	X_GF = X_GB1 * T_B1;
	X_GM = X_GB2 * T_B2;
	X_FM = ~X_GF * X_GM;
	Vec3 p_B1F_G = X_GB1.R() * T_B1.p();   // 15 flops
	Vec3 p_B2M_G = X_GB2.R() * T_B2.p();   // 15 flops
	Vec3 p_FM_G = X_GF.R() * X_FM.p();    // 15 flops

	qr = X_FM.R().convertRotationToBodyFixedXYZ();
	qt = X_FM.p();
	//Vec3 ang2 = X_FM.R().convertRotationToBodyFixedXYZ();
	SpatialVec V_GB1 = body1->getBodyVelocity(*state);
	SpatialVec V_GB2 = body2->getBodyVelocity(*state);

	SpatialVec V_GF = SpatialVec(V_GB1[0], V_GB1[1] + V_GB1[0] % p_B1F_G);
	SpatialVec V_GM = SpatialVec(V_GB2[0], V_GB2[1] + V_GB2[0] % p_B2M_G);

	// This is the velocity of M in F, but with the time derivative
	// taken in the Ground frame.
	const SpatialVec V_FM_G = V_GM - V_GF;

	// To get derivative in F, we must remove the part due to the
	// angular velocity w_GF of F in G.
	SpatialVec V_FM = ~X_GF.R() * SpatialVec(V_FM_G[0],
		V_FM_G[1] - V_GF[0] % p_FM_G);

	// Need angular velocity in M frame for conversion to qdot.
	const Vec3  w_FM_M = ~X_FM.R() * V_FM[0];
	const Mat33 N_FM = Rotation::calcNForBodyXYZInBodyFrame(qr);
	qdotr = N_FM * w_FM_M;
	qdott = V_FM[1];

	Vec3 fkt;
	Vec3 fvt;
	Vec3 ft;
	Vec3 fM;
	Vec3 fMv;
	Vec3 fMk;
	for (int i = 0; i < 3; ++i) {
		fkt[i] = Kt[i] * qt[i];
		fvt[i] = Dt[i] * qdott[i];
		ft[i] = fkt[i] + fvt[i];
		fMk[i] = Kr[i] * qr[i];
		fMv[i] = Dr[i] * qdotr[i];
		fM[i] = fMk[i] + fMv[i];
	}

	Vec3 fB2_q = -fM; // in q basis
	Vec3 fM_F = -ft; // acts at M, but exp. in F frame
	//// Calculate the matrix relating q-space generalized forces to a real-space
	//// moment vector. We know qforce = ~H * moment (where H is the
	//// the hinge matrix for a mobilizer using qdots as generalized speeds).
	//// In that case H would be N^-1, qforce = ~(N^-1)*moment so
	//// moment = ~N*qforce. Caution: our N wants the moment in the outboard
	//// body frame, in this case M.
	//Mat33 N_FM = Rotation::calcNForBodyXYZInBodyFrame(ang2);
	Vec3  mB2_M = ~N_FM * fB2_q; // moment acting on body 2, exp. in M
	Vec3  mB2_G = X_GM.R() * mB2_M; // moment on body 2, now exp. in G
	// Transform force from F frame to ground. This is the force to 
	// apply to body 2 at point OM; -f goes on body 1 at the same
	// spatial location. Here we actually apply it at OF so we have to
	// account for the moment produced by the shift from OM.
	Vec3 fM_G = X_GF.R() * fM_F;
	SpatialVec F_GM = SpatialVec(mB2_G, fM_G);
	SpatialVec F_GF = SpatialVec(-(mB2_G + p_FM_G % fM_G), -fM_G);

	// Shift forces to body origins.
	SpatialVec F_GB2 = SpatialVec(F_GM[0] + p_B2M_G % F_GM[1], F_GM[1]);
	SpatialVec F_GB1 = SpatialVec(F_GF[0] + p_B1F_G % F_GF[1], F_GF[1]);

	BF_Forces[body1->getMobilizedBodyIndex()] = F_GB1;
	BF_Forces[body2->getMobilizedBodyIndex()] = F_GB2;

}

void expressBodyForcesInGround(State* state, OpenSim::Body* body1, OpenSim::Body* body2, OpenSim::PhysicalOffsetFrame* body1_offset, OpenSim::PhysicalOffsetFrame* body2_offset, Array<osim_double_adouble> bushing_values, Vector_<SpatialVec>& appliedBodyForces) {
	Rotation R_G_F1 = body1_offset->getTransformInGround(*state).R();
	Vec3 f1_G = R_G_F1 * Vec3(bushing_values[0], bushing_values[1], bushing_values[2]); // express force in ground
	Vec3 tau1_G = R_G_F1 * Vec3(bushing_values[3], bushing_values[4], bushing_values[5]); //express moment applied to frame origin in ground
	Vec3 pF1_G = body1_offset->getTransformInGround(*state).p(); // origin of frame 1 in ground
	Vec3 pB1_G = body1->getTransformInGround(*state).p();
	Vec3 r_B1_to_F1 = pF1_G - pB1_G;
	Vec3 tau1_aboutBodyOrigin_G = tau1_G + r_B1_to_F1 % f1_G;
	appliedBodyForces[body1->getMobilizedBodyIndex()] +=
		SimTK::SpatialVec(tau1_aboutBodyOrigin_G, f1_G);

	Rotation R_G_F2 = body2_offset->getTransformInGround(*state).R();
	Vec3 f2_G = R_G_F2 * Vec3(bushing_values[6], bushing_values[7], bushing_values[8]); // express force in ground
	Vec3 tau2_G = R_G_F2 * Vec3(bushing_values[9], bushing_values[10], bushing_values[11]); //express moment applied to frame origin in ground
	Vec3 pF2_G = body2_offset->getTransformInGround(*state).p(); // origin of frame 2 in ground
	Vec3 pB2_G = body2->getTransformInGround(*state).p();
	Vec3 r_B2_to_F2 = pF2_G - pB2_G;
	Vec3 tau2_aboutBodyOrigin_G = tau2_G + r_B2_to_F2 % f2_G;
	appliedBodyForces[body2->getMobilizedBodyIndex()] +=
		SimTK::SpatialVec(tau2_aboutBodyOrigin_G, f2_G);

}


template<typename T>
int F_generic(const T** arg, T** res) {

	/*OpenSim::Model* model = new Model("C:/Gil/Collaborations/DarioCazzola_AndreaBraschi/Model2dll/modifiedWrapping.osim");
	std::cout << model->getNumCoordinates() << std::endl;*/

	// Definition of model.
	OpenSim::Model* model;
	model = new OpenSim::Model();

	// Definition of bodies.
	OpenSim::Body* pelvis;
	pelvis = new OpenSim::Body("pelvis", 20.69340000000000046043, Vec3(-0.07206569999999999643, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.22324099999999999500, 0.22324099999999999500, 0.10131900000000000628, 0., 0., 0.));
	model->addBody(pelvis);

	OpenSim::Body* spine;
	spine = new OpenSim::Body("spine", 2.07363304981774021485, Vec3(0.00000000000000000000, 0.24810299999999999021, 0.00000000000000000000), Inertia(0.07332260000000000166, 0.02637329999999999883, 0.07332260000000000166, 0., 0., 0.));
	model->addBody(spine);

	OpenSim::Body* lclavicle;
	lclavicle = new OpenSim::Body("lclavicle", 0.79590000000000005187, Vec3(-0.03406150000000000150, 0.01896089999999999928, -0.11160799999999999887), Inertia(0.00027893400000000000, 0.00030217900000000002, 0.00004648900000000000, 0., 0., 0.));
	model->addBody(lclavicle);

	OpenSim::Body* lscapula;
	lscapula = new OpenSim::Body("lscapula", 3.63840000000000030056, Vec3(-0.05357540000000000219, -0.05914670000000000344, 0.02914749999999999980), Inertia(0.00104933999999999995, 0.00097317599999999997, 0.00115934999999999996, 0., 0., 0.));
	model->addBody(lscapula);

	OpenSim::Body* humerus_l;
	humerus_l = new OpenSim::Body("humerus_l", 3.86580000000000056914, Vec3(0.00000000000000000000, -0.16607099999999999640, 0.00000000000000000000), Inertia(0.02214249999999999899, 0.00763849000000000038, 0.02485429999999999931, 0., 0., 0.));
	model->addBody(humerus_l);

	OpenSim::Body* torso;
	torso = new OpenSim::Body("torso", 23.28146695018225997842, Vec3(0.00000000000000000000, 0.13396299999999999875, 0.00000000000000000000), Inertia(1.50733000000000005869, 1.34272999999999997911, 1.50733000000000005869, 0., 0., 0.));
	model->addBody(torso);

	OpenSim::Body* cerv7;
	cerv7 = new OpenSim::Body("cerv7", 0.90960000000000007514, Vec3(-0.00019934500000000000, 0.00707674999999999971, 0.00000000000000000000), Inertia(0.02760270000000000085, 0.11041099999999999526, 0.02760270000000000085, 0., 0., 0.));
	model->addBody(cerv7);

	OpenSim::Body* cerv6;
	cerv6 = new OpenSim::Body("cerv6", 0.79590000000000005187, Vec3(0.00019934500000000000, 0.00637903999999999988, 0.00000000000000000000), Inertia(0.02822449999999999959, 0.11289799999999999836, 0.02822449999999999959, 0., 0., 0.));
	model->addBody(cerv6);

	OpenSim::Body* cerv5;
	cerv5 = new OpenSim::Body("cerv5", 0.68220000000000002860, Vec3(0.00009967250000000000, 0.00627937000000000026, 0.00000000000000000000), Inertia(0.02922299999999999898, 0.11689199999999999591, 0.02922299999999999898, 0., 0., 0.));
	model->addBody(cerv5);

	OpenSim::Body* cerv4;
	cerv4 = new OpenSim::Body("cerv4", 0.68220000000000002860, Vec3(-0.00029901800000000002, 0.00667806000000000042, 0.00000000000000000000), Inertia(0.02962220000000000131, 0.11848899999999999710, 0.02962220000000000131, 0., 0., 0.));
	model->addBody(cerv4);

	OpenSim::Body* cerv3;
	cerv3 = new OpenSim::Body("cerv3", 0.68220000000000002860, Vec3(-0.00079737999999999999, 0.00707674999999999971, 0.00000000000000000000), Inertia(0.02986700000000000119, 0.11946800000000000475, 0.02986700000000000119, 0., 0., 0.));
	model->addBody(cerv3);

	OpenSim::Body* cerv2;
	cerv2 = new OpenSim::Body("cerv2", 0.90960000000000007514, Vec3(-0.00099672499999999995, 0.00588068000000000011, 0.00000000000000000000), Inertia(0.02743970000000000090, 0.10975899999999999546, 0.02743970000000000090, 0., 0., 0.));
	model->addBody(cerv2);

	OpenSim::Body* cerv1;
	cerv1 = new OpenSim::Body("cerv1", 0.68220000000000002860, Vec3(0.03109779999999999836, 0.00189378000000000007, 0.00000000000000000000), Inertia(0.03066970000000000113, 0.12267899999999999638, 0.03066970000000000113, 0., 0., 0.));
	model->addBody(cerv1);

	OpenSim::Body* skull;
	skull = new OpenSim::Body("skull", 4.66170000000000062101, Vec3(0.00000000000000000000, 0.07558800000000000241, 0.00000000000000000000), Inertia(0.03664930000000000271, 0.01235110000000000027, 0.03664930000000000271, 0., 0., 0.));
	model->addBody(skull);

	OpenSim::Body* jaw;
	jaw = new OpenSim::Body("jaw", 0.68220000000000002860, Vec3(0.00000000000000000000, 0.05590830000000000105, 0.00000000000000000000), Inertia(0.04369809999999999667, 0.01976179999999999951, 0.04369809999999999667, 0., 0., 0.));
	model->addBody(jaw);

	OpenSim::Body* rclavicle;
	rclavicle = new OpenSim::Body("rclavicle", 0.79590000000000005187, Vec3(-0.03394029999999999964, 0.01889340000000000117, 0.11121100000000000430), Inertia(0.00027695200000000001, 0.00030003099999999997, 0.00004615870000000000, 0., 0., 0.));
	model->addBody(rclavicle);

	OpenSim::Body* rscapula;
	rscapula = new OpenSim::Body("rscapula", 3.63840000000000030056, Vec3(-0.04891000000000000208, -0.05399620000000000114, -0.02660929999999999893), Inertia(0.00087454199999999996, 0.00081106699999999995, 0.00096622799999999997, 0., 0., 0.));
	model->addBody(rscapula);

	OpenSim::Body* humerus_r;
	humerus_r = new OpenSim::Body("humerus_r", 3.86580000000000056914, Vec3(0.00000000000000000000, -0.16137899999999999467, 0.00000000000000000000), Inertia(0.02090880000000000169, 0.00721288999999999958, 0.02346950000000000078, 0., 0., 0.));
	model->addBody(humerus_r);

	OpenSim::Body* ulna_l;
	ulna_l = new OpenSim::Body("ulna_l", 1.13700000000000001066, Vec3(0.00000000000000000000, -0.11629200000000000648, 0.00000000000000000000), Inertia(0.00502908999999999995, 0.00104928000000000009, 0.00545525999999999997, 0., 0., 0.));
	model->addBody(ulna_l);

	OpenSim::Body* radius_l;
	radius_l = new OpenSim::Body("radius_l", 1.13700000000000001066, Vec3(0.00000000000000000000, -0.11629200000000000648, 0.00000000000000000000), Inertia(0.00502908999999999995, 0.00104928000000000009, 0.00545525999999999997, 0., 0., 0.));
	model->addBody(radius_l);

	OpenSim::Body* hand_l;
	hand_l = new OpenSim::Body("hand_l", 0.28424999999999317479, Vec3(0.00000000000000000000, -0.07884059999999999679, 0.00000000000000000000), Inertia(0.00176461999999999995, 0.00108212000000000000, 0.00265088999999999993, 0., 0., 0.));
	model->addBody(hand_l);

	OpenSim::Body* ulna_r;
	ulna_r = new OpenSim::Body("ulna_r", 1.13700000000000001066, Vec3(0.00000000000000000000, -0.11966599999999999460, 0.00000000000000000000), Inertia(0.00532513999999999991, 0.00111104999999999993, 0.00577640000000000011, 0., 0., 0.));
	model->addBody(ulna_r);

	OpenSim::Body* radius_r;
	radius_r = new OpenSim::Body("radius_r", 1.13700000000000001066, Vec3(0.00000000000000000000, -0.11966599999999999460, 0.00000000000000000000), Inertia(0.00532513999999999991, 0.00111104999999999993, 0.00577640000000000011, 0., 0., 0.));
	model->addBody(radius_r);

	OpenSim::Body* hand_r;
	hand_r = new OpenSim::Body("hand_r", 0.28424999999999317479, Vec3(0.00000000000000000000, -0.07455109999999999515, 0.00000000000000000000), Inertia(0.00157782999999999990, 0.00096757199999999999, 0.00237029000000000016, 0., 0., 0.));
	model->addBody(hand_r);

	OpenSim::Body* femur_r;
	femur_r = new OpenSim::Body("femur_r", 11.25630000000000130456, Vec3(0.00000000000000000000, -0.17000000000000001221, 0.00000000000000000000), Inertia(0.17274999999999998690, 0.04528000000000000080, 0.18216999999999999860, 0., 0., 0.));
	model->addBody(femur_r);

	OpenSim::Body* tibia_r;
	tibia_r = new OpenSim::Body("tibia_r", 4.54800000000000004263, Vec3(0.00000000000000000000, -0.18670000000000000484, 0.00000000000000000000), Inertia(0.04649999999999999967, 0.00820000000000000069, 0.04649999999999999967, 0., 0., 0.));
	model->addBody(tibia_r);

	OpenSim::Body* talus_r;
	talus_r = new OpenSim::Body("talus_r", 0.11370000000000000939, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00100000000000000002, 0.00100000000000000002, 0.00100000000000000002, 0., 0., 0.));
	model->addBody(talus_r);

	OpenSim::Body* calcn_r;
	calcn_r = new OpenSim::Body("calcn_r", 1.47809999999999996945, Vec3(0.10000000000000000555, 0.02999999999999999889, 0.00000000000000000000), Inertia(0.00178999999999999992, 0.00498999999999999964, 0.00525000000000000033, 0., 0., 0.));
	model->addBody(calcn_r);

	OpenSim::Body* toes_r;
	toes_r = new OpenSim::Body("toes_r", 0.22740000000000001878, Vec3(0.03459999999999999881, 0.00600000000000000012, -0.01750000000000000167), Inertia(0.00013999999999999999, 0.00027999999999999998, 0.00138999999999999996, 0., 0., 0.));
	model->addBody(toes_r);

	OpenSim::Body* femur_l;
	femur_l = new OpenSim::Body("femur_l", 11.25630000000000130456, Vec3(0.00000000000000000000, -0.17000000000000001221, 0.00000000000000000000), Inertia(0.17274999999999998690, 0.04528000000000000080, 0.18216999999999999860, 0., 0., 0.));
	model->addBody(femur_l);

	OpenSim::Body* tibia_l;
	tibia_l = new OpenSim::Body("tibia_l", 4.54800000000000004263, Vec3(0.00000000000000000000, -0.18670000000000000484, 0.00000000000000000000), Inertia(0.04649999999999999967, 0.00820000000000000069, 0.04649999999999999967, 0., 0., 0.));
	model->addBody(tibia_l);

	OpenSim::Body* talus_l;
	talus_l = new OpenSim::Body("talus_l", 0.11370000000000000939, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00100000000000000002, 0.00100000000000000002, 0.00100000000000000002, 0., 0., 0.));
	model->addBody(talus_l);

	OpenSim::Body* calcn_l;
	calcn_l = new OpenSim::Body("calcn_l", 1.47809999999999996945, Vec3(0.10000000000000000555, 0.02999999999999999889, 0.00000000000000000000), Inertia(0.00178999999999999992, 0.00498999999999999964, 0.00525000000000000033, 0., 0., 0.));
	model->addBody(calcn_l);

	OpenSim::Body* toes_l;
	toes_l = new OpenSim::Body("toes_l", 0.22740000000000001878, Vec3(0.03459999999999999881, 0.00600000000000000012, -0.01750000000000000167), Inertia(0.00013999999999999999, 0.00027999999999999998, 0.00138999999999999996, 0., 0., 0.));
	model->addBody(toes_l);

	OpenSim::Body* punching_bag;
	punching_bag = new OpenSim::Body("punching_bag", 50, Vec3(0), Inertia(5.1942899999999996, 0.72250000000000003, 5.1942899999999996, 0, 0, 0));
	model->addBody(punching_bag);

	// Definition of joints.
	SpatialTransform st_pelvisjnt;
	st_pelvisjnt[0].setCoordinateNames(OpenSim::Array<std::string>("gndpitch", 1, 1));
	st_pelvisjnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_pelvisjnt[1].setCoordinateNames(OpenSim::Array<std::string>("gndroll", 1, 1));
	st_pelvisjnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_pelvisjnt[2].setCoordinateNames(OpenSim::Array<std::string>("gndyaw", 1, 1));
	st_pelvisjnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_pelvisjnt[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
	st_pelvisjnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_pelvisjnt[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
	st_pelvisjnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_pelvisjnt[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
	st_pelvisjnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pelvisjnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* pelvisjnt;
	pelvisjnt = new OpenSim::CustomJoint("pelvisjnt", model->getGround(), Vec3(0.08000000000000000167, -0.17299999999999998712, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *pelvis, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_pelvisjnt);

	SpatialTransform st_spine_pelvis;
	st_spine_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("spine_tilt", 1, 1));
	st_spine_pelvis[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_spine_pelvis[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_spine_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("spine_list", 1, 1));
	st_spine_pelvis[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_spine_pelvis[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_spine_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("spine_rotation", 1, 1));
	st_spine_pelvis[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_spine_pelvis[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_spine_pelvis[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.01931653053344994930));
	st_spine_pelvis[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_spine_pelvis[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.01931653053344994930));
	st_spine_pelvis[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_spine_pelvis[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.17794927931372006569));
	st_spine_pelvis[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* spine_pelvis;
	spine_pelvis = new OpenSim::CustomJoint("spine_pelvis", *pelvis, Vec3(-0.08154530000000000112, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *spine, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_spine_pelvis);

	SpatialTransform st_auxSTERNljnt;
	st_auxSTERNljnt[0].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNljnt_r1", 1, 1));
	st_auxSTERNljnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNljnt[0].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNljnt[1].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNljnt_r2", 1, 1));
	st_auxSTERNljnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNljnt[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNljnt[2].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNljnt_r3", 1, 1));
	st_auxSTERNljnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNljnt[2].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxSTERNljnt[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.33748471554110004433));
	st_auxSTERNljnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNljnt[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.11636232671806001626));
	st_auxSTERNljnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNljnt[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.15678104713489005029));
	st_auxSTERNljnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* auxSTERNljnt;
	auxSTERNljnt = new OpenSim::CustomJoint("auxSTERNljnt", *spine, Vec3(-0.01578229999999999911, 0.56744700000000003470, -0.02544920000000000168), Vec3(-0.02000000000000000042, -0.13000000000000000444, -0.10000000000000000555), *lclavicle, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.20000000000000001110, 0.40000000000000002220, 0.00000000000000000000), st_auxSTERNljnt);

	SpatialTransform st_auxSCAPljnt;
	st_auxSCAPljnt[0].setCoordinateNames(OpenSim::Array<std::string>("rollLSCA", 1, 1));
	st_auxSCAPljnt[0].setFunction(new LinearFunction(1.0000, 0.0106));
	st_auxSCAPljnt[0].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPljnt[1].setCoordinateNames(OpenSim::Array<std::string>("pitchLSCA", 1, 1));
	st_auxSCAPljnt[1].setFunction(new LinearFunction(-1.0000, 0.7000));
	st_auxSCAPljnt[1].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPljnt[2].setCoordinateNames(OpenSim::Array<std::string>("yawLSCA", 1, 1));
	st_auxSCAPljnt[2].setFunction(new LinearFunction(1.0000, -0.0389));
	st_auxSCAPljnt[2].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxSCAPljnt[3].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPljnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPljnt[4].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPljnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPljnt[5].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPljnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* auxSCAPljnt;
	auxSCAPljnt = new OpenSim::CustomJoint("auxSCAPljnt", *lclavicle, Vec3(-0.07380000000000000449, -0.02270770000000000080, -0.15895400000000001195), Vec3(-0.78000000000000002665, 0.29999999999999998890, 0.14000000000000001332), *lscapula, Vec3(-0.02906469999999999887, -0.05328530000000000072, 0.00775058999999999991), Vec3(0.45000000000000001110, -1.00000000000000000000, 0.10000000000000000555), st_auxSCAPljnt);

	SpatialTransform st_acromial_l;
	st_acromial_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
	st_acromial_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_acromial_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
	st_acromial_l[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
	st_acromial_l[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_l[2].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.84510758548058095840));
	st_acromial_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.01398270580809990449));
	st_acromial_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.39409739709527991103));
	st_acromial_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* acromial_l;
	acromial_l = new OpenSim::CustomJoint("acromial_l", *lscapula, Vec3(-0.00871941000000000035, -0.03875300000000000272, -0.01937650000000000136), Vec3(0.29999999999999998890, -0.29999999999999998890, 0.29999999999999998890), *humerus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_l);

	OpenSim::WeldJoint* ribcagejnt;
	ribcagejnt = new OpenSim::WeldJoint("ribcagejnt", *spine, Vec3(0.00000000000000000000, 0.19313099999999999712, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *torso, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	SpatialTransform st_auxt1jnt;
	st_auxt1jnt[0].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_r3", 1, 1));
	st_auxt1jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxt1jnt[1].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_r1", 1, 1));
	st_auxt1jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxt1jnt[2].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_r2", 1, 1));
	st_auxt1jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxt1jnt[3].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_t3", 1, 1));
	st_auxt1jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxt1jnt[4].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_t1", 1, 1));
	st_auxt1jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxt1jnt[5].setCoordinateNames(OpenSim::Array<std::string>("auxt1jnt_t2", 1, 1));
	st_auxt1jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxt1jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* auxt1jnt;
	auxt1jnt = new OpenSim::CustomJoint("auxt1jnt", *spine, Vec3(-0.06353050000000000364, 0.59647200000000000220, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv7, Vec3(0.00681760000000000005, -0.00563648000000000014, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_auxt1jnt);

	SpatialTransform st_aux7jnt;
	st_aux7jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_r3", 1, 1));
	st_aux7jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux7jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_r1", 1, 1));
	st_aux7jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux7jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_r2", 1, 1));
	st_aux7jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux7jnt[3].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_t3", 1, 1));
	st_aux7jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux7jnt[4].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_t1", 1, 1));
	st_aux7jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux7jnt[5].setCoordinateNames(OpenSim::Array<std::string>("aux7jnt_t2", 1, 1));
	st_aux7jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux7jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* aux7jnt;
	aux7jnt = new OpenSim::CustomJoint("aux7jnt", *cerv7, Vec3(0.01321259999999999969, 0.01233949999999999977, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv6, Vec3(0.00938814999999999962, -0.00405169000000000015, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux7jnt);

	SpatialTransform st_aux6jnt;
	st_aux6jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_r3", 1, 1));
	st_aux6jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux6jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_r1", 1, 1));
	st_aux6jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux6jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_r2", 1, 1));
	st_aux6jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux6jnt[3].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_t3", 1, 1));
	st_aux6jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux6jnt[4].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_t1", 1, 1));
	st_aux6jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux6jnt[5].setCoordinateNames(OpenSim::Array<std::string>("aux6jnt_t2", 1, 1));
	st_aux6jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux6jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* aux6jnt;
	aux6jnt = new OpenSim::CustomJoint("aux6jnt", *cerv6, Vec3(0.00981774000000000015, 0.00971807000000000053, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv5, Vec3(0.00565142999999999982, -0.00812330999999999991, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux6jnt);

	SpatialTransform st_aux5jnt;
	st_aux5jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_r3", 1, 1));
	st_aux5jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux5jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_r1", 1, 1));
	st_aux5jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux5jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_r2", 1, 1));
	st_aux5jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux5jnt[3].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_t3", 1, 1));
	st_aux5jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux5jnt[4].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_t1", 1, 1));
	st_aux5jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux5jnt[5].setCoordinateNames(OpenSim::Array<std::string>("aux5jnt_t2", 1, 1));
	st_aux5jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux5jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* aux5jnt;
	aux5jnt = new OpenSim::CustomJoint("aux5jnt", *cerv5, Vec3(0.00840238999999999917, 0.00742560000000000008, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv4, Vec3(0.00500355999999999978, -0.01168160000000000032, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux5jnt);

	SpatialTransform st_aux4jnt;
	st_aux4jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_r3", 1, 1));
	st_aux4jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux4jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_r1", 1, 1));
	st_aux4jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux4jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_r2", 1, 1));
	st_aux4jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux4jnt[3].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_t3", 1, 1));
	st_aux4jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux4jnt[4].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_t1", 1, 1));
	st_aux4jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux4jnt[5].setCoordinateNames(OpenSim::Array<std::string>("aux4jnt_t2", 1, 1));
	st_aux4jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux4jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* aux4jnt;
	aux4jnt = new OpenSim::CustomJoint("aux4jnt", *cerv4, Vec3(0.00729603000000000022, 0.00635911000000000028, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv3, Vec3(0.00469457000000000034, -0.00899046000000000047, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux4jnt);

	SpatialTransform st_aux3jnt;
	st_aux3jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_r3", 1, 1));
	st_aux3jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux3jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_r1", 1, 1));
	st_aux3jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux3jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_r2", 1, 1));
	st_aux3jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux3jnt[3].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_t3", 1, 1));
	st_aux3jnt[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[3].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux3jnt[4].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_t1", 1, 1));
	st_aux3jnt[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[4].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux3jnt[5].setCoordinateNames(OpenSim::Array<std::string>("aux3jnt_t2", 1, 1));
	st_aux3jnt[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux3jnt[5].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	OpenSim::CustomJoint* aux3jnt;
	aux3jnt = new OpenSim::CustomJoint("aux3jnt", *cerv3, Vec3(0.00523281000000000007, 0.00504342999999999979, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv2, Vec3(0.00541221999999999970, -0.01119320000000000037, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux3jnt);

	SpatialTransform st_aux2jnt;
	st_aux2jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux2jnt_r3", 1, 1));
	st_aux2jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux2jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux2jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux2jnt_r1", 1, 1));
	st_aux2jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux2jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux2jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux2jnt_r2", 1, 1));
	st_aux2jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux2jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux2jnt[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux2jnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux2jnt[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux2jnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux2jnt[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux2jnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* aux2jnt;
	aux2jnt = new OpenSim::CustomJoint("aux2jnt", *cerv2, Vec3(0.00686743999999999983, 0.02121030000000000137, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *cerv1, Vec3(0.04231100000000000139, 0.00350847000000000008, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux2jnt);

	SpatialTransform st_aux1jnt;
	st_aux1jnt[0].setCoordinateNames(OpenSim::Array<std::string>("aux1jnt_r3", 1, 1));
	st_aux1jnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux1jnt[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_aux1jnt[1].setCoordinateNames(OpenSim::Array<std::string>("aux1jnt_r1", 1, 1));
	st_aux1jnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux1jnt[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux1jnt[2].setCoordinateNames(OpenSim::Array<std::string>("aux1jnt_r2", 1, 1));
	st_aux1jnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_aux1jnt[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux1jnt[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux1jnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_aux1jnt[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux1jnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_aux1jnt[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99672504620965995947));
	st_aux1jnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* aux1jnt;
	aux1jnt = new OpenSim::CustomJoint("aux1jnt", *cerv1, Vec3(0.04320800000000000335, 0.01487109999999999994, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *skull, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_aux1jnt);

	OpenSim::WeldJoint* jawjnt;
	jawjnt = new OpenSim::WeldJoint("jawjnt", *skull, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *jaw, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	SpatialTransform st_auxSTERNrjnt;
	st_auxSTERNrjnt[0].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNrjnt_r1", 1, 1));
	st_auxSTERNrjnt[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNrjnt[0].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNrjnt[1].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNrjnt_r2", 1, 1));
	st_auxSTERNrjnt[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNrjnt[1].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNrjnt[2].setCoordinateNames(OpenSim::Array<std::string>("auxSTERNrjnt_r3", 1, 1));
	st_auxSTERNrjnt[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_auxSTERNrjnt[2].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxSTERNrjnt[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.33748471554110004433));
	st_auxSTERNrjnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNrjnt[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.11636232671806001626));
	st_auxSTERNrjnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxSTERNrjnt[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.15678104713489005029));
	st_auxSTERNrjnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* auxSTERNrjnt;
	auxSTERNrjnt = new OpenSim::CustomJoint("auxSTERNrjnt", *spine, Vec3(-0.01578229999999999911, 0.56744700000000003470, 0.02544920000000000168), Vec3(0.02000000000000000042, 0.13000000000000000444, -0.10000000000000000555), *rclavicle, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(-0.20000000000000001110, -0.40000000000000002220, 0.00000000000000000000), st_auxSTERNrjnt);

	SpatialTransform st_auxSCAPrjnt;
	st_auxSCAPrjnt[0].setCoordinateNames(OpenSim::Array<std::string>("rollRSCA", 1, 1));
	st_auxSCAPrjnt[0].setFunction(new LinearFunction(-1.0000, -0.0106));
	st_auxSCAPrjnt[0].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPrjnt[1].setCoordinateNames(OpenSim::Array<std::string>("pitchRSCA", 1, 1));
	st_auxSCAPrjnt[1].setFunction(new LinearFunction(1.0000, -0.7000));
	st_auxSCAPrjnt[1].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPrjnt[2].setCoordinateNames(OpenSim::Array<std::string>("yawRSCA", 1, 1));
	st_auxSCAPrjnt[2].setFunction(new LinearFunction(1.0000, -0.0389));
	st_auxSCAPrjnt[2].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_auxSCAPrjnt[3].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPrjnt[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPrjnt[4].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPrjnt[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_auxSCAPrjnt[5].setFunction(new Constant(0.00000000000000000000));
	st_auxSCAPrjnt[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* auxSCAPrjnt;
	auxSCAPrjnt = new OpenSim::CustomJoint("auxSCAPrjnt", *rclavicle, Vec3(-0.07353729999999999989, -0.02262689999999999840, 0.15838800000000000101), Vec3(0.78000000000000002665, -0.29999999999999998890, 0.14000000000000001332), *rscapula, Vec3(-0.02653379999999999975, -0.04864519999999999955, -0.00707567000000000005), Vec3(-0.45000000000000001110, 1.00000000000000000000, 0.10000000000000000555), st_auxSCAPrjnt);

	SpatialTransform st_acromial_r;
	st_acromial_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
	st_acromial_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_acromial_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
	st_acromial_r[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
	st_acromial_r[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_acromial_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.77151548957347304558));
	st_acromial_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.92568493897224202183));
	st_acromial_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.27269918565626993789));
	st_acromial_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* acromial_r;
	acromial_r = new OpenSim::CustomJoint("acromial_r", *rscapula, Vec3(-0.00796012999999999929, -0.03537830000000000141, 0.01768919999999999867), Vec3(-0.29999999999999998890, 0.29999999999999998890, 0.29999999999999998890), *humerus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_r);

	SpatialTransform st_elbow_l;
	st_elbow_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
	st_elbow_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_elbow_l[0].setAxis(Vec3(-0.22604695999999999123, -0.02226900000000000060, 0.97386183000000003940));
	st_elbow_l[1].setFunction(new Constant(0.00000000000000000000));
	st_elbow_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[2].setFunction(new Constant(0.00000000000000000000));
	st_elbow_l[2].setAxis(Vec3(0.97386183000000003940, -0.00000000000000000000, 0.22604695999999999123));
	st_elbow_l[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.00954054932107006870));
	st_elbow_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.09666873211731008553));
	st_elbow_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.00954054932107006870));
	st_elbow_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* elbow_l;
	elbow_l = new OpenSim::CustomJoint("elbow_l", *humerus_l, Vec3(0.01326940000000000064, -0.31394699999999997608, 0.00968654000000000043), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_l);

	SpatialTransform st_radioulnar_l;
	st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
	st_radioulnar_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radioulnar_l[0].setAxis(Vec3(-0.05639803000000000177, -0.99840645999999999560, 0.00195199999999999996));
	st_radioulnar_l[1].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[2].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_l[2].setAxis(Vec3(0.00195199999999999996, -0.00000000000000000000, 0.05639803000000000177));
	st_radioulnar_l[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.97641107638446500427));
	st_radioulnar_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.96486984786590102026));
	st_radioulnar_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.97641107638446500427));
	st_radioulnar_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radioulnar_l;
	radioulnar_l = new OpenSim::CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00656831999999999989, -0.01255009999999999980, -0.02546769999999999937), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_l);

	SpatialTransform st_radius_hand_l;
	st_radius_hand_l[0].setCoordinateNames(OpenSim::Array<std::string>("wrist_flex_l", 1, 1));
	st_radius_hand_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radius_hand_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_radius_hand_l[1].setCoordinateNames(OpenSim::Array<std::string>("wrist_dev_l", 1, 1));
	st_radius_hand_l[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radius_hand_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_l[2].setFunction(new Constant(0.00000000000000000000));
	st_radius_hand_l[2].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_l[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.97641107638446500427));
	st_radius_hand_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_l[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.96486984786590102026));
	st_radius_hand_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_l[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.97641107638446500427));
	st_radius_hand_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radius_hand_l;
	radius_hand_l = new OpenSim::CustomJoint("radius_hand_l", *radius_l, Vec3(-0.00858948999999999997, -0.22755600000000000827, -0.01328900000000000046), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radius_hand_l);

	SpatialTransform st_elbow_r;
	st_elbow_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
	st_elbow_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_elbow_r[0].setAxis(Vec3(0.22604695999999999123, 0.02226900000000000060, 0.97386183000000003940));
	st_elbow_r[1].setFunction(new Constant(0.00000000000000000000));
	st_elbow_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[2].setFunction(new Constant(0.00000000000000000000));
	st_elbow_r[2].setAxis(Vec3(0.97386183000000003940, 0.00000000000000000000, -0.22604695999999999123));
	st_elbow_r[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.98101267768981403883));
	st_elbow_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.06567876858103005588));
	st_elbow_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.98101267768981403883));
	st_elbow_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* elbow_r;
	elbow_r = new OpenSim::CustomJoint("elbow_r", *humerus_r, Vec3(0.01289440000000000031, -0.30507499999999998508, -0.00941282000000000056), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_r);

	SpatialTransform st_radioulnar_r;
	st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
	st_radioulnar_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radioulnar_r[0].setAxis(Vec3(0.05639803000000000177, 0.99840645999999999560, 0.00195199999999999996));
	st_radioulnar_r[1].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[2].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_r[2].setAxis(Vec3(0.00195199999999999996, 0.00000000000000000000, -0.05639803000000000177));
	st_radioulnar_r[3].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[4].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[5].setFunction(new Constant(0.00000000000000000000));
	st_radioulnar_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radioulnar_r;
	radioulnar_r = new OpenSim::CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00675888000000000010, -0.01291420000000000068, 0.02620660000000000003), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_r);

	SpatialTransform st_radius_hand_r;
	st_radius_hand_r[0].setCoordinateNames(OpenSim::Array<std::string>("wrist_flex_r", 1, 1));
	st_radius_hand_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radius_hand_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_radius_hand_r[1].setCoordinateNames(OpenSim::Array<std::string>("wrist_dev_r", 1, 1));
	st_radius_hand_r[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_radius_hand_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_r[2].setFunction(new Constant(0.00000000000000000000));
	st_radius_hand_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_r[3].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.00473974325811998831));
	st_radius_hand_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_r[4].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 0.99286366845818996296));
	st_radius_hand_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radius_hand_r[5].setFunction(new MultiplierFunction(new Constant(0.00000000000000000000), 1.00473974325811998831));
	st_radius_hand_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radius_hand_r;
	radius_hand_r = new OpenSim::CustomJoint("radius_hand_r", *radius_r, Vec3(-0.00883869999999999968, -0.23415800000000000503, 0.01367449999999999916), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radius_hand_r);

	SpatialTransform st_hip_r;
	st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
	st_hip_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
	st_hip_r[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
	st_hip_r[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_r[3].setFunction(new Constant(0.00000000000000000000));
	st_hip_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_r[4].setFunction(new Constant(0.00000000000000000000));
	st_hip_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_r[5].setFunction(new Constant(0.00000000000000000000));
	st_hip_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* hip_r;
	hip_r = new OpenSim::CustomJoint("hip_r", *pelvis, Vec3(-0.07069999999999999896, -0.06610000000000000597, 0.08350000000000000477), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_r);

	SpatialTransform st_knee_r;
	st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_knee_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_r[1].setFunction(new Constant(0.00000000000000000000));
	st_knee_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[2].setFunction(new Constant(0.00000000000000000000));
	st_knee_r[2].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_knee_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_knee_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[5].setFunction(new Constant(0.00000000000000000000));
	st_knee_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_r;
	knee_r = new OpenSim::CustomJoint("knee_r", *femur_r, Vec3(-0.00501372522863820051, -0.39658566558528168811, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_r);

	SpatialTransform st_ankle_r;
	st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
	st_ankle_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_ankle_r[0].setAxis(Vec3(-0.10501354999999999718, -0.17402244999999999520, 0.97912631999999999444));
	st_ankle_r[1].setFunction(new Constant(0.00000000000000000000));
	st_ankle_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[2].setFunction(new Constant(0.00000000000000000000));
	st_ankle_r[2].setAxis(Vec3(0.97912631999999999444, -0.00000000000000000000, 0.10501354999999999718));
	st_ankle_r[3].setFunction(new Constant(0.00000000000000000000));
	st_ankle_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[4].setFunction(new Constant(0.00000000000000000000));
	st_ankle_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[5].setFunction(new Constant(0.00000000000000000000));
	st_ankle_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* ankle_r;
	ankle_r = new OpenSim::CustomJoint("ankle_r", *tibia_r, Vec3(0.00000000000000000000, -0.42999999999999999334, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_r);

	SpatialTransform st_subtalar_r;
	st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
	st_subtalar_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_subtalar_r[0].setAxis(Vec3(0.78717961000000002958, 0.60474746000000001445, -0.12094949000000000672));
	st_subtalar_r[1].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[2].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_r[2].setAxis(Vec3(-0.12094949000000000672, 0.00000000000000000000, -0.78717961000000002958));
	st_subtalar_r[3].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[4].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[5].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* subtalar_r;
	subtalar_r = new OpenSim::CustomJoint("subtalar_r", *talus_r, Vec3(-0.04877000000000000085, -0.04195000000000000118, 0.00791999999999999996), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_r);

	SpatialTransform st_mtp_r;
	st_mtp_r[0].setCoordinateNames(OpenSim::Array<std::string>("mtp_angle_r", 1, 1));
	st_mtp_r[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_mtp_r[0].setAxis(Vec3(-0.58095439999999998193, 0.00000000000000000000, 0.81393610999999999045));
	st_mtp_r[1].setFunction(new Constant(0.00000000000000000000));
	st_mtp_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_mtp_r[2].setFunction(new Constant(0.00000000000000000000));
	st_mtp_r[2].setAxis(Vec3(0.81393610999999999045, -0.00000000000000000000, 0.58095439999999998193));
	st_mtp_r[3].setFunction(new Constant(0.00000000000000000000));
	st_mtp_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_mtp_r[4].setFunction(new Constant(0.00000000000000000000));
	st_mtp_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_mtp_r[5].setFunction(new Constant(0.00000000000000000000));
	st_mtp_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* mtp_r;
	mtp_r = new OpenSim::CustomJoint("mtp_r", *calcn_r, Vec3(0.17879999999999998672, -0.00200000000000000004, 0.00108000000000000001), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_mtp_r);

	SpatialTransform st_hip_l;
	st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
	st_hip_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
	st_hip_l[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
	st_hip_l[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_hip_l[2].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_hip_l[3].setFunction(new Constant(0.00000000000000000000));
	st_hip_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_l[4].setFunction(new Constant(0.00000000000000000000));
	st_hip_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_l[5].setFunction(new Constant(0.00000000000000000000));
	st_hip_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* hip_l;
	hip_l = new OpenSim::CustomJoint("hip_l", *pelvis, Vec3(-0.07069999999999999896, -0.06610000000000000597, -0.08350000000000000477), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_l);

	SpatialTransform st_knee_l;
	st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
	st_knee_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_knee_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_l[1].setFunction(new Constant(0.00000000000000000000));
	st_knee_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[2].setFunction(new Constant(0.00000000000000000000));
	st_knee_l[2].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[3].setFunction(new Constant(0.00000000000000000000));
	st_knee_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[4].setFunction(new Constant(0.00000000000000000000));
	st_knee_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[5].setFunction(new Constant(0.00000000000000000000));
	st_knee_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_l;
	knee_l = new OpenSim::CustomJoint("knee_l", *femur_l, Vec3(-0.00501372522863820051, -0.39658566558528168811, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_l);

	SpatialTransform st_ankle_l;
	st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
	st_ankle_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_ankle_l[0].setAxis(Vec3(0.10501354999999999718, 0.17402244999999999520, 0.97912631999999999444));
	st_ankle_l[1].setFunction(new Constant(0.00000000000000000000));
	st_ankle_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[2].setFunction(new Constant(0.00000000000000000000));
	st_ankle_l[2].setAxis(Vec3(0.97912631999999999444, 0.00000000000000000000, -0.10501354999999999718));
	st_ankle_l[3].setFunction(new Constant(0.00000000000000000000));
	st_ankle_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[4].setFunction(new Constant(0.00000000000000000000));
	st_ankle_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[5].setFunction(new Constant(0.00000000000000000000));
	st_ankle_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* ankle_l;
	ankle_l = new OpenSim::CustomJoint("ankle_l", *tibia_l, Vec3(0.00000000000000000000, -0.42999999999999999334, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_l);

	SpatialTransform st_subtalar_l;
	st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
	st_subtalar_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_subtalar_l[0].setAxis(Vec3(-0.78717961000000002958, -0.60474746000000001445, -0.12094949000000000672));
	st_subtalar_l[1].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[2].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_l[2].setAxis(Vec3(-0.12094949000000000672, 0.00000000000000000000, 0.78717961000000002958));
	st_subtalar_l[3].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[4].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[5].setFunction(new Constant(0.00000000000000000000));
	st_subtalar_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* subtalar_l;
	subtalar_l = new OpenSim::CustomJoint("subtalar_l", *talus_l, Vec3(-0.04877000000000000085, -0.04195000000000000118, -0.00791999999999999996), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_l);

	SpatialTransform st_mtp_l;
	st_mtp_l[0].setCoordinateNames(OpenSim::Array<std::string>("mtp_angle_l", 1, 1));
	st_mtp_l[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_mtp_l[0].setAxis(Vec3(0.58095439999999998193, 0.00000000000000000000, 0.81393610999999999045));
	st_mtp_l[1].setFunction(new Constant(0.00000000000000000000));
	st_mtp_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_mtp_l[2].setFunction(new Constant(0.00000000000000000000));
	st_mtp_l[2].setAxis(Vec3(0.81393610999999999045, 0.00000000000000000000, -0.58095439999999998193));
	st_mtp_l[3].setFunction(new Constant(0.00000000000000000000));
	st_mtp_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_mtp_l[4].setFunction(new Constant(0.00000000000000000000));
	st_mtp_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_mtp_l[5].setFunction(new Constant(0.00000000000000000000));
	st_mtp_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* mtp_l;
	mtp_l = new OpenSim::CustomJoint("mtp_l", *calcn_l, Vec3(0.17879999999999998672, -0.00200000000000000004, -0.00108000000000000001), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_mtp_l);

	SpatialTransform st_pbground;
	st_pbground[0].setCoordinateNames(OpenSim::Array<std::string>("global_rotZ", 1, 1));
	st_pbground[0].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[0].setAxis(Vec3(0,0,1));
	st_pbground[1].setCoordinateNames(OpenSim::Array<std::string>("global_rotX", 1, 1));
	st_pbground[1].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[1].setAxis(Vec3(1,0,0));
	st_pbground[2].setCoordinateNames(OpenSim::Array<std::string>("global_rotY", 1, 1));
	st_pbground[2].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[2].setAxis(Vec3(0,1,0));
	st_pbground[3].setCoordinateNames(OpenSim::Array<std::string>("global_tx", 1, 1));
	st_pbground[3].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[3].setAxis(Vec3(1,0,0));
	st_pbground[4].setCoordinateNames(OpenSim::Array<std::string>("global_ty", 1, 1));
	st_pbground[4].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[4].setAxis(Vec3(0,1,0));
	st_pbground[5].setCoordinateNames(OpenSim::Array<std::string>("global_tz", 1, 1));
	st_pbground[5].setFunction(new LinearFunction(1.0000, 0.0000));
	st_pbground[5].setAxis(Vec3(0,0,1));
	OpenSim::CustomJoint* pbground;
	pbground = new OpenSim::CustomJoint("pbground", model->getGround(), Vec3(0,0,0), Vec3(0,0,0), *punching_bag, Vec3(0,0,0), Vec3(0,0,0), st_pbground);



	model->addJoint(pelvisjnt);
	model->addJoint(spine_pelvis);
	model->addJoint(auxSTERNljnt);
	model->addJoint(auxSCAPljnt);
	model->addJoint(acromial_l);
	model->addJoint(ribcagejnt);
	model->addJoint(auxt1jnt);
	model->addJoint(aux7jnt);
	model->addJoint(aux6jnt);
	model->addJoint(aux5jnt);
	model->addJoint(aux4jnt);
	model->addJoint(aux3jnt);
	model->addJoint(aux2jnt);
	model->addJoint(aux1jnt);
	model->addJoint(jawjnt);
	model->addJoint(auxSTERNrjnt);
	model->addJoint(auxSCAPrjnt);
	model->addJoint(acromial_r);
	model->addJoint(elbow_l);
	model->addJoint(radioulnar_l);
	model->addJoint(radius_hand_l);
	model->addJoint(elbow_r);
	model->addJoint(radioulnar_r);
	model->addJoint(radius_hand_r);
	model->addJoint(hip_r);
	model->addJoint(knee_r);
	model->addJoint(ankle_r);
	model->addJoint(subtalar_r);
	model->addJoint(mtp_r);
	model->addJoint(hip_l);
	model->addJoint(knee_l);
	model->addJoint(ankle_l);
	model->addJoint(subtalar_l);
	model->addJoint(mtp_l);
	model->addJoint(pbground);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_heel_r;
	SmoothSphereHalfSpaceForce_heel_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_heel_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_heel_r_location(0.01000000000000000021, 0.00692291751087808875, -0.00499720000000000025);
	SmoothSphereHalfSpaceForce_heel_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_heel_r_location);
	double SmoothSphereHalfSpaceForce_heel_r_radius = (0.03200000000000000067);
	SmoothSphereHalfSpaceForce_heel_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_heel_r_radius);
	SmoothSphereHalfSpaceForce_heel_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_heel_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_heel_r->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_heel_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_heel_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_heel_r->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_heel_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_heel_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_hindfoot_lat_r;
	SmoothSphereHalfSpaceForce_hindfoot_lat_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_hindfoot_lat_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_r_location(0.05999999999999999778, 0.01192291751087808972, 0.02000100000000000142);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_hindfoot_lat_r_location);
	double SmoothSphereHalfSpaceForce_hindfoot_lat_r_radius = (0.03200000000000000067);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_hindfoot_lat_r_radius);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_hindfoot_lat_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_hindfoot_lat_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_forefoot_lat_r;
	SmoothSphereHalfSpaceForce_forefoot_lat_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_forefoot_lat_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_r_location(0.16420000000000001261, 0.00200000000000000004, 0.02029999999999999860);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_forefoot_lat_r_location);
	double SmoothSphereHalfSpaceForce_forefoot_lat_r_radius = (0.02299999999999999961);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_forefoot_lat_r_radius);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_forefoot_lat_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_forefoot_lat_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_forefoot_med_r;
	SmoothSphereHalfSpaceForce_forefoot_med_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_forefoot_med_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_r_location(0.16420000000000001261, 0.00200000000000000004, -0.01080000000000000057);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_forefoot_med_r_location);
	double SmoothSphereHalfSpaceForce_forefoot_med_r_radius = (0.02100000000000000130);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_forefoot_med_r_radius);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_forefoot_med_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_forefoot_med_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_hallux_r;
	SmoothSphereHalfSpaceForce_hallux_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_hallux_r", *toes_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_hallux_r_location(0.05315399999999999986, -0.00153854124456095573, -0.00341729999999999985);
	SmoothSphereHalfSpaceForce_hallux_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_hallux_r_location);
	double SmoothSphereHalfSpaceForce_hallux_r_radius = (0.01600000000000000033);
	SmoothSphereHalfSpaceForce_hallux_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_hallux_r_radius);
	SmoothSphereHalfSpaceForce_hallux_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_hallux_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_hallux_r->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hallux_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hallux_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_hallux_r->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_r->connectSocket_sphere_frame(*toes_r);
	SmoothSphereHalfSpaceForce_hallux_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_hallux_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_heel_l;
	SmoothSphereHalfSpaceForce_heel_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_heel_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_heel_l_location(0.01000000000000000021, 0.00692291751087808875, 0.00499720000000000025);
	SmoothSphereHalfSpaceForce_heel_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_heel_l_location);
	double SmoothSphereHalfSpaceForce_heel_l_radius = (0.03200000000000000067);
	SmoothSphereHalfSpaceForce_heel_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_heel_l_radius);
	SmoothSphereHalfSpaceForce_heel_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_heel_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_heel_l->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_heel_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_heel_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_heel_l->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_heel_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_heel_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_heel_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_hindfoot_lat_l;
	SmoothSphereHalfSpaceForce_hindfoot_lat_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_hindfoot_lat_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_l_location(0.05999999999999999778, 0.01192291751087808972, -0.02000100000000000142);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_hindfoot_lat_l_location);
	double SmoothSphereHalfSpaceForce_hindfoot_lat_l_radius = (0.03200000000000000067);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_hindfoot_lat_l_radius);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_hindfoot_lat_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_hindfoot_lat_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_forefoot_lat_l;
	SmoothSphereHalfSpaceForce_forefoot_lat_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_forefoot_lat_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_l_location(0.16420000000000001261, 0.00200000000000000004, -0.02029999999999999860);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_forefoot_lat_l_location);
	double SmoothSphereHalfSpaceForce_forefoot_lat_l_radius = (0.02299999999999999961);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_forefoot_lat_l_radius);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_forefoot_lat_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_forefoot_lat_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_forefoot_med_l;
	SmoothSphereHalfSpaceForce_forefoot_med_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_forefoot_med_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_l_location(0.16420000000000001261, 0.00200000000000000004, 0.01080000000000000057);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_forefoot_med_l_location);
	double SmoothSphereHalfSpaceForce_forefoot_med_l_radius = (0.02100000000000000130);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_forefoot_med_l_radius);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_forefoot_med_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_forefoot_med_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_forefoot_med_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_hallux_l;
	SmoothSphereHalfSpaceForce_hallux_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_hallux_l", *toes_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_hallux_l_location(0.05315399999999999986, -0.00153854124456095573, 0.00341729999999999985);
	SmoothSphereHalfSpaceForce_hallux_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_hallux_l_location);
	double SmoothSphereHalfSpaceForce_hallux_l_radius = (0.01600000000000000033);
	SmoothSphereHalfSpaceForce_hallux_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_hallux_l_radius);
	SmoothSphereHalfSpaceForce_hallux_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_hallux_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079633000000007392));
	SmoothSphereHalfSpaceForce_hallux_l->set_stiffness(10000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hallux_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_hallux_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_hallux_l->set_constant_contact_force(0.00001000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->set_hertz_smoothing(300.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->set_hunt_crossley_smoothing(50.00000000000000000000);
	SmoothSphereHalfSpaceForce_hallux_l->connectSocket_sphere_frame(*toes_l);
	SmoothSphereHalfSpaceForce_hallux_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_hallux_l);

	
	//-->
	// Define Bushing Forces
	Vec3 Kt_auxt1jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_auxt1jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_auxt1jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_auxt1jnt = Vec3(1.5, 1.5, 1.5);
	Rotation R;
	//OpenSim::PhysicalOffsetFrame* spine_offset;
	auto* spine_offset = new PhysicalOffsetFrame("spine_offset", *spine, Transform(R, Vec3(-0.063530500000000004, 0.596472, 0)));
	auto* cerv7_offset = new PhysicalOffsetFrame("cerv7_offset", *cerv7, Transform(R, Vec3(0.0068176, -0.0056364800000000001, 0)));
	model->addComponent(spine_offset);
	model->addComponent(cerv7_offset);
	auto* auxt1jnt_bushing = new BushingForce("auxt1jnt_bushing",
		"spine_offset", "cerv7_offset",
		Kt_auxt1jnt, Kr_auxt1jnt, Dt_auxt1jnt, Dr_auxt1jnt);
	model->addComponent(auxt1jnt_bushing);

	Vec3 Kt_aux7jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux7jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_aux7jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux7jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv7_offset_aux7 = new PhysicalOffsetFrame("cerv7_offset_aux7", *cerv7, Transform(R, Vec3(0.0132126, 0.0123395, 0)));
	auto* cerv6_offset = new PhysicalOffsetFrame("cerv6_offset", *cerv6, Transform(R, Vec3(0.0093881499999999996, - 0.0040516900000000002, 0)));
	model->addComponent(cerv7_offset_aux7);
	model->addComponent(cerv6_offset);
	auto* aux7jnt_bushing = new BushingForce("aux7jnt_bushing",
		"cerv7_offset_aux7", "cerv6_offset",
		Kt_aux7jnt, Kr_aux7jnt, Dt_aux7jnt, Dr_aux7jnt);
	model->addComponent(aux7jnt_bushing);

	Vec3 Kt_aux6jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux6jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_aux6jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux6jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv6_offset_aux6 = new PhysicalOffsetFrame("cerv6_offset_aux6", *cerv6, Transform(R, Vec3(0.0098177400000000001, 0.0097180700000000005, 0)));
	auto* cerv5_offset = new PhysicalOffsetFrame("cerv5_offset", *cerv5, Transform(R, Vec3(0.0056514299999999998, -0.0081233099999999999, 0)));
	model->addComponent(cerv6_offset_aux6);
	model->addComponent(cerv5_offset);
	auto* aux6jnt_bushing = new BushingForce("aux6jnt_bushing",
		"cerv6_offset_aux6", "cerv5_offset",
		Kt_aux6jnt, Kr_aux6jnt, Dt_aux6jnt, Dr_aux6jnt);
	model->addComponent(aux6jnt_bushing);

	Vec3 Kt_aux5jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux5jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_aux5jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux5jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv5_offset_aux5 = new PhysicalOffsetFrame("cerv5_offset_aux5", *cerv5, Transform(R, Vec3(0.0084023899999999992, 0.0074256000000000001, 0)));
	auto* cerv4_offset = new PhysicalOffsetFrame("cerv4_offset", *cerv4, Transform(R, Vec3(0.0050035599999999998, -0.0116816, 0)));
	model->addComponent(cerv5_offset_aux5);
	model->addComponent(cerv4_offset);
	auto* aux5jnt_bushing = new BushingForce("aux5jnt_bushing",
		"cerv5_offset_aux5", "cerv4_offset",
		Kt_aux5jnt, Kr_aux5jnt, Dt_aux5jnt, Dr_aux5jnt);
	model->addComponent(aux5jnt_bushing);

	Vec3 Kt_aux4jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux4jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_aux4jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux4jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv4_offset_aux4 = new PhysicalOffsetFrame("cerv4_offset_aux4", *cerv4, Transform(R, Vec3(0.0072960300000000002, 0.0063591100000000003, 0)));
	auto* cerv3_offset = new PhysicalOffsetFrame("cerv3_offset", *cerv3, Transform(R, Vec3(0.0046945700000000003, - 0.0089904600000000005, 0)));
	model->addComponent(cerv4_offset_aux4);
	model->addComponent(cerv3_offset);
	auto* aux4jnt_bushing = new BushingForce("aux4jnt_bushing",
		"cerv4_offset_aux4", "cerv3_offset",
		Kt_aux4jnt, Kr_aux4jnt, Dt_aux4jnt, Dr_aux4jnt);
	model->addComponent(aux4jnt_bushing);

	Vec3 Kt_aux3jnt = Vec3(75900, 23100000, 73000); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux3jnt = Vec3(1400, 4300, 1000);
	Vec3 Kr_aux3jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux3jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv3_offset_aux3 = new PhysicalOffsetFrame("cerv3_offset_aux3", *cerv3, Transform(R, Vec3(0.0052328100000000001, 0.0050434299999999998, 0)));
	auto* cerv2_offset = new PhysicalOffsetFrame("cerv2_offset", *cerv2, Transform(R, Vec3(0.0054122199999999997, - 0.0111932, 0)));
	model->addComponent(cerv3_offset_aux3);
	model->addComponent(cerv2_offset);
	auto* aux3jnt_bushing = new BushingForce("aux3jnt_bushing",
		"cerv3_offset_aux3", "cerv2_offset",
		Kt_aux3jnt, Kr_aux3jnt, Dt_aux3jnt, Dr_aux3jnt);
	model->addComponent(aux3jnt_bushing);

	Vec3 Kt_aux2jnt = Vec3(0); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux2jnt = Vec3(0);
	Vec3 Kr_aux2jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux2jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv2_offset_aux2 = new PhysicalOffsetFrame("cerv2_offset_aux2", *cerv2, Transform(R, Vec3(0.0068674399999999998, 0.021210300000000001, 0)));
	auto* cerv1_offset = new PhysicalOffsetFrame("cerv1_offset", *cerv1, Transform(R, Vec3(0.042311000000000001, 0.0035084700000000001, 0)));
	model->addComponent(cerv2_offset_aux2);
	model->addComponent(cerv1_offset);
	auto* aux2jnt_bushing = new BushingForce("aux2jnt_bushing",
		"cerv2_offset_aux2", "cerv1_offset",
		Kt_aux2jnt, Kr_aux2jnt, Dt_aux2jnt, Dr_aux2jnt);
	model->addComponent(aux2jnt_bushing);

	Vec3 Kt_aux1jnt = Vec3(0); //need to modify those 4 values (stiffness and damping values)
	Vec3 Dt_aux1jnt = Vec3(0);
	Vec3 Kr_aux1jnt = Vec3(18.91, 18.579999999999998, 24.059999999999999);
	Vec3 Dr_aux1jnt = Vec3(1.5, 1.5, 1.5);
	auto* cerv1_offset_aux1 = new PhysicalOffsetFrame("cerv1_offset_aux1", *cerv1, Transform(R, Vec3(0.20000000000000001, 0.20000000000000001, 0.20000000000000001)));
	auto* skull_offset = new PhysicalOffsetFrame("skull_offset", *skull, Transform(R, Vec3(0)));
	model->addComponent(cerv1_offset_aux1);
	model->addComponent(skull_offset);
	auto* aux1jnt_bushing = new BushingForce("aux1jnt_bushing",
		"cerv1_offset_aux1", "skull_offset",
		Kt_aux1jnt, Kr_aux1jnt, Dt_aux1jnt, Dr_aux1jnt);
	model->addComponent(aux1jnt_bushing);
	//<--

	// Initialize system.
	SimTK::State* state;
	state = new State(model->initSystem());

	std::cout << state->getNQ() << std::endl;

	// Read inputs.
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

	// States and controls.
	T ua[nCoordinates];
	Vector QsUs(NX);
	/// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	/// Controls
	/// OpenSim and Simbody have different state orders.
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < nCoordinates; ++i) ua[i] = u[indicesOSInSimbody[i]];

	// Set state variables and realize.
	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// get state variable names for debugging purposes
	Array<std::string> statevarnames = model->getStateVariableNames();
	for (int i = 0; i < statevarnames.size(); i++) {
		std::cout << i << " " << statevarnames[i] << std::endl;
	}

	// Compute residual forces.
	/// Set appliedMobilityForces (# mobilities).
	//Vector appliedMobilityForces(nCoordinates);
	Vector appliedMobilityForces;
	//appliedMobilityForces.setToZero();
	/// Set appliedBodyForces (# bodies + ground).
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1;
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Set gravity.
	Vec3 gravity(0);
	gravity[0] = 0.00000000000000000000;
	gravity[1] = -9.80664999999999942304;
	gravity[2] = 0.00000000000000000000;
	/// Add weights to appliedBodyForces.
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
			model->getBodySet().get(i).getMobilizedBodyIndex(),
			model->getBodySet().get(i).getMassCenter(),
			model->getBodySet().get(i).getMass() * gravity, appliedBodyForces);
	}
	/// Add contact forces to appliedBodyForces.
	Array<osim_double_adouble> Force_0 = SmoothSphereHalfSpaceForce_heel_r->getRecordValues(*state);
	SpatialVec GRF_0;
	GRF_0[0] = Vec3(Force_0[3], Force_0[4], Force_0[5]);
	GRF_0[1] = Vec3(Force_0[0], Force_0[1], Force_0[2]);
	int c_idx_0 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_0] += GRF_0;

	Array<osim_double_adouble> Force_1 = SmoothSphereHalfSpaceForce_hindfoot_lat_r->getRecordValues(*state);
	SpatialVec GRF_1;
	GRF_1[0] = Vec3(Force_1[3], Force_1[4], Force_1[5]);
	GRF_1[1] = Vec3(Force_1[0], Force_1[1], Force_1[2]);
	int c_idx_1 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_1] += GRF_1;

	Array<osim_double_adouble> Force_2 = SmoothSphereHalfSpaceForce_forefoot_lat_r->getRecordValues(*state);
	SpatialVec GRF_2;
	GRF_2[0] = Vec3(Force_2[3], Force_2[4], Force_2[5]);
	GRF_2[1] = Vec3(Force_2[0], Force_2[1], Force_2[2]);
	int c_idx_2 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_2] += GRF_2;

	Array<osim_double_adouble> Force_3 = SmoothSphereHalfSpaceForce_forefoot_med_r->getRecordValues(*state);
	SpatialVec GRF_3;
	GRF_3[0] = Vec3(Force_3[3], Force_3[4], Force_3[5]);
	GRF_3[1] = Vec3(Force_3[0], Force_3[1], Force_3[2]);
	int c_idx_3 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_3] += GRF_3;

	Array<osim_double_adouble> Force_4 = SmoothSphereHalfSpaceForce_hallux_r->getRecordValues(*state);
	SpatialVec GRF_4;
	GRF_4[0] = Vec3(Force_4[3], Force_4[4], Force_4[5]);
	GRF_4[1] = Vec3(Force_4[0], Force_4[1], Force_4[2]);
	int c_idx_4 = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_4] += GRF_4;

	Array<osim_double_adouble> Force_5 = SmoothSphereHalfSpaceForce_heel_l->getRecordValues(*state);
	SpatialVec GRF_5;
	GRF_5[0] = Vec3(Force_5[3], Force_5[4], Force_5[5]);
	GRF_5[1] = Vec3(Force_5[0], Force_5[1], Force_5[2]);
	int c_idx_5 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_5] += GRF_5;

	Array<osim_double_adouble> Force_6 = SmoothSphereHalfSpaceForce_hindfoot_lat_l->getRecordValues(*state);
	SpatialVec GRF_6;
	GRF_6[0] = Vec3(Force_6[3], Force_6[4], Force_6[5]);
	GRF_6[1] = Vec3(Force_6[0], Force_6[1], Force_6[2]);
	int c_idx_6 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_6] += GRF_6;

	Array<osim_double_adouble> Force_7 = SmoothSphereHalfSpaceForce_forefoot_lat_l->getRecordValues(*state);
	SpatialVec GRF_7;
	GRF_7[0] = Vec3(Force_7[3], Force_7[4], Force_7[5]);
	GRF_7[1] = Vec3(Force_7[0], Force_7[1], Force_7[2]);
	int c_idx_7 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_7] += GRF_7;

	Array<osim_double_adouble> Force_8 = SmoothSphereHalfSpaceForce_forefoot_med_l->getRecordValues(*state);
	SpatialVec GRF_8;
	GRF_8[0] = Vec3(Force_8[3], Force_8[4], Force_8[5]);
	GRF_8[1] = Vec3(Force_8[0], Force_8[1], Force_8[2]);
	int c_idx_8 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_8] += GRF_8;

	Array<osim_double_adouble> Force_9 = SmoothSphereHalfSpaceForce_hallux_l->getRecordValues(*state);
	SpatialVec GRF_9;
	GRF_9[0] = Vec3(Force_9[3], Force_9[4], Force_9[5]);
	GRF_9[1] = Vec3(Force_9[0], Force_9[1], Force_9[2]);
	int c_idx_9 = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_9] += GRF_9;

	//-->
	// compute and apply bushing forces
	// Note that the output will have 12 values,
	// the first 6 are the 3 forces applied to the body1 
	// and then 3 moments applied to the body 1 (at the origin of frame1), expressed in body1 frame
	// We want to insert them into appliedBody forces, which are SpatialVecs containing 
	// a vector of a moment applied to the body origin (not frame origin) and a vector of force, 
	// both expressed in ground, so we need to apply a transformation
	Array<osim_double_adouble> auxt1jnt_bushing_values = auxt1jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, spine, cerv7, spine_offset, cerv7_offset, auxt1jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux7jnt_bushing_values = aux7jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv7, cerv6, cerv7_offset, cerv6_offset, aux7jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux6jnt_bushing_values = aux6jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv6, cerv5, cerv6_offset, cerv5_offset, aux6jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux5jnt_bushing_values = aux5jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv5, cerv4, cerv5_offset, cerv4_offset, aux5jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux4jnt_bushing_values = aux4jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv4, cerv3, cerv4_offset, cerv3_offset, aux4jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux3jnt_bushing_values = aux3jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv3, cerv2, cerv3_offset, cerv2_offset, aux3jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux2jnt_bushing_values = aux2jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv2, cerv1, cerv2_offset, cerv1_offset, aux2jnt_bushing_values, appliedBodyForces);
	Array<osim_double_adouble> aux1jnt_bushing_values = aux1jnt_bushing->getRecordValues(*state);
	expressBodyForcesInGround(state, cerv1, skull, cerv1_offset, skull_offset, aux1jnt_bushing_values, appliedBodyForces);
	//<--
	std::cout << appliedBodyForces[cerv1->getMobilizedBodyIndex()] << std::endl;
	std::cout << "appliedMobilityForces.size()" << appliedMobilityForces.size() << std::endl;
	std::cout << "NU--> " << model->getMatterSubsystem().getNU(*state) << std::endl;

	/// knownUdot.
	const auto& matter = model->getMatterSubsystem();
	/*Vector knownUdot(nCoordinates);*/
	const int nu = matter.getNU(*state);
	SimTK::Vector knownUdot(nu);

	knownUdot.setToZero();
	for (int i = 0; i < nCoordinates; ++i) knownUdot[i] = ua[i];

	/// Calculate residual forces.
	Vector residualMobilityForces(nCoordinates);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);

	/// Ground reaction forces.
	SpatialVec GRF_r;
	SpatialVec GRF_l;
	GRF_r.setToZero();
	GRF_l.setToZero();

	GRF_r += GRF_0;
	GRF_r += GRF_1;
	GRF_r += GRF_2;
	GRF_r += GRF_3;
	GRF_r += GRF_4;
	GRF_l += GRF_5;
	GRF_l += GRF_6;
	GRF_l += GRF_7;
	GRF_l += GRF_8;
	GRF_l += GRF_9;

	/// Ground reaction moments.
	SpatialVec GRM_r;
	SpatialVec GRM_l;
	GRM_r.setToZero();
	GRM_l.setToZero();
	Vec3 normal(0, 1, 0);

	SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
	Vec3 SmoothSphereHalfSpaceForce_heel_r_location_G = calcn_r->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_heel_r_location);
	Vec3 SmoothSphereHalfSpaceForce_heel_r_locationCP_G = SmoothSphereHalfSpaceForce_heel_r_location_G - SmoothSphereHalfSpaceForce_heel_r_radius * normal;
	Vec3 locationCP_G_adj_0 = SmoothSphereHalfSpaceForce_heel_r_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_heel_r_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_heel_r_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_0, *calcn_r);
	Vec3 GRM_0 = (TR_GB_calcn_r * SmoothSphereHalfSpaceForce_heel_r_locationCP_B) % GRF_0[1];
	GRM_r += GRM_0;

	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_r_location_G = calcn_r->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_hindfoot_lat_r_location);
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_r_locationCP_G = SmoothSphereHalfSpaceForce_hindfoot_lat_r_location_G - SmoothSphereHalfSpaceForce_hindfoot_lat_r_radius * normal;
	Vec3 locationCP_G_adj_1 = SmoothSphereHalfSpaceForce_hindfoot_lat_r_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_hindfoot_lat_r_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_r_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_1, *calcn_r);
	Vec3 GRM_1 = (TR_GB_calcn_r * SmoothSphereHalfSpaceForce_hindfoot_lat_r_locationCP_B) % GRF_1[1];
	GRM_r += GRM_1;

	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_r_location_G = calcn_r->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_forefoot_lat_r_location);
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_r_locationCP_G = SmoothSphereHalfSpaceForce_forefoot_lat_r_location_G - SmoothSphereHalfSpaceForce_forefoot_lat_r_radius * normal;
	Vec3 locationCP_G_adj_2 = SmoothSphereHalfSpaceForce_forefoot_lat_r_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_forefoot_lat_r_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_r_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_2, *calcn_r);
	Vec3 GRM_2 = (TR_GB_calcn_r * SmoothSphereHalfSpaceForce_forefoot_lat_r_locationCP_B) % GRF_2[1];
	GRM_r += GRM_2;

	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_r_location_G = calcn_r->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_forefoot_med_r_location);
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_r_locationCP_G = SmoothSphereHalfSpaceForce_forefoot_med_r_location_G - SmoothSphereHalfSpaceForce_forefoot_med_r_radius * normal;
	Vec3 locationCP_G_adj_3 = SmoothSphereHalfSpaceForce_forefoot_med_r_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_forefoot_med_r_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_r_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_3, *calcn_r);
	Vec3 GRM_3 = (TR_GB_calcn_r * SmoothSphereHalfSpaceForce_forefoot_med_r_locationCP_B) % GRF_3[1];
	GRM_r += GRM_3;

	SimTK::Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);
	Vec3 SmoothSphereHalfSpaceForce_hallux_r_location_G = toes_r->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_hallux_r_location);
	Vec3 SmoothSphereHalfSpaceForce_hallux_r_locationCP_G = SmoothSphereHalfSpaceForce_hallux_r_location_G - SmoothSphereHalfSpaceForce_hallux_r_radius * normal;
	Vec3 locationCP_G_adj_4 = SmoothSphereHalfSpaceForce_hallux_r_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_hallux_r_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_hallux_r_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_4, *toes_r);
	Vec3 GRM_4 = (TR_GB_toes_r * SmoothSphereHalfSpaceForce_hallux_r_locationCP_B) % GRF_4[1];
	GRM_r += GRM_4;

	SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
	Vec3 SmoothSphereHalfSpaceForce_heel_l_location_G = calcn_l->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_heel_l_location);
	Vec3 SmoothSphereHalfSpaceForce_heel_l_locationCP_G = SmoothSphereHalfSpaceForce_heel_l_location_G - SmoothSphereHalfSpaceForce_heel_l_radius * normal;
	Vec3 locationCP_G_adj_5 = SmoothSphereHalfSpaceForce_heel_l_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_heel_l_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_heel_l_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_5, *calcn_l);
	Vec3 GRM_5 = (TR_GB_calcn_l * SmoothSphereHalfSpaceForce_heel_l_locationCP_B) % GRF_5[1];
	GRM_l += GRM_5;

	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_l_location_G = calcn_l->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_hindfoot_lat_l_location);
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_l_locationCP_G = SmoothSphereHalfSpaceForce_hindfoot_lat_l_location_G - SmoothSphereHalfSpaceForce_hindfoot_lat_l_radius * normal;
	Vec3 locationCP_G_adj_6 = SmoothSphereHalfSpaceForce_hindfoot_lat_l_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_hindfoot_lat_l_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_l_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_6, *calcn_l);
	Vec3 GRM_6 = (TR_GB_calcn_l * SmoothSphereHalfSpaceForce_hindfoot_lat_l_locationCP_B) % GRF_6[1];
	GRM_l += GRM_6;

	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_l_location_G = calcn_l->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_forefoot_lat_l_location);
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_l_locationCP_G = SmoothSphereHalfSpaceForce_forefoot_lat_l_location_G - SmoothSphereHalfSpaceForce_forefoot_lat_l_radius * normal;
	Vec3 locationCP_G_adj_7 = SmoothSphereHalfSpaceForce_forefoot_lat_l_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_forefoot_lat_l_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_l_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_7, *calcn_l);
	Vec3 GRM_7 = (TR_GB_calcn_l * SmoothSphereHalfSpaceForce_forefoot_lat_l_locationCP_B) % GRF_7[1];
	GRM_l += GRM_7;

	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_l_location_G = calcn_l->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_forefoot_med_l_location);
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_l_locationCP_G = SmoothSphereHalfSpaceForce_forefoot_med_l_location_G - SmoothSphereHalfSpaceForce_forefoot_med_l_radius * normal;
	Vec3 locationCP_G_adj_8 = SmoothSphereHalfSpaceForce_forefoot_med_l_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_forefoot_med_l_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_l_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_8, *calcn_l);
	Vec3 GRM_8 = (TR_GB_calcn_l * SmoothSphereHalfSpaceForce_forefoot_med_l_locationCP_B) % GRF_8[1];
	GRM_l += GRM_8;

	SimTK::Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);
	Vec3 SmoothSphereHalfSpaceForce_hallux_l_location_G = toes_l->findStationLocationInGround(*state, SmoothSphereHalfSpaceForce_hallux_l_location);
	Vec3 SmoothSphereHalfSpaceForce_hallux_l_locationCP_G = SmoothSphereHalfSpaceForce_hallux_l_location_G - SmoothSphereHalfSpaceForce_hallux_l_radius * normal;
	Vec3 locationCP_G_adj_9 = SmoothSphereHalfSpaceForce_hallux_l_locationCP_G - 0.5 * SmoothSphereHalfSpaceForce_hallux_l_locationCP_G[1] * normal;
	Vec3 SmoothSphereHalfSpaceForce_hallux_l_locationCP_B = model->getGround().findStationLocationInAnotherFrame(*state, locationCP_G_adj_9, *toes_l);
	Vec3 GRM_9 = (TR_GB_toes_l * SmoothSphereHalfSpaceForce_hallux_l_locationCP_B) % GRF_9[1];
	GRM_l += GRM_9;

	/// Contact spheres deformation power.
	Vec3 SmoothSphereHalfSpaceForce_heel_r_velocity_G = calcn_r->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_heel_r_location);
	osim_double_adouble P_HC_y_0 = SmoothSphereHalfSpaceForce_heel_r_velocity_G[1] * GRF_0[1][1];
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_r_velocity_G = calcn_r->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_hindfoot_lat_r_location);
	osim_double_adouble P_HC_y_1 = SmoothSphereHalfSpaceForce_hindfoot_lat_r_velocity_G[1] * GRF_1[1][1];
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_r_velocity_G = calcn_r->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_forefoot_lat_r_location);
	osim_double_adouble P_HC_y_2 = SmoothSphereHalfSpaceForce_forefoot_lat_r_velocity_G[1] * GRF_2[1][1];
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_r_velocity_G = calcn_r->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_forefoot_med_r_location);
	osim_double_adouble P_HC_y_3 = SmoothSphereHalfSpaceForce_forefoot_med_r_velocity_G[1] * GRF_3[1][1];
	Vec3 SmoothSphereHalfSpaceForce_hallux_r_velocity_G = toes_r->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_hallux_r_location);
	osim_double_adouble P_HC_y_4 = SmoothSphereHalfSpaceForce_hallux_r_velocity_G[1] * GRF_4[1][1];
	Vec3 SmoothSphereHalfSpaceForce_heel_l_velocity_G = calcn_l->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_heel_l_location);
	osim_double_adouble P_HC_y_5 = SmoothSphereHalfSpaceForce_heel_l_velocity_G[1] * GRF_5[1][1];
	Vec3 SmoothSphereHalfSpaceForce_hindfoot_lat_l_velocity_G = calcn_l->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_hindfoot_lat_l_location);
	osim_double_adouble P_HC_y_6 = SmoothSphereHalfSpaceForce_hindfoot_lat_l_velocity_G[1] * GRF_6[1][1];
	Vec3 SmoothSphereHalfSpaceForce_forefoot_lat_l_velocity_G = calcn_l->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_forefoot_lat_l_location);
	osim_double_adouble P_HC_y_7 = SmoothSphereHalfSpaceForce_forefoot_lat_l_velocity_G[1] * GRF_7[1][1];
	Vec3 SmoothSphereHalfSpaceForce_forefoot_med_l_velocity_G = calcn_l->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_forefoot_med_l_location);
	osim_double_adouble P_HC_y_8 = SmoothSphereHalfSpaceForce_forefoot_med_l_velocity_G[1] * GRF_8[1][1];
	Vec3 SmoothSphereHalfSpaceForce_hallux_l_velocity_G = toes_l->findStationVelocityInGround(*state, SmoothSphereHalfSpaceForce_hallux_l_location);
	osim_double_adouble P_HC_y_9 = SmoothSphereHalfSpaceForce_hallux_l_velocity_G[1] * GRF_9[1][1];
	/// Outputs.
	/// Residual forces (OpenSim and Simbody have different state orders).
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < nCoordinates; ++i) res[0][i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
	/// Ground reaction forces.
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 0] = value<T>(GRF_r[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 3] = value<T>(GRF_l[1][i]);
	/// Separate Ground reaction forces.
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 6] = value<T>(GRF_0[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 9] = value<T>(GRF_1[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 12] = value<T>(GRF_2[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 15] = value<T>(GRF_3[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 18] = value<T>(GRF_4[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 21] = value<T>(GRF_5[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 24] = value<T>(GRF_6[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 27] = value<T>(GRF_7[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 30] = value<T>(GRF_8[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 33] = value<T>(GRF_9[1][i]);
	/// Ground reaction moments.
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 36] = value<T>(GRM_r[1][i]);
	for (int i = 0; i < 3; ++i) res[0][i + nCoordinates + 39] = value<T>(GRM_l[1][i]);
	/// Contact spheres deformation power.
	res[0][nCoordinates + 42] = value<T>(P_HC_y_0);
	res[0][nCoordinates + 43] = value<T>(P_HC_y_1);
	res[0][nCoordinates + 44] = value<T>(P_HC_y_2);
	res[0][nCoordinates + 45] = value<T>(P_HC_y_3);
	res[0][nCoordinates + 46] = value<T>(P_HC_y_4);
	res[0][nCoordinates + 47] = value<T>(P_HC_y_5);
	res[0][nCoordinates + 48] = value<T>(P_HC_y_6);
	res[0][nCoordinates + 49] = value<T>(P_HC_y_7);
	res[0][nCoordinates + 50] = value<T>(P_HC_y_8);
	res[0][nCoordinates + 51] = value<T>(P_HC_y_9);

	return 0;
}

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

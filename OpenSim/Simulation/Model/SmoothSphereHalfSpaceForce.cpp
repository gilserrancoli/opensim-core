/* -------------------------------------------------------------------------- *
 *                   OpenSim: SmoothSphereHalfSpaceForce.cpp                  *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2017-19 Stanford University and the Authors                  *
 *                                                                            *
 * Author(s): Antoine Falisse, Gil Serrancoli                                 *
 * Contributors: Peter Eastman                                                *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0          *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SmoothSphereHalfSpaceForce.h"

#include <simbody/internal/SmoothSphereHalfSpaceForce.h>

#include <OpenSim/Simulation/Model/Model.h>

using namespace OpenSim;

//=============================================================================
//  SMOOTH SPHERE HALF SPACE FORCE
//=============================================================================
// Uses default (compiler-generated) destructor, copy constructor, copy
// assignment operator.

// Default constructor.
SmoothSphereHalfSpaceForce::SmoothSphereHalfSpaceForce() {
    constructProperties();
}

// Take over ownership of supplied object.
SmoothSphereHalfSpaceForce::SmoothSphereHalfSpaceForce(const std::string& name,
        const PhysicalFrame& contactSphereBodyFrame,
        const PhysicalFrame& contactHalfSpaceBodyFrame) {
    setName(name);
    connectSocket_sphere_frame(contactSphereBodyFrame);
    connectSocket_half_space_frame(contactHalfSpaceBodyFrame);

    constructProperties();
}

// TODO: unclear why it does not work
//SmoothSphereHalfSpaceForce::SmoothSphereHalfSpaceForce(const std::string& name,
//        const PhysicalFrame& contactSphereBodyFrame,
//        SimTK::Vec3 contactSphereLocation,
//        osim_double_adouble contactSphereRadius,
//        const PhysicalFrame& contactHalfSpaceBodyFrame,
//        SimTK::Vec3 contactHalfSpaceLocation,
//        SimTK::Vec3 contactHalfSpaceOrientation) {
//    setName(name);
//	set_contact_sphere_location(contactSphereLocation);
//	set_contact_sphere_radius(contactSphereRadius);
//	set_contact_half_space_location(contactHalfSpaceLocation);
//	set_contact_half_space_orientation(contactHalfSpaceOrientation);
//
//    connectSocket_sphere_frame(contactSphereBodyFrame);
//    connectSocket_half_space_frame(contactHalfSpaceBodyFrame);  
//
//    constructProperties();
//}

void SmoothSphereHalfSpaceForce::extendAddToSystem(
        SimTK::MultibodySystem& system) const {

    Super::extendAddToSystem(system);

    const SimTK::Vec3& contactSphereLocation = get_contact_sphere_location();
    const osim_double_adouble& contactSphereRadius =
        get_contact_sphere_radius();

    osim_double_adouble stiffness = get_stiffness();
    osim_double_adouble dissipation = get_dissipation();
    osim_double_adouble staticFriction = get_static_friction();
    osim_double_adouble dynamicFriction = get_dynamic_friction();
    osim_double_adouble viscousFriction = get_viscous_friction();
    osim_double_adouble transitionVelocity = get_transition_velocity();
    osim_double_adouble cf = get_constant_contact_force();
    osim_double_adouble bd = get_hertz_smoothing();
    osim_double_adouble bv = get_hunt_crossley_smoothing();

    SimTK::SmoothSphereHalfSpaceForce force(_model->updForceSubsystem());

    const auto& sphereFrame = getConnectee<PhysicalFrame>("sphere_frame");
    const auto& halfSpaceFrame =
        getConnectee<PhysicalFrame>("half_space_frame");

    force.setStiffness(stiffness);
    force.setDissipation(dissipation);
    force.setStaticFriction(staticFriction);
    force.setDynamicFriction(dynamicFriction);
    force.setViscousFriction(viscousFriction);
    force.setTransitionVelocity(transitionVelocity);
    force.setConstantContactForce(cf);
    force.setHertzSmoothing(bd);
    force.setHuntCrossleySmoothing(bv);

    force.setContactSphereBody(sphereFrame.getMobilizedBody());
    force.setContactSphereLocationInBody(
            sphereFrame.findTransformInBaseFrame() *
            contactSphereLocation);
    force.setContactSphereRadius(contactSphereRadius);

    force.setContactHalfSpaceBody(halfSpaceFrame.getMobilizedBody());

    SimTK::Transform halfSpaceFrameTransform(
        SimTK::Rotation(SimTK::BodyRotationSequence,
            get_contact_half_space_orientation()[0], SimTK::XAxis,
            get_contact_half_space_orientation()[1], SimTK::YAxis,
            get_contact_half_space_orientation()[2], SimTK::ZAxis),
        get_contact_half_space_location());

    force.setContactHalfSpaceFrame(
            halfSpaceFrame.findTransformInBaseFrame() *
            halfSpaceFrameTransform);

    auto* mutableThis = const_cast<SmoothSphereHalfSpaceForce*>(this);
    mutableThis->_index = force.getForceIndex();
}

//void OpenSim::SmoothSphereHalfSpaceForce::extendRealizeInstance(
//        const SimTK::State& state) const {
//    Super::extendRealizeInstance(state);
//    if (!getProperty_force_visualization_scale_factor().empty()) {
//        m_forceVizScaleFactor = get_force_visualization_scale_factor();
//    } else {
//        const Model& model = getModel();
//        const osim_double_adouble mass = model.getTotalMass(state);
//        const osim_double_adouble weight = mass * model.getGravity().norm();
//        m_forceVizScaleFactor = 1 / weight;
//    }
//}

void SmoothSphereHalfSpaceForce::constructProperties() {
    constructProperty_stiffness(1.0);
    constructProperty_dissipation(0.0);
    constructProperty_static_friction(0.0);
    constructProperty_dynamic_friction(0.0);
    constructProperty_viscous_friction(0.0);
    constructProperty_transition_velocity(0.01);
    constructProperty_constant_contact_force(1e-5);
    constructProperty_hertz_smoothing(300.0);
    constructProperty_hunt_crossley_smoothing(50.0);
    //constructProperty_force_visualization_radius(0.01);
    //constructProperty_force_visualization_scale_factor();
    constructProperty_contact_sphere_location(SimTK::Vec3(0));
    constructProperty_contact_sphere_radius(0.0);
    constructProperty_contact_half_space_location(SimTK::Vec3(0));
    constructProperty_contact_half_space_orientation(
        SimTK::Vec3(0,0,-0.5*SimTK::Pi));
}

//=============================================================================
//  REPORTING
//=============================================================================
// Provide names of the quantities (column labels) of the force value(s)
OpenSim::Array<std::string>
SmoothSphereHalfSpaceForce::getRecordLabels() const {
    OpenSim::Array<std::string> labels("");

    labels.append(getName() + ".Sphere" + ".force.X");
    labels.append(getName() + ".Sphere" + ".force.Y");
    labels.append(getName() + ".Sphere" + ".force.Z");
    labels.append(getName() + ".Sphere" + ".torque.X");
    labels.append(getName() + ".Sphere" + ".torque.Y");
    labels.append(getName() + ".Sphere" + ".torque.Z");

    labels.append(getName() + ".HalfSpace" + ".force.X");
    labels.append(getName() + ".HalfSpace" + ".force.Y");
    labels.append(getName() + ".HalfSpace" + ".force.Z");
    labels.append(getName() + ".HalfSpace" + ".torque.X");
    labels.append(getName() + ".HalfSpace" + ".torque.Y");
    labels.append(getName() + ".HalfSpace" + ".torque.Z");

    return labels;
}

// Provide the value(s) to be reported that correspond to the labels
OpenSim::Array<osim_double_adouble> SmoothSphereHalfSpaceForce::getRecordValues(
        const SimTK::State& state) const {

    OpenSim::Array<osim_double_adouble> values(1);

    const auto& sphereFrame = getConnectee<PhysicalFrame>("sphere_frame");
    const auto sphereIdx = sphereFrame.getMobilizedBodyIndex();

    const auto& halfSpaceFrame =
        getConnectee<PhysicalFrame>("half_space_frame");
    const auto halfSpaceIdx = halfSpaceFrame.getMobilizedBodyIndex();

    const Model& model = getModel();
    const auto& forceSubsys = model.getForceSubsystem();
    const SimTK::Force& abstractForce = forceSubsys.getForce(_index);
    const auto& simtkForce =
            static_cast<const SimTK::SmoothSphereHalfSpaceForce&>(
                    abstractForce);

    SimTK::Vector_<SimTK::SpatialVec> bodyForces(0);
    SimTK::Vector_<SimTK::Vec3> particleForces(0);
    SimTK::Vector mobilityForces(0);

    simtkForce.calcForceContribution(
            state, bodyForces, particleForces, mobilityForces);

    // On sphere
    const auto& thisBodyForce1 = bodyForces(sphereIdx);
    SimTK::Vec3 forces1 = thisBodyForce1[1];
    SimTK::Vec3 torques1 = thisBodyForce1[0];
    values.append(3, &forces1[0]);
    values.append(3, &torques1[0]);

    // On plane
    const auto& thisBodyForce2 = bodyForces(halfSpaceIdx);
    SimTK::Vec3 forces2 = thisBodyForce2[1];
    SimTK::Vec3 torques2 = thisBodyForce2[0];
    values.append(3, &forces2[0]);
    values.append(3, &torques2[0]);

    return values;
}

//void SmoothSphereHalfSpaceForce::generateDecorations(bool fixed,
//        const ModelDisplayHints& hints, const SimTK::State& state,
//        SimTK::Array_<SimTK::DecorativeGeometry>& geometry) const {
//    Super::generateDecorations(fixed, hints, state, geometry);
//
//    if (!fixed && (state.getSystemStage() >= SimTK::Stage::Dynamics) &&
//            hints.get_show_forces()) {
//        // Get the underlying SimTK force element.
//        const Model& model = getModel();
//        const auto& forceSubsystem = model.getForceSubsystem();
//        const SimTK::Force& abstractForce = forceSubsystem.getForce(_index);
//        const auto& simtkForce =
//                static_cast<const SimTK::SmoothSphereHalfSpaceForce&>(
//                        abstractForce);
//
//        // Compute the body forces.
//        SimTK::Vector_<SimTK::SpatialVec> bodyForces(0);
//        SimTK::Vector_<SimTK::Vec3> particleForces(0);
//        SimTK::Vector mobilityForces(0);
//        simtkForce.calcForceContribution(
//                state, bodyForces, particleForces, mobilityForces);
//
//        // Get the index to the associated contact sphere.
//        const auto& sphere = getConnectee<ContactSphere>("sphere");
//        const auto& sphereIdx = sphere.getFrame().getMobilizedBodyIndex();
//
//        // Get the translational force for the contact sphere associated with
//        // this force element.
//        const auto& sphereForce = bodyForces(sphereIdx)[1];
//
//        // Scale the contact force vector and compute the cylinder length.
//        const auto& scaledContactForce =
//                m_forceVizScaleFactor * sphereForce;
//        const SimTK::Real length(scaledContactForce.norm());
//
//        // Compute the force visualization transform.
//        const SimTK::Vec3 contactSpherePosition =
//                sphere.getFrame().findStationLocationInGround(
//                        state, sphere.get_location());
//        const SimTK::Transform forceVizTransform(
//                SimTK::Rotation(SimTK::UnitVec3(scaledContactForce),
//                        SimTK::YAxis),
//                contactSpherePosition + scaledContactForce / 2.0);
//
//        // Construct the force decoration and add it to the list of geometries.
//        SimTK::DecorativeCylinder forceViz(get_force_visualization_radius(),
//                0.5 * length);
//        forceViz.setTransform(forceVizTransform);
//        forceViz.setColor(SimTK::Vec3(0.0, 0.6, 0.0));
//        geometry.push_back(forceViz);
//    }
//}
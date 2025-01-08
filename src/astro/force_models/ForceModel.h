//
// Created by Douglas Garza on 8/25/24.
//

#pragma once
#include <Eigen/Dense>
#include "../SpaceVehicle/SpaceVehicle.h"
#include "ExponentialDragLookupTable.h"
#include "../../math/UniversalConstants.h"
#include <cmath>


class ForceModel {
public:
    virtual ~ForceModel() = default;
    virtual std::string getName() const = 0;
    virtual Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const = 0;
    Eigen::Vector3d accelPerturb;
};

class TwoBody : public ForceModel {
public:
    TwoBody();
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
};

class J2Force : public ForceModel {
public:
    J2Force();
    ~J2Force() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
};

class DragExpForce : public ForceModel {
public:
    DragExpForce();
    ~DragExpForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
};


class SimpleSRPForce : public ForceModel {
public:
    SimpleSRPForce();
    ~SimpleSRPForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
};


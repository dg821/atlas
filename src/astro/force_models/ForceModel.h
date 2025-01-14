//
// Created by Douglas Garza on 8/25/24.
//

#pragma once
#include <Eigen/Dense>
#include "../SpaceVehicle/SpaceVehicle.h"
#include "ExponentialDragLookupTable.h"
#include "../../math/UniversalConstants.h"
#include <cmath>
#include "../../input_ouput/planetary_ephemerides/meeusEphemeris.h"
#include "../geodetic_model/GeodeticModel.h"
#include <vector>
#include <stdexcept>


// Each new force model inherits the methods of the abstract class ForceModel

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


class LunisolarThirdBodyForce : public ForceModel {
public:
    LunisolarThirdBodyForce();
    ~LunisolarThirdBodyForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
};



class NonSphericalGravity : public ForceModel {
private:
    static constexpr int DEFAULT_ORDER = 4;
    static constexpr int MAX_ORDER = 50;

    int order_n;  // max degree
    int order_m;  // max order
    std::vector<std::vector<double>> Cnm;  // normalized cosine coefficients
    std::vector<std::vector<double>> Snm;  // normalized sine coefficients

    // stored calculations to avoid excessive computations
    std::vector<std::vector<double>> Pnm;      // associated Legendre polynomials
    std::vector<std::vector<double>> dPnm;     // derivatives of Legendre polynomials
    std::vector<double> factorials;            // pre-computed factorials

    void initializeCoefficients();
    void computeLegendrePolynomials(double phi);

public:
    NonSphericalGravity(int n = DEFAULT_ORDER, int m = DEFAULT_ORDER);
    Eigen::Vector3d computeAcceleration(double t, const SpaceVehicle& sv) const override;
    void setOrders(int n, int m);
    std::string getName() const override;
};
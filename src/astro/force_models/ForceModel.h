//
// Created by Douglas Garza on 8/25/24.
//

#pragma once
#include <Eigen/Dense>
#include "../SpaceVehicle/SpaceVehicle.h"
#include "ExponentialDragLookupTable.h"
#include "../../math/UniversalConstants.h"
#include "../../math/utilities/NumericalDiff.h"
#include <cmath>
#include "../../input_output/planetary_ephemerides/meeusEphemeris.h"
#include "../geodetic_model/geodeticModel.h"
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include "../../math/utilities/mathFunctions.h"

// Each new force model inherits the methods of the abstract class ForceModel

class ForceModel {
public:
    virtual ~ForceModel() = default;
    virtual std::string getName() const = 0;
    virtual Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const = 0;
    Eigen::Vector3d accelPerturb;
};

class TwoBody : public ForceModel {
public:
    TwoBody();
    ~TwoBody() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
};




class DragExpForce : public ForceModel {
public:
    DragExpForce();
    ~DragExpForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
};


class SimpleSRPForce : public ForceModel {
public:
    SimpleSRPForce();
    ~SimpleSRPForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
};


class LunisolarThirdBodyForce : public ForceModel {
public:
    LunisolarThirdBodyForce();
    ~LunisolarThirdBodyForce() override;
    std::string getName() const override;
    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
};


class NonSphericalGravity : public ForceModel {
private:
    friend class J2;

    int degree_l;
    int order_m;
    static constexpr int DEFAULT_ORDER = 4;
    static constexpr int DEFAULT_DEGREE = 4;
    static constexpr int MAX_ORDER = 20;
    static constexpr int MAX_DEGREE = 20;
    static constexpr double numDiffTol = 1e-4;

    std::string coefficients_path;

    // Gravity coefficients
    std::vector<std::vector<double>> Clm;
    std::vector<std::vector<double>> Slm;

    void initializeCoefficients();
    bool loadCoefficientsFromFile();

    std::vector<std::vector<double>> computeLegendrePolynomials(double phi) const;
    double normalizeLegendreFunction(double Plm, double l, double m, double k) const;

public:
    NonSphericalGravity(int l = DEFAULT_DEGREE, int m = DEFAULT_ORDER,
                       const std::string& coeff_path = "../../../data/egm96_20x20.csv");
    ~NonSphericalGravity() override;

    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
    void setDegreeOrder(int l, int m);
    std::string getCoefficientsPath() const { return coefficients_path; }

    std::string getName() const override;
};


class J2 : public ForceModel {
private:
    NonSphericalGravity nsg;

public:
    J2() : nsg(2, 0) {};

    Eigen::Vector3d computeAcceleration(double t, const Eigen::Vector3d& r, const Eigen::Vector3d& v, const SpaceVehicle& sv) const override;
    std::string getName() const override;

};
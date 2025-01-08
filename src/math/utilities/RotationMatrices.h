//
// Created by Douglas Garza on 8/25/24.
//

#pragma once

#include <Eigen/Dense>

class RotationMatrices {
  public:
    static Eigen::Matrix3d rotX(double theta);
    static Eigen::Matrix3d rotY(double theta);
    static Eigen::Matrix3d rotZ(double theta);
};

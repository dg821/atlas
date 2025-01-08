//
// Created by Douglas Garza on 8/25/24.
//

#include "RotationMatrices.h"
#include <Eigen/Dense>
#include <cmath>

Eigen::Matrix3d RotationMatrices::rotX(const double theta) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);

  Eigen::Matrix3d rot_mat;
  
  rot_mat(0,0) = 1.0;   rot_mat(0,1) = 0.0;   rot_mat(0,2) = 0.0;
  rot_mat(1,0) = 0.0;   rot_mat(1,1) = c;   rot_mat(1,2) = s;
  rot_mat(2,0) = 0.0;   rot_mat(2,1) = -s;  rot_mat(2,2) = c;

  return rot_mat;
}


Eigen::Matrix3d RotationMatrices::rotY(const double theta) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);

  Eigen::Matrix3d rot_mat;

  rot_mat(0,0) = c;   rot_mat(0,1) = 0.0;   rot_mat(0,2) = -s;
  rot_mat(1,0) = 0.0;   rot_mat(1,1) = 1.0;   rot_mat(1,2) = 0.0;
  rot_mat(2,0) = s;   rot_mat(2,1) = 0.0;  rot_mat(2,2) = c;

  return rot_mat;
}


Eigen::Matrix3d RotationMatrices::rotZ(const double theta) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);

  Eigen::Matrix3d rot_mat;

  rot_mat(0,0) = c;   rot_mat(0,1) = s;   rot_mat(0,2) = 0.0;
  rot_mat(1,0) = -s;   rot_mat(1,1) = c;   rot_mat(1,2) = 0.0;
  rot_mat(2,0) = 0.0;   rot_mat(2,1) = 0.0;  rot_mat(2,2) = 1.0;

  return rot_mat;
}
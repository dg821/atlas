//
// Created by Douglas Garza on 10/15/24.
//

#pragma once

#include "../../math/UniversalConstants.h"
#include <cmath>
#include <vector>
#include <functional>

class FundamentalArguments {
public:

  static double getMeanAnomalyMoon(double t_tt) {
    return (485868.249036 + 1717915923.2178 * t_tt + 31.8792 * pow(t_tt, 2) +
           0.051635 * pow(t_tt, 3) - 0.00024470 * pow(t_tt, 4)) * UniversalConstants::ARCSECONDS_2_RAD;
  }

  static double getMeanAnomalySun(double t_tt) {
    return (1287104.79305 + 129596581.0481 * t_tt - 0.5532 * pow(t_tt, 2) +
           0.000136 * pow(t_tt, 3) - 0.00001149 * pow(t_tt, 4)) * UniversalConstants::ARCSECONDS_2_RAD;
  }

  static double getMeanArgLatMoon(double t_tt) {
    return (335779.526232 + 1739527262.8478 * t_tt - 12.7512 * pow(t_tt, 2) -
           0.001037 * pow(t_tt, 3) + 0.00000417 * pow(t_tt, 4)) * UniversalConstants::ARCSECONDS_2_RAD;
  }

  static double getMeanElongationSun(double t_tt) {
    return (1072260.70369 + 1602961601.2090 * t_tt - 6.3706 * pow(t_tt, 2) +
           0.006593 * pow(t_tt, 3) - 0.00003169 * pow(t_tt, 4)) * UniversalConstants::ARCSECONDS_2_RAD;
  }

  static double getNodeMoon(double t_tt) {
    return (450160.398036 - 6962890.5431 * t_tt + 7.4722 * pow(t_tt, 2) +
           0.007702 * pow(t_tt, 3) - 0.00005939 * pow(t_tt, 4)) * UniversalConstants::ARCSECONDS_2_RAD;
  }

  static double getMeanLongMercury(double t_tt) {
    return 4.402608842 + 2608.7903141574 * t_tt;
  }

  static double getMeanLongVenus(double t_tt) {
    return 3.176146697 + 1021.3285546211 * t_tt;
  }

  static double getMeanLongEarth(double t_tt) {
    return 1.753470314 + 628.3075849991 * t_tt;
  }

  static double getMeanLongMars(double t_tt) {
    return 6.203480913 + 334.0612426700 * t_tt;
  }

  static double getMeanLongJupiter(double t_tt) {
    return 0.599546497 + 52.9690962641 * t_tt;
  }

  static double getMeanLongSaturn(double t_tt) {
    return 0.874016757 + 21.3299104960 * t_tt;
  }

  static double getMeanLongUranus(double t_tt) {
    return 5.481293872 + 7.478159856 * t_tt;
  }

  static double getMeanLongNeptune(double t_tt) {
    return 5.311886287 + 3.813303563 * t_tt;
  }

  static double getLongPrecession(double t_tt) {
    return 0.02438175 * t_tt + 0.00000538691 * t_tt * t_tt;
  }

  static std::vector<std::function<double(double t_tt)>> getAllFundamentalArguments(double t_tt) {
    getMeanAnomalyMoon(t_tt), getMeanAnomalySun(t_tt), getMeanArgLatMoon(t_tt), getMeanElongationSun(t_tt), getNodeMoon(t_tt),
    getMeanLongMercury(t_tt), getMeanLongVenus(t_tt), getMeanLongEarth(t_tt), getMeanLongMars(t_tt), getMeanLongJupiter(t_tt),
    getMeanLongSaturn(t_tt), getMeanLongUranus(t_tt), getMeanLongNeptune(t_tt), getLongPrecession(t_tt);
  }



};
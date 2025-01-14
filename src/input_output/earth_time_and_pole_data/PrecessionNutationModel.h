//
// Created by Douglas Garza on 10/14/24.
//

#pragma once

#include <cstdint>

class PrecessionNutationTable {
public:

  // Structure to hold each row of the non-polynomial part table
  struct NonPolynomialTerm {
    double C_s;
    double C_c;
    int8_t meanAnomMoon, meanAnomSun, argLatMoon, meanElongSun, nodeMoon;
    int8_t meanLong_Me, meanLong_Ve, meanLong_E, meanLong_Ma, meanLong_Ju, meanLong_Sa, meanLong_U, meanLong_Ne, longPrec;
  };


};

//
// Created by Douglas Garza on 10/15/24.
//

#pragma once

#include <array>
#include "PrecessionNutationModel.h"


class PrecessionNutation_sData : public PrecessionNutationTable {
public:

  // Polynomial coefficients for s + XY/2 (unit: microarcsecond)
  static constexpr std::array<double, 6> polynomialCoeffs = {
      94.0 * UniversalConstants::ARCSECONDS_2_RAD,
      3808.65 * UniversalConstants::ARCSECONDS_2_RAD,
      -122.68 * UniversalConstants::ARCSECONDS_2_RAD,
      -72574.11 * UniversalConstants::ARCSECONDS_2_RAD,
      27.98 * UniversalConstants::ARCSECONDS_2_RAD,
      15.62 * UniversalConstants::ARCSECONDS_2_RAD
  };

  // Non-polynomial part for j = 0
  static constexpr std::array<NonPolynomialTerm, 33> nonPolynomialJ0 = {{
      {-2640.73, 0.39, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-63.53, 0.02, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-11.75, -0.01, 0, 0, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-11.21, -0.01, 0, 0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {4.57, 0.00, 0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-2.02, 0.00, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-1.98, 0.00, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.72, 0.00, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.41, 0.01, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.26, 0.01, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.63, 0.00, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.63, 0.00, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.46, 0.00, 0, 1, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.45, 0.00, 0, 1, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.36, 0.00, 0, 0, 4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.24, 0.12, 0, 0, 1, -1, 1, 0, -8, 12, 0, 0, 0, 0, 0, 0},
      {-0.32, 0.00, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.28, 0.00, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.27, 0.00, 1, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.26, 0.00, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.21, 0.00, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.19, 0.00, 0, 1, -2, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.18, 0.00, 0, 1, -2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.10, -0.05, 0, 0, 0, 0, 0, 0, 8, -13, 0, 0, 0, 0, 0, -1},
      {-0.15, 0.00, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.14, 0.00, 2, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.14, 0.00, 0, 1, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.14, 0.00, 1, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.14, 0.00, 1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.13, 0.00, 0, 0, 4, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.11, 0.00, 0, 0, 2, -2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.11, 0.00, 1, 0, -2, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.11, 0.00, 1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  }};

  // Non-polynomial part for j = 1
  static constexpr std::array<NonPolynomialTerm, 3> nonPolynomialJ1 = {{
      {-0.07, 3.57, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.73, -0.03, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.00, 0.48, 0, 0, 2, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  }};

  // Non-polynomial part for j = 2
  static constexpr std::array<NonPolynomialTerm, 25> nonPolynomialJ2 = {{
      {743.52, -0.17, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {56.91, 0.06, 0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {9.84, -0.01, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-8.85, 0.01, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-6.38, -0.05, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-3.07, 0.00, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {2.23, 0.00, 0, 1, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.67, 0.00, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.30, 0.00, 1, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.93, 0.00, 0, 1, -2, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.68, 0.00, 1, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.55, 0.00, 0, 0, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.53, 0.00, 1, 0, -2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.27, 0.00, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.27, 0.00, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.26, 0.00, 1, 0, -2, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.25, 0.00, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.22, 0.00, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.21, 0.00, 2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.20, 0.00, 2, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.17, 0.00, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.13, 0.00, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.13, 0.00, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.12, 0.00, 1, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.11, 0.00, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  }};

  // Non-polynomial part for j = 3
  static constexpr std::array<NonPolynomialTerm, 4> nonPolynomialJ3 = {{
      {0.30, -23.42, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.03, -1.46, 0, 0, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {-0.01, -0.25, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0.00, 0.23, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  }};

  // Non-polynomial part for j = 4
  static constexpr std::array<NonPolynomialTerm, 1> nonPolynomialJ4 = {{
    {-0.26, -0.01, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  }};

  // auto nonPolynomialCoeffs = std::make_tuple(nonPolynomialJ0, nonPolynomialJ1, nonPolynomialJ2, nonPolynomialJ3, nonPolynomialJ4);

};
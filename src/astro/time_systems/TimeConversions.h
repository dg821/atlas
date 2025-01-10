//
// Created by Douglas Garza on 9/23/24.
//

#pragma once

#include <unordered_set>
#include <cmath>
#include "../../input_ouput/earth_time_and_pole_data/ExtractTimeAndPoleData.h"

namespace timeConversions {
    struct GregorianDate {
        int year{}, month{}, day{}, hour{}, minute{};
        double second{};
    };

    const std::unordered_set<int> leap_second_mjds = {
        41498, 41682, 42047, 42412, 42777, 43143, 43508, 43873, 44238,
        44785, 45150, 45515, 46246, 47160, 47891, 48256, 48803, 49168,
        49533, 50082, 50629, 51178, 53735, 54831, 56108, 57203, 57753
    };

    double date_to_jd(double year, double month, double day, double hour, double minute, double second);
    double jd_to_mjd(double jd);
    double jd_to_timej2k(double jd);
    double date_to_timej2k( double year, double month, double day, double hour, double minute, double second);

    GregorianDate jd_to_date(double jd);
    double timej2k_to_jd(double timej2k);
    GregorianDate timej2k_to_date(double timej2k);

    double jd_to_besselian_year(double jd_tt);
    double get_delta_ut2(double besselianYear);
    double getAuxiliaryTerrestrialTime(double timej2k);
    double getJulianCenturiesOfTT(double timej2k);

    double T_to_JD(double T);
    double secondsToT(double timej2k)
}
//
// Created by Douglas Garza on 9/23/24.
//

#pragma once

#include "../../input_output/earth_time_and_pole_data/ExtractTimeAndPoleData.h"
#include <map>
#include <cmath>

namespace timeConversions {
    struct GregorianDate {
        int year{}, month{}, day{}, hour{}, minute{};
        double second{};

        bool operator<(const GregorianDate& other) const {
            if (year != other.year) return year < other.year;
            if (month != other.month) return month < other.month;
            if (day != other.day) return day < other.day;
            if (hour != other.hour) return hour < other.hour;
            if (minute != other.minute) return minute < other.minute;
            return second < other.second;
        }
    };

    double date_to_jd(int yearIn, int monthIn, int dayIn, int hourIn, int minuteIn, double secondIn);
    double jd_to_mjd(double jd);
    double jd_to_timej2k(double jd);
    double date_to_timej2k( int year, int month, int day, int hour, int minute, double second);

    GregorianDate jd_to_date(double jd);
    double timej2k_to_jd(double timej2k);
    GregorianDate timej2k_to_date(double timej2k);

    double jd_to_besselian_year(double jd);
    double getJulianCenturiesOfTT(double JD);
    double getDeltaUT1(double timej2k);

    GregorianDate getDateTT2UTC(GregorianDate TT);
    GregorianDate getDateUTC2TT(GregorianDate UTC);

    double T_to_JD(double T);
    double secondsToT(double timej2k);
    int getLeapSeconds(const GregorianDate& date);

    double getDaysSinceJ2000(double timej2k);
}
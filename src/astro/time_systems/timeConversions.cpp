//
// Created by Douglas Garza on 9/23/24.
//

#include "timeConversions.h"

#include <math/UniversalConstants.h>
namespace timeConversions {

    const std::map<GregorianDate, int> leap_seconds = {
        {{1972, 1, 1, 0, 0, 0}, 10},
        {{1972, 7, 1, 0, 0, 0}, 11},
        {{1973, 1, 1, 0, 0, 0}, 12},
        {{1974, 1, 1, 0, 0, 0}, 13},
        {{1975, 1, 1, 0, 0, 0}, 14},
        {{1976, 1, 1, 0, 0, 0}, 15},
        {{1977, 1, 1, 0, 0, 0}, 16},
        {{1978, 1, 1, 0, 0, 0}, 17},
        {{1979, 1, 1, 0, 0, 0}, 18},
        {{1980, 1, 1, 0, 0, 0}, 19},
        {{1981, 7, 1, 0, 0, 0}, 20},
        {{1982, 7, 1, 0, 0, 0}, 21},
        {{1983, 7, 1, 0, 0, 0}, 22},
        {{1985, 7, 1, 0, 0, 0}, 23},
        {{1988, 1, 1, 0, 0, 0}, 24},
        {{1990, 1, 1, 0, 0, 0}, 25},
        {{1991, 1, 1, 0, 0, 0}, 26},
        {{1992, 7, 1, 0, 0, 0}, 27},
        {{1993, 7, 1, 0, 0, 0}, 28},
        {{1994, 7, 1, 0, 0, 0}, 29},
        {{1996, 1, 1, 0, 0, 0}, 30},
        {{1997, 7, 1, 0, 0, 0}, 31},
        {{1999, 1, 1, 0, 0, 0}, 32},
        {{2006, 1, 1, 0, 0, 0}, 33},
        {{2009, 1, 1, 0, 0, 0}, 34},
        {{2012, 7, 1, 0, 0, 0}, 35},
        {{2015, 7, 1, 0, 0, 0}, 36},
        {{2017, 1, 1, 0, 0, 0}, 37}
    };

    GregorianDate getDateUTC2TT(GregorianDate UTC) {
        GregorianDate TT = UTC;
        int deltaAT = getLeapSeconds(UTC);

        TT.second = UTC.second + deltaAT + UniversalConstants::TAI_2_TT;

        return TT;
    }

    GregorianDate getDateTT2UTC(GregorianDate TT) {

        GregorianDate UTC = TT;
        int deltaAT = getLeapSeconds(TT);

        UTC.second = TT.second - deltaAT - UniversalConstants::TAI_2_TT;

        return UTC;
    }

    double date_to_jd(int yearIn, int monthIn, const int dayIn, const int hourIn, const int minuteIn, const double secondIn) {
        // input date in UTC, output jd in tt
        GregorianDate utcDate = {yearIn, monthIn, dayIn, hourIn, minuteIn, secondIn};
        GregorianDate ttDate = getDateUTC2TT(utcDate);

        int year = ttDate.year; int month = ttDate.month; int day = ttDate.day; int hour = ttDate.hour; int minute = ttDate.minute; double second = ttDate.second;

        // For epoch 1900 to 2100
        if (month == 1 || month == 2) {
            year -= 1;
            month += 12;
        }

        const double B = 2.0 - std::trunc(year / 100.0) + std::trunc(std::trunc(year / 100.0) / 4.0);
        double C = ((second / 60.0 + minute) / 60.0 + hour) / 24.0;

        const double jd_day = std::trunc(365.25 * (year + 4716.0)) + std::trunc(30.6001 * (month + 1)) + day + B - 1524.5;
        double jd_day_fraction = C;

        double jd = jd_day + jd_day_fraction;

        GregorianDate dateComp {year, month, day, 0, 0, 0};
        if (leap_seconds.contains(dateComp)) {
            C = ((second / 61.0 + minute) / 60.0 + hour) / 24.0;
            jd_day_fraction = C;

            jd = jd_day + jd_day_fraction;
        }

        return jd;
    }


    double jd_to_mjd(const double jd) {
        const double mjd = jd - 2400000.5;
        return mjd;
    }

    double jd_to_timej2k(const double jd) {

        constexpr double j2000_epoch = 2451545.0;           // TT
        constexpr double seconds_per_day = 86400.0;

        const double timej2k = (jd - j2000_epoch) * seconds_per_day;              // seconds
        return timej2k;
    }


    double date_to_timej2k(int year, int month, const int day, const int hour, const int minute, const double second) {
        // input in UTC, output in TT
        const double jd = date_to_jd(year, month, day, hour, minute, second);           // input in UTC, output in TT
        const double timej2k = jd_to_timej2k(jd);

        return timej2k;
    };



   GregorianDate jd_to_date(const double jd) {
        // input in TT, output in UTC
       GregorianDate date;

        // Add 0.5 to JD and take the integer part
        const int Z = static_cast<int>(std::floor(jd + 0.5));

        // Calculate the fractional part of the day
        const double F = (jd + 0.5) - Z;

        double second = F * 86400.0;
        const int hour = static_cast<int>(second / 3600);
        second -= hour * 3600;
        const int minute = static_cast<int>(second / 60);
        second -= minute * 60;

        date.hour = hour;
        date.minute = minute;
        date.second = second;

        // Calculate alpha
        const int alpha = static_cast<int>((Z - 1867216.25) / 36524.25);

        // Calculate A
        const int A = Z + 1 + alpha - static_cast<int>(alpha / 4);

        // Calculate B, C, D
        const int B = A + 1524;
        const int C = static_cast<int>((B - 122.1) / 365.25);
        const int D = static_cast<int>(365.25 * C);

        // Calculate E, month
        const int E = static_cast<int>((B - D) / 30.6001);
        date.month = (E < 14) ? E - 1 : E - 13;

        // Calculate year, day
        date.year = (date.month > 2) ? C - 4716 : C - 4715;
        date.day = B - D - static_cast<int>(30.6001 * E);

        GregorianDate outDate = getDateTT2UTC(date);  // input in TT, output in UTC

        return outDate;
    }

   double timej2k_to_jd(const double timej2k) {


        constexpr double j2000_epoch = 2451545.0;       // TT
        constexpr double seconds_per_day = 86400.0;

        const double jd = timej2k / seconds_per_day + j2000_epoch;

        return jd;
   }


    GregorianDate timej2k_to_date(const double timej2k) {
       // input in TT, output in UTC
        const double jd = timej2k_to_jd(timej2k);

        const GregorianDate date = jd_to_date(jd);  // input in TT, output in UTC

        return date;
    }

    double jd_to_besselian_year(double jd) {
        double mjd = jd_to_mjd(jd);
        double B = 2000.0 + (mjd - 51544.03) / 365.242198781;
        return B;
    }

    double getJulianCenturiesOfTT(double jd) {
        double T = (jd - 2451545.0) / 36525.0;
        return T;
    }

    // convert Julian centuries to Julian Date
    double T_to_JD(double T) {
        return T * 36525.0 + 2451545.0;
    }

    // convert j2k seconds to Julian centuries
    double secondsToT(double timej2k) {
        return timej2k / (86400.0 * 36525.0);
    }

    int getLeapSeconds(const GregorianDate& date) {
       int leap_seconds_count = 10;  // initialize to first value (1972)

       for (const auto& [leap_date, count] : leap_seconds) {
           if (!(date < leap_date)) {
               leap_seconds_count = count;
           } else {
               break;
           }
       }

       return leap_seconds_count;
   }


    double getDaysSinceJ2000(double timej2k) {

       double jd = timej2k_to_jd(timej2k);

       // subtract J2000 epoch (2451545.0)
       return jd - 2451545.0;
    }

}


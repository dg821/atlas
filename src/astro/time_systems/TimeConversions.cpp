//
// Created by Douglas Garza on 9/23/24.
//

#include "TimeConversions.h"

#include <math/UniversalConstants.h>

double date_to_jd(double year, double month, const double day, const double hour, const double minute, const double second) {

    if (month == 1 || month == 2) {
        year -= 1;
        month += 12;
    }

    const double B = 2.0 - std::trunc(year / 100.0) + std::trunc(std::trunc(year / 100.0) / 4.0);
    double C = ((second / 60.0 + minute) / 60.0 + hour) / 24.0;

    const double jd_day = std::trunc(365.25 * (year + 4716.0)) + std::trunc(30.6001 * (month + 1)) + day + B - 1524.5;
    double jd_day_fraction = C;

    double jd = jd_day + jd_day_fraction;

    double mjd = timeConversions::jd_to_mjd(jd);

    if (timeConversions::leap_second_mjds.contains(std::floor(mjd))) {
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
    constexpr double j2000_epoch = 2451545.0;
    constexpr double seconds_per_day = 86400.0;

    const double timej2k = (jd - j2000_epoch) * seconds_per_day;              // seconds
    return timej2k;
}


double date_to_timej2k(const double year, const double month, const double day, const double hour, const double minute, const double second) {

    const double jd = date_to_jd(year, month, day, hour, minute, second);
    const double timej2k = jd_to_timej2k(jd);

    return timej2k;
};



timeConversions::GregorianDate jd_to_date(const double jd) {
    timeConversions::GregorianDate date;

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

    return date;
}

double timej2k_to_jd(const double timej2k) {
    constexpr double j2000_epoch = 2451545.0;
    constexpr double seconds_per_day = 86400.0;

    const double jd = timej2k / seconds_per_day + j2000_epoch;

    return jd;
}


timeConversions::GregorianDate timej2k_to_date(const double timej2k) {
    const double jd = timej2k_to_jd(timej2k);

    const timeConversions::GregorianDate date = jd_to_date(jd);

    return date;
}

double jd_to_besselian_year(double jd_tt) {
    double B = 1900.0 + (jd_tt - 241502.31352) / 365.242198781;
    return B;
}

double get_delta_ut2(double besselianYear) {
    double delta_ut2ut1 = 0.022 * std::sin(2.0 * M_PI * besselianYear) - 0.012 * std::cos(2.0 * M_PI * besselianYear) - 0.006 * std::sin(4.0 * M_PI * besselianYear) + 0.007 * std::cos(4.0 * M_PI * besselianYear);
    return delta_ut2ut1;
}

double get_delta_ut1(double timej2k) {
    double timej2k_tt = timej2k + UniversalConstants::UTC_2_TAI + UniversalConstants::TAI_2_TT;
    double jd_tt = timej2k_to_jd(timej2k_tt);
    double mjd = jd_to_mjd(timej2k_to_jd(timej2k));

    double besselianYear = jd_to_besselian_year(jd_tt);
    double delta_ut2 = get_delta_ut2(besselianYear);
    double delta_ut1 = 0.5309 - 0.00123(mjd - 57808) - delta_ut2;

    return delta_ut1;
}


std::tuple<double, timeConversions::GregorianDate> getAuxiliaryTerrestrialTime(const double timej2k) {
    // UTC, UT1, GPS, TAI, TT

    timeConversions::GregorianDate UTC = timej2k_to_date(timej2k);

    TimeAndPoleData tpdata;
    timeConversions::GregorianDate dateToday = tpdata.getCurrentUTCTime();

    double timej2kToday = date_to_timej2k(dateToday);

    double deltaUT1;
    double deltaAT;
    // if (timej2kToday > timej2k) {
    //     // get UT1 from database
    //     // WCAFTL
    //     TimeAndPoleData tpdata;
    //     std::string deltaUT1Data = tpdata.fetchData(tpdata.getDeltaUt1andPoleUrl());
    //     deltaUT1 = tpdata.parseDeltaUT1(deltaUT1Data);
    //
    //     std::string deltaATData = tpdata.fetchData(tpdata.getDeltaATUrl());
    //     deltaAT = tpdata.parseDeltaAT(deltaATData);
    //
    // } else {
        // get UT1 from calculation
        deltaUT1 = get_delta_ut1(timej2k);
        deltaAT = UniversalConstants::UTC_2_TAI;
    // }

    timeConversions::GregorianDate UT1 = UTC;
    timeConversions::GregorianDate TAI = UTC;
    timeConversions::GregorianDate TT = UTC;
    timeConversions::GregorianDate GPS = UTC;
    UT1.second = UTC.second + deltaUT1;
    TAI.second = UTC.second + deltaAT;
    GPS.second = UTC.second + UniversalConstants::UTC_2_GPS;
    TT.second = TAI.second + UniversalConstants::TAI_2_TT;

    double jd_tt = date_to_jd(TT.year, TT.month, TT.day, TT.hour, TT.minute, TT.second);

    double t_tt = (jd_tt - 2451545.0) / 36525.0;

    return std::make_tuple(t_tt, UT1);
}

double getJulianCenturiesOfTT(double timej2k) {
    std::tuple<double, timeConversions::GregorianDate> auxTuple = getAuxiliaryTerrestrialTime(timej2k);

    double t_tt = std::get<0>(auxTuple);
    return t_tt;
}

// convert Julian centuries to Julian Date
double T_to_JD(double T) {
    return T * 36525.0 + 2451545.0;
}

// convert j2k seconds to Julian centuries
double secondsToT(double timej2k) {
    return timej2k / (86400.0 * 36525.0);
}
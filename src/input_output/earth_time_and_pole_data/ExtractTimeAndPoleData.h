//
// Created by Douglas Garza on 10/10/24.
//

#pragma once
#include <iostream>
#include <string>
#include <curl/curl.h>
#include <regex>
#include <vector>
#include <sstream>
#include <cmath>
#include "../../astro/time_systems/TimeConversions.h"


// WCAFTL: needs work

class TimeAndPoleData {

    public:

    struct PoleData {
        static double mjd;
        static double x_pole;
        static double y_pole;
    };

        size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output);
        std::string fetchData(const std::string& url);
        double parseDeltaUT1(const std::string& data);
        int parseDeltaAT(const std::string& data);
        std::vector<PoleData> parseIERSPoleData(const std::string& data);
        timeConversions::GregorianDate getCurrentUTCTime();

        std::string getDeltaUt1andPoleUrl() const { return deltaUT1adnPoleUrl; }
        std::string getDeltaATUrl() const { return deltaATUrl; }
        std::string getPoleDataURL() const {return deltaUT1adnPoleUrl; }

    private:
        const std::string deltaUT1adnPoleUrl = "https://datacenter.iers.org/data/latestVersion/6_BULLETIN_A_V2013_016.txt";
        const std::string deltaATUrl = "https://maia.usno.navy.mil/ser7/tai-utc.dat";
};
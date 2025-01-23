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

class TimeAndPoleData {
public:
    struct PoleData {
        double mjd;
        double x_pole;
        double y_pole;
    };

    static double getDeltaUT1FromData(double timej2k);

private:
    static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output);
    static std::string fetchData(const std::string& url);
};
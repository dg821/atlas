//
// Created by Douglas Garza on 10/10/24.
//

#include "ExtractTimeAndPoleData.h"


// Callback function for cURL to write received data
size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

// Function to fetch data from a URL
std::string fetchData(const std::string& url) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);

        if (res != CURLE_OK) {
            std::cerr << "cURL error: " << curl_easy_strerror(res) << std::endl;
            return "";
        }
    }

    return readBuffer;
}

// Function to parse Delta UT1 from IERS data
double parseDeltaUT1(const std::string& data) {
    std::regex pattern(R"(\d{4}\s+\d+\s+MJD\s+\d+\.\d+\s+(\S+))");
    std::smatch matches;
    if (std::regex_search(data, matches, pattern) && matches.size() > 1) {
        return std::stod(matches[1].str());
    }
    return 0.0;
}

// Function to parse Delta AT from USNO data
int parseDeltaAT(const std::string& data) {
    std::regex pattern(R"((\d{4} \w+ \d+) = (\d+))");
    std::smatch matches;
    if (std::regex_search(data, matches, pattern) && matches.size() > 2) {
        return std::stoi(matches[2].str());
    }
    return 0;
}

timeConversions::GregorianDate getCurrentUTCTime() {
    // Get current time point
    auto now = std::chrono::system_clock::now();

    // Get duration since epoch
    auto duration = now.time_since_epoch();

    // Get seconds and fractional seconds
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto fractional_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(duration) -
                              std::chrono::duration_cast<std::chrono::duration<double>>(seconds);

    // Convert to std::time_t
    auto timer = std::chrono::system_clock::to_time_t(now);

    // Convert to UTC
    std::tm utc_tm = *std::gmtime(&timer);

    // Populate GregorianDate struct
    timeConversions::GregorianDate date;
    date.year = utc_tm.tm_year + 1900; // tm_year is years since 1900
    date.month = utc_tm.tm_mon + 1;    // tm_mon is 0-11
    date.day = utc_tm.tm_mday;
    date.hour = utc_tm.tm_hour;
    date.minute = utc_tm.tm_min;
    date.second = utc_tm.tm_sec + fractional_seconds.count();

    return date;
}


std::vector<TimeAndPoleData::PoleData> parseIERSPoleData(const std::string& data) {
    std::vector<TimeAndPoleData::PoleData> result;
    std::istringstream dataStream(data);
    std::string line;

    // Regular expression to match the data lines
    std::regex dataRegex(R"(\s*(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+))");

    // Skip header lines
    while (std::getline(dataStream, line)) {
        if (line.find("x pol") != std::string::npos) {
            break;
        }
    }

    // Parse data lines
    while (std::getline(dataStream, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, dataRegex)) {
            TimeAndPoleData::PoleData pd;
            int year = std::stoi(matches[1]);
            int month = std::stoi(matches[2]);
            int day = std::stoi(matches[3]);
            pd.mjd = 367 * year - 7 * (year + (month + 9) / 12) / 4 + 275 * month / 9 + day - 678987;
            pd.x_pole = std::stod(matches[4]);
            pd.y_pole = std::stod(matches[5]);
            result.push_back(pd);
        }
    }

    return result;
}

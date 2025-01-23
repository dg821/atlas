#include "ExtractTimeAndPoleData.h"

size_t TimeAndPoleData::WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t totalSize = size * nmemb;
    output->append(static_cast<char*>(contents), totalSize);
    return totalSize;
}

std::string TimeAndPoleData::fetchData(const std::string& url) {
    CURL* curl = curl_easy_init();
    std::string readBuffer;

    if (!curl) {
        throw std::runtime_error("Failed to initialize CURL");
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);

    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        throw std::runtime_error(std::string("CURL error: ") + curl_easy_strerror(res));
    }

    return readBuffer;
}

double TimeAndPoleData::getDeltaUT1FromData(double timej2k) {
    // Convert J2000 seconds to MJD
    double daysSinceJ2000 = timej2k / (24.0 * 3600.0);
    double targetMJD = 51544.5 + daysSinceJ2000;

    // Fetch the latest IERS data
    try {
        std::string iersUrl = "https://datacenter.iers.org/data/latestVersion/finals2000A.all";
        std::string iersData = fetchData(iersUrl);

        if (iersData.empty()) {
            throw std::runtime_error("Empty response from IERS server");
        }

        // regular expression to match MJD and UT1-UTC values
        std::regex pattern(R"(\d{4}\s+\d+\s+MJD\s+(\d+\.\d+)\s+(\S+))");
        std::smatch matches;
        std::string::const_iterator searchStart(iersData.cbegin());

        double closestMJD = 0.0;
        double closestDeltaUT1 = 0.0;
        double smallestDiff = std::numeric_limits<double>::max();

        // search through all matches to find the closest MJD
        while (std::regex_search(searchStart, iersData.cend(), matches, pattern)) {
            double mjd = std::stod(matches[1].str());
            double deltaUT1 = std::stod(matches[2].str());

            double diff = std::abs(mjd - targetMJD);
            if (diff < smallestDiff) {
                smallestDiff = diff;
                closestMJD = mjd;
                closestDeltaUT1 = deltaUT1;
            }

            searchStart = matches.suffix().first;
        }

        if (closestMJD == 0.0) {
            throw std::runtime_error("No valid UT1-UTC data found");
        }

        return closestDeltaUT1;

    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error in getDeltaUT1FromData: ") + e.what());
    }
}
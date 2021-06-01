#ifdef __JSGEODA__
#include <emscripten/bind.h>
#endif

#include <iostream>
#include <rapidjson/stringbuffer.h>
#include "../libgeoda_src/weights/GalWeight.h"
#include "../libgeoda_src/GenUtils.h"
#include "geojson.h"
#include "jsgeoda.h"

extern std::map<std::string, GdaGeojson*> geojson_maps;

std::vector<double> natural_breaks(int k, const std::vector<double>& data,
                                    const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::NaturalBreaks(k, data, undefs_b);
}

std::vector<double> quantile_breaks(int k, const std::vector<double>& data,
                                   const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::QuantileBreaks(k, data, undefs_b);
}

std::vector<double> percentile_breaks(const std::vector<double>& data,
                                  const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::PercentileBreaks(data, undefs_b);
}

std::vector<double> stddev_breaks(const std::vector<double>& data,
                                  const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::StddevBreaks(data, undefs_b);
}

std::vector<double> hinge15_breaks(const std::vector<double>& data,
                                  const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::Hinge15Breaks(data, undefs_b);
}

std::vector<double> hinge30_breaks(const std::vector<double>& data,
                                   const std::vector<int>& undefs)
{
    std::vector<bool> undefs_b(undefs.size());
    for (size_t i=0; i<undefs.size(); ++i) {
        undefs_b[i] = undefs[i] == 0 ? false : true;
    }
    return GenUtils::Hinge30Breaks(data, undefs_b);
}
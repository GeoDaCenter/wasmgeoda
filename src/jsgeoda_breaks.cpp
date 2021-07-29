//
// Created by Xun Li on 2/1/21. <lixun910@gmail.com>
//

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

std::vector<double> excess_risk(const std::vector<double>& event_data,
                                   const std::vector<double>& base_data)
{
    int obs = event_data.size();
    std::vector<bool> undefs(obs, false);
    double* P = new double[obs];
    double* E = new double[obs];
    double* r = new double[obs];

    for (int i=0; i<obs; ++i) {
        P[i] = event_data[i];
        E[i] = base_data[i];
        r[i] = 0;
    }
    GdaAlgs::RateSmoother_ExcessRisk(obs, P, E, r, undefs);

    std::vector<double> result(obs);
    for (int i=0; i<obs; ++i) {
        result[i] = r[i];
    }

    return result;
}

std::vector<double> eb_risk(const std::vector<double>& event_data,
                                const std::vector<double>& base_data)
{
    int obs = event_data.size();
    std::vector<bool> undefs(obs, false);
    double* P = new double[obs];
    double* E = new double[obs];
    double* r = new double[obs];

    for (int i=0; i<obs; ++i) {
        P[i] = event_data[i];
        E[i] = base_data[i];
        r[i] = 0;
    }
    GdaAlgs::RateSmoother_EBS(obs, P, E, r, undefs);

    std::vector<double> result(obs);
    for (int i=0; i<obs; ++i) {
        result[i] = r[i];
    }

    delete[] P;
    delete[] E;
    delete[] r;

    return result;
}

std::vector<double> spatial_lag(const std::string map_uid,
                                const std::string& weight_uid,
                                const std::vector<double>& data, bool is_binary,
                                bool row_stand, bool inc_diag)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int obs = (int)data.size();
            std::vector<double> result(obs);
            std::vector<bool> undefs(obs, false);

            for (int i=0; i<obs; ++i) {
                std::vector<long> nbrs = w->GetNeighbors(i);
                std::vector<double> wvals = w->GetNeighborWeights(i);
                double lag = 0;
                if (is_binary || wvals.empty()) {
                    for (int j=0; j < nbrs.size(); ++j) {
                        if (nbrs[j] != i || inc_diag) {
                            lag += data[nbrs[j]];
                        }
                    }
                    if (nbrs.empty() == false && row_stand) {
                        lag = lag / nbrs.size();
                    }
                } else {
                    double sumW = 0;
                    for (int j=0; j < nbrs.size(); ++j) {
                        if (nbrs[j] != i || inc_diag) {
                            sumW += wvals[j];
                        }
                    }
                    if (sumW ==0) {
                        lag = 0;
                    } else {
                        for (int j = 0; j < nbrs.size(); ++j) {
                            if (nbrs[j] != i || inc_diag) {
                                lag += data[nbrs[j]] * wvals[j] / sumW;
                            }
                        }
                    }
                }
                result[i] = lag;
            }
            return result;
        }
    }
    return std::vector<double>();
}

std::vector<double> spatial_rate(const std::vector<double>& event_data,
                                 const std::vector<double>& base_data,
                                 const std::string map_uid,
                                 const std::string& weight_uid)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int obs = (int) event_data.size();
            std::vector<bool> undefs(obs, false);
            double* P = new double[obs];
            double* E = new double[obs];
            double* r = new double[obs];

            for (int i=0; i<obs; ++i) {
                P[i] = event_data[i];
                E[i] = base_data[i];
                r[i] = 0;
            }

            GdaAlgs::RateSmoother_SRS(obs, w, E, P, r, undefs);

            std::vector<double> result(obs);
            for (int i=0; i<obs; ++i) {
                result[i] = r[i];
            }

            delete[] P;
            delete[] E;
            delete[] r;

            return result;
        }
    }

    return std::vector<double>();
}

std::vector<double> spatial_eb(const std::vector<double>& event_data,
                                 const std::vector<double>& base_data,
                                 const std::string map_uid,
                                 const std::string& weight_uid)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int obs = (int) event_data.size();
            std::vector<bool> undefs(obs, false);
            double* P = new double[obs];
            double* E = new double[obs];
            double* r = new double[obs];

            for (int i=0; i<obs; ++i) {
                P[i] = event_data[i];
                E[i] = base_data[i];
                r[i] = 0;
            }

            GdaAlgs::RateSmoother_SEBS(obs, w, E, P, r, undefs);

            std::vector<double> result(obs);
            for (int i=0; i<obs; ++i) {
                result[i] = r[i];
            }

            delete[] P;
            delete[] E;
            delete[] r;

            return result;
        }
    }

    return std::vector<double>();
}
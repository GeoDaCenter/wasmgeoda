//
// Created by Xun Li on 2/1/21. <lixun910@gmail.com>
//

#include "../libgeoda_src/sa/UniLocalMoran.h"
#include "../libgeoda_src/sa/UniG.h"
#include "../libgeoda_src/gda_sa.h"
#include "../libgeoda_src/weights/GalWeight.h"
#include "../libgeoda_src/GenUtils.h"
#include "geojson.h"
#include "jsgeoda.h"

extern std::map<std::string, GdaGeojson*> geojson_maps;


void set_lisa_content(LISA* lisa, LisaResult& rst)
{
    rst.is_valid = true;
    rst.sig_local_vec = lisa->GetLocalSignificanceValues();
    rst.sig_cat_vec = lisa->GetSigCatIndicators();
    rst.cluster_vec = lisa->GetClusterIndicators();
    rst.lag_vec = lisa->GetSpatialLagValues();
    rst.lisa_vec = lisa->GetLISAValues();
    rst.nn_vec = lisa->GetNumNeighbors();
    rst.labels = lisa->GetLabels();
    rst.colors = lisa->GetColors();
}

LisaResult local_moran(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(vals.size(), false);
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_localmoran(w, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                        permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_g(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                   const std::vector<int>& undefs, double significance_cutoff, int permutations,
                   const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(undefs.size());
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_localg(w, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_gstar(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                   const std::vector<int>& undefs, double significance_cutoff, int permutations,
                   const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(undefs.size());
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_localgstar(w, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_geary(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(undefs.size());
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_localgeary(w, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_joincount(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(undefs.size());
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_localjoincount(w, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                            permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult quantile_lisa(const std::string map_uid, const std::string weight_uid, int k, int quantile,
                         const std::vector<double>& vals, const std::vector<int>& undefs, double significance_cutoff,
                         int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<bool> undefs_b(undefs.size());
            for (size_t i=0; i<undefs.size(); ++i) {
                undefs_b[i] = undefs[i] == 0 ? false : true;
            }
            LISA* lisa = gda_quantilelisa(w, k, quantile, vals, undefs_b, significance_cutoff, nCPUs, permutations,
                                            permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_moran_eb(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<double>& base_values, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            LISA* lisa = gda_localmoran_eb(w, vals, base_values, significance_cutoff, nCPUs, permutations,
                                        permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_multigeary(const std::string map_uid, const std::string weight_uid,
                            const std::vector<std::vector<double> > &data,
                            const std::vector<std::vector<int> > &undefs, double significance_cutoff,
                            int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<std::vector<bool> > undefs_b(data.size());
            for (size_t i=0; i<data.size(); ++i) {
                undefs_b[i].resize(data[i].size(), false);
                for (size_t j=0; j < data[i].size(); ++j) {
                    undefs_b[i][j] = undefs[i][j] == 0 ? false : true;
                }
            }
            LISA* lisa = gda_localmultigeary(w, data, undefs_b, significance_cutoff, nCPUs, permutations,
                                        permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_multijoincount(const std::string map_uid, const std::string weight_uid,
                            const std::vector<std::vector<double> > &data,
                            const std::vector<std::vector<int> > &undefs, double significance_cutoff,
                            int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<std::vector<bool> > undefs_b(data.size());
            for (size_t i=0; i<data.size(); ++i) {
                undefs_b[i].resize(data[i].size(), false);
                for (size_t j=0; j < data[i].size(); ++j) {
                    undefs_b[i][j] = undefs[i][j] == 0 ? false : true;
                }
            }
            LISA* lisa = gda_localmultijoincount(w, data, undefs_b, significance_cutoff, nCPUs, permutations,
                                             permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult multi_quantile_lisa(const std::string map_uid, const std::string weight_uid,
                               const std::vector<int> &k_s, const std::vector<int> &quantile_s,
                               const std::vector<std::vector<double> > &data,
                               const std::vector<std::vector<int> > &undefs,
                               double significance_cutoff, int permutations,
                               const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            std::vector<std::vector<bool> > undefs_b(data.size());
            for (size_t i=0; i<data.size(); ++i) {
                undefs_b[i].resize(data[i].size(), false);
                for (size_t j=0; j < data[i].size(); ++j) {
                    undefs_b[i][j] = undefs[i][j] == 0 ? false : true;
                }
            }
            LISA* lisa = gda_multiquantilelisa(w, k_s, quantile_s, data, undefs_b, significance_cutoff, nCPUs, permutations,
                                                 permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

std::vector<std::vector<double> > neighbor_match_test(const std::string map_uid,
                                                      int knn, double power, bool is_inverse, bool is_arc,
                                                      bool is_mile, const std::vector<std::vector<double> >& data,
                                                      const std::string& scale_method, const std::string& dist_type)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::vector<std::vector<double> > r = gda_neighbor_match_test(json_map, knn, power, is_inverse, is_arc, is_mile,
                                                                    data, scale_method, dist_type);
        return r;
    }
    return std::vector<std::vector<double> >();
}
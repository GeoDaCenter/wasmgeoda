//
// Created by Xun Li on 6/3/21. <lixun910@gmail.com>
//

#include "../libgeoda_src/gda_clustering.h"
#include "../libgeoda_src/GenUtils.h"
#include "geojson.h"
#include "jsgeoda.h"

extern std::map<std::string, GdaGeojson*> geojson_maps;

ClusteringResult redcap(const std::string map_uid, const std::string weight_uid, int k, const std::string &method,
                        const std::vector<std::vector<double> > &data,
                        const std::vector<double>& bound_vals, double min_bound,
                        const std::string& scale_method, const std::string &distance_method)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int nCPUs = 1;
            int seed = 123456789;// not used
            std::vector<std::vector<int> > cluster_ids = gda_redcap(k, w, data, scale_method, method, distance_method,
                                                                    bound_vals, min_bound, seed, nCPUs, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult schc(const std::string map_uid, const std::string weight_uid, int k, const std::string &method,
                        const std::vector<std::vector<double> > &data,
                        const std::vector<double>& bound_vals, double min_bound,
                        const std::string& scale_method, const std::string &distance_method)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<std::vector<int> > cluster_ids = gda_schc(k, w, data, scale_method, method, distance_method,
                                                                  bound_vals, min_bound, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult azp_greedy(const std::string map_uid, const std::string weight_uid, int k,
                            const std::vector<std::vector<double> > &data, int inits, const std::vector<int>& init_regions,
                            const std::string& scale_method, const std::string &distance_method,
                            const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                            const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                            int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }

            std::vector<std::vector<int> > cluster_ids = gda_azp_greedy(k, w, data, scale_method, inits,
                                                                        in_min_bounds, in_max_bounds, init_regions,
                                                                        distance_method, seed, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult azp_sa(const std::string map_uid, const std::string weight_uid, int k, double cooling_rate, int sa_maxit,
                        const std::vector<std::vector<double> > &data, int inits, const std::vector<int>& init_regions,
                        const std::string& scale_method, const std::string &distance_method,
                        const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                        const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                        int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }

            std::vector<std::vector<int> > cluster_ids = gda_azp_sa(k, w, data, scale_method, inits, cooling_rate, sa_maxit,
                                                                    in_min_bounds, in_max_bounds, init_regions,
                                                                    distance_method, seed, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult azp_tabu(const std::string map_uid, const std::string weight_uid, int k, int tabu_length, int conv_tabu,
                          const std::vector<std::vector<double> > &data, int inits, const std::vector<int>& init_regions,
                          const std::string& scale_method, const std::string &distance_method,
                          const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                          const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                          int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }

            std::vector<std::vector<int> > cluster_ids = gda_azp_tabu(k, w, data, scale_method, inits, tabu_length, conv_tabu,
                                                                      in_min_bounds, in_max_bounds, init_regions,
                                                                      distance_method, seed, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult maxp_greedy(const std::string map_uid, const std::string weight_uid,
                             const std::vector<std::vector<double> > &data, int iterations,
                             const std::string& scale_method, const std::string &distance_method,
                             const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                             const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                             int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }
            std::vector<int> init_regions;
            std::vector<std::vector<int> > cluster_ids = gda_maxp_greedy(w, data, scale_method, iterations,
                                                                      in_min_bounds, in_max_bounds, init_regions,
                                                                      distance_method, seed, 1, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult maxp_sa(const std::string map_uid, const std::string weight_uid,
                             const std::vector<std::vector<double> > &data, int iterations, double cooling_rate, int sa_maxit,
                             const std::string& scale_method, const std::string &distance_method,
                             const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                             const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                             int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }
            std::vector<int> init_regions;
            std::vector<std::vector<int> > cluster_ids = gda_maxp_sa(w, data, scale_method, iterations, cooling_rate, sa_maxit,
                                                                         in_min_bounds, in_max_bounds, init_regions,
                                                                         distance_method, seed, 1, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}

ClusteringResult maxp_tabu(const std::string map_uid, const std::string weight_uid,
                         const std::vector<std::vector<double> > &data, int iterations, int tabu_length, int conv_tabu,
                         const std::string& scale_method, const std::string &distance_method,
                         const std::vector<std::vector<double> >& min_bounds_values, const std::vector<double>& min_bounds,
                         const std::vector<std::vector<double> >& max_bounds_values, const std::vector<double>& max_bounds,
                         int seed)
{
    ClusteringResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            int seed = 123456789;// not used

            std::vector<std::pair<double, std::vector<double> > > in_min_bounds, in_max_bounds;

            for (int i=0; i<min_bounds.size(); ++i) {
                in_min_bounds.push_back(std::make_pair(min_bounds[i], min_bounds_values[i]));
            }
            for (int i=0; i<max_bounds.size(); ++i) {
                in_max_bounds.push_back(std::make_pair(max_bounds[i], max_bounds_values[i]));
            }
            std::vector<int> init_regions;
            std::vector<std::vector<int> > cluster_ids = gda_maxp_tabu(w, data, scale_method, iterations, tabu_length, conv_tabu,
                                                                     in_min_bounds, in_max_bounds, init_regions,
                                                                     distance_method, seed, 1, NULL);

            rst.is_valid = true;
            rst.between_ss = gda_betweensumofsquare(cluster_ids, data);
            rst.total_ss = gda_totalsumofsquare(data);
            rst.ratio = rst.between_ss / rst.total_ss;
            rst.within_ss = gda_withinsumofsquare(cluster_ids, data);
            rst.cluster_vec = GenUtils::flat_2dclusters(w->num_obs, cluster_ids);
        }
    }
    return rst;
}
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

void set_weights_content(GeoDaWeight* w, const string& map_uid, WeightsResult& rst)
{
    if (w) {
        rst.map_uid = map_uid;
        rst.is_valid = true;
        rst.weight_type = w->weight_type;
        rst.is_symmetric = w->is_symmetric;
        rst.max_nbrs = w->max_nbrs;
        rst.mean_nbrs =  w->mean_nbrs;
        rst.median_nbrs =  w->median_nbrs;
        rst.min_nbrs =  w->min_nbrs;
        rst.num_obs = w->num_obs;
        rst.sparsity = w->sparsity;
        rst.uid = w->uid;
    }
}

WeightsResult queen_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold)
{
	std::cout << "queen_weights(0)" << map_uid << std::endl;
    WeightsResult rst;
    rst.is_valid = false;
	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
	    std::cout << "enter queen_weights();" << std::endl;
		GeoDaWeight *w = json_map->CreateQueenWeights(order, include_lower_order, precision_threshold);
	    set_weights_content(w, map_uid, rst);
	}
	return rst;
}

WeightsResult rook_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold)
{
    std::cout << "rook_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateRookWeights(order, include_lower_order, precision_threshold);
        set_weights_content(w, map_uid, rst);
    }
    return rst;
}

WeightsResult knn_weights(std::string map_uid, int k, double power, bool is_inverse, bool is_arc, bool is_mile)
{
    std::cout << "knn_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateKnnWeights(k, power, is_inverse, is_arc, is_mile);
        set_weights_content(w, map_uid, rst);
    }
    return rst;
}

WeightsResult dist_weights(std::string map_uid, double dist_thres, double power, bool is_inverse, bool is_arc, bool is_mile)
{
    std::cout << "distance_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateDistanceWeights(dist_thres, power, is_inverse, is_arc, is_mile);
        set_weights_content(w, map_uid, rst);
    }
    return rst;
}

WeightsResult kernel_weights(std::string map_uid, int k, std::string kernel, bool adaptive_bandwidth,
        bool use_kernel_diagonals, double power, bool is_inverse, bool is_arc, bool is_mile)
{
    std::cout << "kernel_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelKnnWeights(k, kernel_str, adaptive_bandwidth, use_kernel_diagonals,
                power, is_inverse, is_arc, is_mile);
        set_weights_content(w, map_uid, rst);
    }
    return rst;
}

WeightsResult kernel_bandwidth_weights(std::string map_uid, double dist_thres, std::string kernel,
    bool use_kernel_diagonals, double power, bool is_inverse, bool is_arc, bool is_mile)
{
    std::cout << "kernel_bandwidth_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelWeights(dist_thres, kernel_str, use_kernel_diagonals, power, is_inverse,
                                                       is_arc, is_mile);
        set_weights_content(w, map_uid, rst);
    }
    return rst;
}

double get_min_dist_threshold(std::string map_uid, bool is_arc, bool is_mile)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        return json_map->GetMinDistanceThreshold(is_arc, is_mile);
    }
    return 0;
}

#ifdef __JSGEODA__
#include <emscripten/bind.h>
#endif
#include <iostream>
#include <sstream>
#include <vector>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>
#include <boost/algorithm/string.hpp>

#include "../libgeoda_src/sa/UniLocalMoran.h"
#include "../libgeoda_src/sa/UniG.h"
#include "../libgeoda_src/sa/UniGstar.h"
#include "../libgeoda_src/sa/UniGeary.h"
#include "../libgeoda_src/sa/UniJoinCount.h"
#include "../libgeoda_src/clustering/DorlingCartogram.h"
#include "../libgeoda_src/weights/GalWeight.h"
#include "../libgeoda_src/GenUtils.h"
#include "../libgeoda_src/gda_weights.h"
#include "../libgeoda_src/gda_sa.h"
#include "../libgeoda_src/gda_clustering.h"

#include "geojson.h"

extern "C" {
	void print_json(char* content);
    //void new_geojsonmap(const char* file_name, const char* addr, const size_t & len);
    void new_geojsonmap1(const char* file_name, uint8_t* data, size_t len);
}

auto geojson_maps = std::map<std::string, GdaGeojson*>(); 

void free_geojsonmap() 
{
	std::map<std::string, GdaGeojson*>::iterator it;
	for (it = geojson_maps.begin(); it != geojson_maps.end(); ++it) {
		delete it->second;
	}
	geojson_maps.clear();
}

/**
 * Create a geojson map in memory using the user drag-n-dropped 
 * *.geojson file via webpage
 * 
 * @param int addr The pointer address of the byte array
 * @param size_t len The length of the member of byte array
 * 
 * @return std::string Return a unique string represents as a key
 * 			to access the geojson map. Used as a reference to 
 * 			all jsgeoda APIs. e.g. queen_weights(map_uid)
 * 
 */
void new_geojsonmap(std::string file_name, const int & addr, const size_t & len) {
	//We get out pointer as a plain int from javascript
	//We use a reinterpret_cast to turn our plain int into a uint8_t pointer. After
	//which we can play with the data just like we would normally.
	char* data = reinterpret_cast<char*>(addr);
    //std::cout << "new_geojsonmap():" << file_name << len << "," << data[0] << std::endl;

	GdaGeojson *json_map = new GdaGeojson(file_name.c_str(), data);

	// store globally, has to be release by calling free_geojsonmap()
	geojson_maps[file_name] = json_map;
}

void new_geojsonmap1(const char* file_name, uint8_t* in, size_t len) {
    char* data = new char[len+1];
    for (size_t i=0; i<len;++i) {
        data[i] = in[i];
        if (i == len-1)
        std::cout << data[i];
    }
    data[len] = '\0';
    GdaGeojson *json_map = new GdaGeojson(file_name, data);
    geojson_maps[std::string(file_name)] = json_map;
}

int get_num_obs(std::string map_uid) {
	std::cout << "get_num_obs()" << map_uid << std::endl;
	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		return json_map->GetNumObs();
	}
	return 0;
}


int get_map_type(std::string map_uid) {
	std::cout << "get_map_type()" << map_uid << std::endl;
	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		return json_map->GetMapType();
	}
	return 0;
}

struct CCentroids {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> get_x() { return x;}
    std::vector<double> get_y() { return y;}
};

struct WeightsResult {
    int weight_type;
    int num_obs;
    bool is_symmetric;
    double sparsity;
    int max_nbrs;
    int min_nbrs;
    double median_nbrs;
    double mean_nbrs;
    std::string uid;

    int get_weight_type() {return weight_type;}
    int get_num_obs() {return num_obs;}
    bool get_is_symmetric() { return is_symmetric;}
    double get_sparsity() { return sparsity;}
    int get_max_nbrs() {return max_nbrs;}
    int get_min_nbrs() {return min_nbrs;}
    double get_median_nbrs() { return median_nbrs;}
    double get_mean_nbrs() { return mean_nbrs;}
    std::string get_uid() { return uid;}
};

void set_weights_content(GeoDaWeight* w, WeightsResult& rst)
{
    if (w) {
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

CCentroids get_centroids(std::string map_uid)
{
    CCentroids ct;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        const std::vector<gda::PointContents*>& pts = json_map->GetCentroids();
        size_t num_obs = json_map->GetNumObs();
        ct.x.resize(num_obs);
        ct.y.resize(num_obs);

        for (size_t i=0; i<num_obs; ++i) {
            ct.x[i] = pts[i]->x;
            ct.y[i] = pts[i]->y;
        }
    }
    return ct;
}

WeightsResult queen_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold)
{
	std::cout << "queen_weights()" << map_uid << std::endl;
    WeightsResult rst;

	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		GeoDaWeight *w = json_map->CreateQueenWeights(order, include_lower_order, precision_threshold);
	    set_weights_content(w, rst);
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
        set_weights_content(w, rst);
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
        set_weights_content(w, rst);
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
        set_weights_content(w, rst);
    }
    return rst;
}

WeightsResult kernel_weights(std::string map_uid, int k, std::string kernel, bool adaptive_bandwidth,
        bool use_kernel_diagonals, bool is_arc, bool is_mile)
{
    std::cout << "kernel_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelWeights(k, kernel_str, adaptive_bandwidth, use_kernel_diagonals, is_arc, is_mile);
        set_weights_content(w, rst);
    }
    return rst;
}

WeightsResult kernel_bandwidth_weights(std::string map_uid, double dist_thres, std::string kernel,
    bool use_kernel_diagonals, bool is_arc, bool is_mile)
{
    std::cout << "kernel_bandwidth_weights()" << map_uid << std::endl;
    WeightsResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelWeights(dist_thres, kernel_str, use_kernel_diagonals, is_arc, is_mile);
        set_weights_content(w, rst);
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

std::vector<double> get_numeric_col(std::string map_uid, std::string col_name) {
    std::cout << "get_numeric_col()" << map_uid << std::endl;
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        return json_map->GetNumericCol(col_name);
    }
    return std::vector<double>();
}

// represent lisa results
struct LisaResult {
    bool is_valid;
    std::vector<double> sig_local_vec;
    std::vector<int> sig_cat_vec;
    std::vector<int> cluster_vec;
    std::vector<double> lag_vec;
    std::vector<double> lisa_vec;
    std::vector<int> nn_vec;
    std::vector<std::string> labels;
    std::vector<std::string> colors;

    bool valid() { return is_valid; }
    std::vector<double> get_sig_local() { return sig_local_vec;}
    std::vector<int> get_sig_cat() { return sig_cat_vec;}
    std::vector<int>  get_cluster() { return cluster_vec;}
    std::vector<double>  get_lag() { return lag_vec;}
    std::vector<double>  get_lisa() { return lisa_vec;}
    std::vector<int>  get_nn() { return nn_vec;}
    std::vector<std::string>  get_labels() { return labels;}
    std::vector<std::string>  get_colors() { return colors;}
};

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

LisaResult local_moran(const std::string map_uid, const std::string weight_uid, std::string col_name, double significance_cutoff,
                       int nCPUs, int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<double> values = json_map->GetNumericCol(col_name);
            std::vector<bool> undefs;
            LISA* lisa = gda_localmoran(w, values, undefs, significance_cutoff, nCPUs, permutations,
                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_moran1(const std::string map_uid, const std::string weight_uid, std::vector<double> vals,
        std::vector<bool> undefs, double significance_cutoff, int nCPUs, int permutations, const std::string&
        permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            LISA* lisa = gda_localmoran(w, vals, undefs, significance_cutoff, nCPUs, permutations,
                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_g(const std::string map_uid, const std::string weight_uid, std::string col_name, double significance_cutoff,
                   int nCPUs, int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<double> values = json_map->GetNumericCol(col_name);
            std::vector<bool> undefs;
            LISA* lisa = gda_localg(w, values, undefs, significance_cutoff, nCPUs, permutations,
                                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_gstar(const std::string map_uid, const std::string weight_uid, std::string col_name, double significance_cutoff,
                       int nCPUs, int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<double> values = json_map->GetNumericCol(col_name);
            std::vector<bool> undefs;
            LISA* lisa = gda_localgstar(w, values, undefs, significance_cutoff, nCPUs, permutations,
                                        permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_gstar1(const std::string map_uid, const std::string weight_uid, std::vector<double> vals,
        std::vector<bool> undefs, double significance_cutoff, int nCPUs, int permutations, const std::string&
        permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            LISA* lisa = gda_localgstar(w, vals, undefs, significance_cutoff, nCPUs, permutations,
                                        permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_geary(const std::string map_uid, const std::string weight_uid, std::string col_name, double significance_cutoff,
                       int nCPUs, int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<double> values = json_map->GetNumericCol(col_name);
            std::vector<bool> undefs;
            LISA* lisa = gda_localgeary(w, values, undefs, significance_cutoff, nCPUs, permutations,
                    permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

LisaResult local_joincount(const std::string map_uid, const std::string weight_uid, std::string col_name, double significance_cutoff,
                           int nCPUs, int permutations, const std::string& permutation_method, int last_seed_used)
{
    LisaResult rst;
    rst.is_valid = false;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            std::vector<double> values = json_map->GetNumericCol(col_name);
            std::vector<bool> undefs;
            LISA* lisa = gda_localjoincount(w, values, undefs, significance_cutoff, nCPUs, permutations,
                                       permutation_method, last_seed_used);
            set_lisa_content(lisa, rst);
            delete lisa;
        }
    }
    return rst;
}

std::vector<int> get_neighbors(const std::string map_uid, const std::string weight_uid, int id)
{
    std::vector<int> nbrs;
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            const  std::vector<long>& _nbrs = w->GetNeighbors(id);
            for (size_t i=0; i<_nbrs.size(); ++i) {
                nbrs.push_back(_nbrs[i]);
            }
        }
    }
    return nbrs;
}

std::vector<int> redcap(const std::string map_uid, const std::string weight_uid,
        int k, std::vector<std::vector<double> > data, std::vector<double> bound_vals, double min_bound,
        std::string method, std::string scale_method, std::string distance_method, int seed, int cpu_threads)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            if (scale_method.empty()) {
                scale_method = "standardize";
            }
            if (distance_method.empty()) {
                distance_method = "euclidean";
            }
            if ( bound_vals.empty()) {
                min_bound = -1;
            }
            std::vector<std::vector<int> > clusters = gda_redcap(k, w, data, scale_method, method, distance_method,
                    bound_vals, min_bound, seed, cpu_threads);
            return GenUtils::flat_2dclusters(w->num_obs, clusters);
        }
    }
    return std::vector<int>();
}

std::vector<double> custom_breaks(const std::string map_uid, int k, std::string col_name,
        std::string break_name)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::vector<bool> undef;
        std::vector<double> data = json_map->GetNumericCol(col_name);
        if (boost::iequals(break_name, "natural_breaks")) {
            return GenUtils::NaturalBreaks(k, data, undef);
        } else if (boost::iequals(break_name, "quantile_breaks")) {
            return GenUtils::QuantileBreaks(k, data, undef);
        } else if (boost::iequals(break_name, "stddev_breaks")) {
            return GenUtils::StddevBreaks(data, undef);
        } else if (boost::iequals(break_name, "hinge15_breaks")) {
            return GenUtils::Hinge15Breaks(data, undef);
        } else if (boost::iequals(break_name, "hinge30_breaks")) {
            return GenUtils::Hinge30Breaks(data, undef);
        }
    }
    return std::vector<double>();
}

std::vector<double> custom_breaks1(const std::string map_uid, int k,
                                  std::string break_name, std::vector<double> data)
{
    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::vector<bool> undef;
        if (boost::iequals(break_name, "natural_breaks")) {
            return GenUtils::NaturalBreaks(k, data, undef);
        } else if (boost::iequals(break_name, "quantile_breaks")) {
            return GenUtils::QuantileBreaks(k, data, undef);
        } else if (boost::iequals(break_name, "stddev_breaks")) {
            return GenUtils::StddevBreaks(data, undef);
        } else if (boost::iequals(break_name, "hinge15_breaks")) {
            return GenUtils::Hinge15Breaks(data, undef);
        } else if (boost::iequals(break_name, "hinge30_breaks")) {
            return GenUtils::Hinge30Breaks(data, undef);
        }
    }
    return std::vector<double>();
}

struct CartogramResult {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> r;

    std::vector<double> get_x() { return x;}
    std::vector<double> get_y() { return y;}
    std::vector<double> get_radius() { return r;}
};

CartogramResult cartogram(const std::string map_uid, std::vector<double> values)
{
    CartogramResult r;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GalWeight* w = (GalWeight*)json_map->CreateRookWeights(1, false, 0);
        int num_obs = w->num_obs;
        const std::vector<gda::PointContents*>& cents = json_map->GetCentroids();
        CartNbrInfo* cart_nbr_info = new CartNbrInfo(w->gal, num_obs);
        std::vector<double> x(num_obs), y(num_obs);
        double orig_data_min=values[0], orig_data_max = values[0];
        for (size_t i=0; i<num_obs; ++i) {
            x[i] = cents[i]->x;
            y[i] = cents[i]->y;
            if (orig_data_min > values[i]) orig_data_min = values[i];
            if (orig_data_max < values[i]) orig_data_max = values[i];
        }
        DorlingCartogram* cartogram = new DorlingCartogram(cart_nbr_info, x, y, values, orig_data_min, orig_data_max);
        cartogram->improve(100);


        for (size_t i=0; i<num_obs; i++) {
            r.x.push_back(cartogram->output_x[i]);
            r.y.push_back(cartogram->output_y[i]);
            r.r.push_back(cartogram->output_radius[i]);
        }
    }
    return r;
}

#ifdef __JSGEODA__
//Using this command to compile
//  emcc --bind -O3 readFile.cpp -s WASM=1 -s TOTAL_MEMORY=268435456 -o api.js --std=c++11
//Note that you need to make sure that there's enough memory available to begin with.
//I got only 16mb without passing the TOTAL_MEMORY setting.
EMSCRIPTEN_BINDINGS(my_module) {

    emscripten::register_vector<std::string>("VectorString");
    emscripten::register_vector<int>("VectorInt");
    emscripten::register_vector<std::vector<int>>("VecVecInt");
    //emscripten::register_vector<std::string>("VectorString");
    emscripten::register_vector<float>("VectorFloat");
    //emscripten::register_vector<float>("vector<float>");
    emscripten::register_vector<double>("VectorDouble");

    emscripten::register_map<std::string, std::vector<float> >("map<string, vector<float>>");

    emscripten::class_<CCentroids>("CCentroids")
        .function("get_x", &CCentroids::get_x)
        .function("get_y", &CCentroids::get_y)
        ;

    emscripten::class_<CartogramResult>("CartogramResult")
        .function("get_x", &CartogramResult::get_x)
        .function("get_y", &CartogramResult::get_y)
        .function("get_radius", &CartogramResult::get_radius)
        ;

    emscripten::class_<LisaResult>("LisaResult")
        .function("significances", &LisaResult::get_sig_local)
        .function("sig_categories", &LisaResult::get_sig_cat)
        .function("clusters", &LisaResult::get_cluster)
        .function("spatial_lags", &LisaResult::get_lag)
        .function("lisa_values", &LisaResult::get_lisa)
        .function("nn", &LisaResult::get_nn)
        .function("labels", &LisaResult::get_labels)
        .function("colors", &LisaResult::get_colors)
        ;

    emscripten::class_<WeightsResult>("WeightsResult")
        .function("get_weight_type", &WeightsResult::get_weight_type)
        .function("get_num_obs", &WeightsResult::get_num_obs)
        .function("get_is_symmetric", &WeightsResult::get_is_symmetric)
        .function("valid", &WeightsResult::valid)
        .function("get_sparsity", &WeightsResult::get_sparsity)
        .function("get_max_nbrs", &WeightsResult::get_max_nbrs)
        .function("get_min_nbrs", &WeightsResult::get_min_nbrs)
        .function("get_median_nbrs", &WeightsResult::get_median_nbrs)
        .function("get_mean_nbrs", &WeightsResult::get_mean_nbrs)
        .function("get_uid", &WeightsResult::get_uid)
        ;

    emscripten::function("new_geojsonmap", &new_geojsonmap);
    emscripten::function("get_num_obs", &get_num_obs);
    emscripten::function("get_map_type", &get_map_type);
    emscripten::function("get_numeric_col", &get_numeric_col);

    emscripten::function("min_distance_threshold", &get_min_dist_threshold);
    emscripten::function("queen_weights", &queen_weights);
    emscripten::function("rook_weights", &rook_weights);
    emscripten::function("knn_weights", &knn_weights);
    emscripten::function("dist_weights", &dist_weights);
    emscripten::function("kernel_weights", &kernel_weights);
    emscripten::function("kernel_bandwidth_weights", &kernel_bandwidth_weights);

    emscripten::function("local_moran", &local_moran);
    emscripten::function("local_moran1", &local_moran1);
    emscripten::function("local_g", &local_g);
    emscripten::function("local_gstar", &local_gstar);
    emscripten::function("local_gstar1", &local_gstar1);
    emscripten::function("local_geary", &local_geary);
    emscripten::function("local_joincount", &local_joincount);

    emscripten::function("redcap", &redcap);

    emscripten::function("custom_breaks", &custom_breaks);
    emscripten::function("custom_breaks1", &custom_breaks1);
    emscripten::function("get_centroids", &get_centroids);
    emscripten::function("cartogram", &cartogram);
    emscripten::function("get_neighbors", &get_neighbors);
}


#endif

int main() {
    std::cout << "print_json" << std::endl;
    std::string file_path = "../data/natregimes.geojson";
    GdaGeojson gda(file_path.c_str());
    geojson_maps["natregimes.geojson"] = &gda;
    WeightsResult r = queen_weights("natregimes.geojson",1,0,0);
    CCentroids c = get_centroids("natregimes.geojson");
    //std::vector<std::string> col_names = {"Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants", "Suicids"};
    std::vector<std::string> col_names = {"hr60", "po60"};
    std::vector<std::vector<int> > clt = redcap("natregimes.geojson", r.uid, 4, col_names, "", -1,"firstorder-singlelinkage");
    return 0;
}

void print_json(char* content) {
	std::cout << "print_json" << std::endl;
	std::cout << content << std::endl;
}
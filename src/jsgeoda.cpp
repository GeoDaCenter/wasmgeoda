#ifdef __JSGEODA__
#include <emscripten/bind.h>
#endif
#include <iostream>
#include <sstream>
#include <vector>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

#include "sa/UniLocalMoran.h"
#include "geojson.h"
#include "gda_weights.h"
#include "gda_sa.h"

extern "C" {
	void print_json(char* content);
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
void new_geojsonmap(std::string file_name, const int & addr, const size_t & len){
	//We get out pointer as a plain int from javascript
	//We use a reinterpret_cast to turn our plain int into a uint8_t pointer. After
	//which we can play with the data just like we would normally.
	char* data = reinterpret_cast<char*>(addr);

	std::cout << "new_geojsonmap():" << file_name << std::endl;
	GdaGeojson *json_map = new GdaGeojson(file_name.c_str(), data);

	// store globally, has to be release by calling free_geojsonmap()
	geojson_maps[file_name] = json_map;
}


extern "C" 
int get_num_obs(char* _map_uid) {
	std::string map_uid(_map_uid);
	std::cout << "get_num_obs()" << map_uid << std::endl;
	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		return json_map->GetNumObs();
	}
	return 0;
}


extern "C"
int get_map_type(char* _map_uid) {
	std::string map_uid(_map_uid);
	std::cout << "get_map_type()" << map_uid << std::endl;
	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		return json_map->GetMapType();
	}
	return 0;
}

extern "C" 
const char* queen_weights(char* _map_uid,
	int order,
	int include_lower_order,
 	float precision_threshold)
{
	std::string map_uid(_map_uid);
	std::cout << "queen_weights()" << map_uid << std::endl;

	GdaGeojson *json_map = geojson_maps[map_uid];
	if (json_map) {
		GeoDaWeight *w = json_map->CreateQueenWeights(order, include_lower_order, precision_threshold);
		if (w) return w->uid.c_str();
	}
	return "";
}

extern "C"
const char* rook_weights(char* _map_uid,
                          int order,
                          int include_lower_order,
                          float precision_threshold)
{
    std::string map_uid(_map_uid);
    std::cout << "rook_weights()" << map_uid << std::endl;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateRookWeights(order, include_lower_order, precision_threshold);
        double sparsity =w->GetSparsity();
        std::cout << "sparsity:" << sparsity << std::endl;
        std::stringstream rst;
        rst << sparsity;
        return rst.str().c_str();
    }
    return "";
}

extern "C"
const char* knn_weights(char* _map_uid,
                        unsigned int k,
                        float power,
                        bool is_inverse,
                        bool is_arc,
                        bool is_mile)
{
    std::string map_uid(_map_uid);
    std::cout << "knn_weights()" << map_uid << std::endl;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateKnnWeights(k, power, is_inverse, is_arc, is_mile);
        double sparsity =w->GetSparsity();
        std::cout << "sparsity:" << sparsity << std::endl;
        std::stringstream rst;
        rst << sparsity;
        return rst.str().c_str();
    }
    return "";
}

extern "C"
const char* dist_weights(char* _map_uid,
                         float dist_thres,
                        float power,
                        bool is_inverse,
                        bool is_arc,
                        bool is_mile)
{
    std::string map_uid(_map_uid);
    std::cout << "distance_weights()" << map_uid << std::endl;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->CreateDistanceWeights(dist_thres, power, is_inverse, is_arc, is_mile);
        double sparsity =w->GetSparsity();
        std::cout << "sparsity:" << sparsity << std::endl;
        std::stringstream rst;
        rst << sparsity;
        return rst.str().c_str();
    }
    return "";
}

extern "C"
const char* kernel_weights(char* _map_uid,
        unsigned int k,
        char* kernel,
        bool adaptive_bandwidth,
        bool use_kernel_diagonals,
        bool is_arc,
        bool is_mile)
{
    std::string map_uid(_map_uid);
    std::cout << "kernel_weights()" << map_uid << std::endl;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelWeights(k, kernel_str, adaptive_bandwidth, use_kernel_diagonals, is_arc, is_mile);
        double sparsity =w->GetSparsity();
        std::cout << "sparsity:" << sparsity << std::endl;
        std::stringstream rst;
        rst << sparsity;
        return rst.str().c_str();
    }
    return "";
}

extern "C"
const char* kernel_bandwidth_weights(char* _map_uid,
                                     float dist_thres,
                           char* kernel,
                           bool use_kernel_diagonals,
                           bool is_arc,
                           bool is_mile)
{
    std::string map_uid(_map_uid);
    std::cout << "kernel_bandwidth_weights()" << map_uid << std::endl;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        std::string kernel_str(kernel);
        GeoDaWeight *w = json_map->CreateKernelWeights(dist_thres, kernel_str, use_kernel_diagonals, is_arc, is_mile);
        double sparsity =w->GetSparsity();
        std::cout << "sparsity:" << sparsity << std::endl;
        std::stringstream rst;
        rst << sparsity;
        return rst.str().c_str();
    }
    return "";
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
    std::vector<double> sig_local_vec;
    std::vector<int> sig_cat_vec;
    std::vector<int> cluster_vec;
    std::vector<double> lag_vec;
    std::vector<double> lisa_vec;
    std::vector<int> nn_vec;
    std::vector<std::string> labels;
    std::vector<std::string> colors;

    std::vector<double> get_sig_local() { return sig_local_vec;}
    std::vector<int> get_sig_cat() { return sig_cat_vec;}
    std::vector<int>  get_cluster() { return cluster_vec;}
    std::vector<double>  get_lag() { return lag_vec;}
    std::vector<double>  get_lisa() { return lisa_vec;}
    std::vector<int>  get_nn() { return nn_vec;}
    std::vector<std::string>  get_labels() { return labels;}
    std::vector<std::string>  get_colors() { return colors;}
};

LisaResult local_moran(const std::string map_uid,
        const std::string weight_uid, std::vector<double> values)
{
    LisaResult rst;

    GdaGeojson *json_map = geojson_maps[map_uid];
    if (json_map) {
        GeoDaWeight *w = json_map->GetWeights(weight_uid);
        if (w) {
            UniLocalMoran* moran = gda_localmoran(w, values);
            rst.sig_local_vec = moran->GetLocalSignificanceValues();
            rst.sig_cat_vec = moran->GetSigCatIndicators();
            rst.cluster_vec = moran->GetClusterIndicators();
            delete moran;
        }
    }

    return rst;
}


#ifdef __JSGEODA__
//Using this command to compile
//  emcc --bind -O3 readFile.cpp -s WASM=1 -s TOTAL_MEMORY=268435456 -o api.js --std=c++11
//Note that you need to make sure that there's enough memory available to begin with.
//I got only 16mb without passing the TOTAL_MEMORY setting.
EMSCRIPTEN_BINDINGS(my_module) {

    emscripten::register_vector<std::string>("VectorString");
    emscripten::register_vector<int>("VectorInt");
    //emscripten::register_vector<std::string>("VectorString");
    emscripten::register_vector<float>("VectorFloat");
    //emscripten::register_vector<float>("vector<float>");
    emscripten::register_vector<double>("VectorDouble");

    emscripten::register_map<std::string, std::vector<float> >("map<string, vector<float>>");

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

    emscripten::function("new_geojsonmap", &new_geojsonmap);
    emscripten::function("local_moran", &local_moran);
    emscripten::function("get_numeric_col", &get_numeric_col);
}

int main() {
    std::cout << "print_json" << std::endl;
	return 0;
}
#endif

void print_json(char* content) {
	std::cout << "print_json" << std::endl;
	std::cout << content << std::endl;
}
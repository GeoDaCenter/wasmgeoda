#ifdef __JSGEODA__
#include <emscripten/bind.h>
#endif
#include <iostream>
#include <sstream>
#include <vector>

#include <rapidjson/document.h>
#include <rapidjson/writer.h>
#include <rapidjson/stringbuffer.h>

#include "geojson.h"
#include "gda_weights.h"

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

#ifdef __JSGEODA__
//Using this command to compile
//  emcc --bind -O3 readFile.cpp -s WASM=1 -s TOTAL_MEMORY=268435456 -o api.js --std=c++11
//Note that you need to make sure that there's enough memory available to begin with.
//I got only 16mb without passing the TOTAL_MEMORY setting.
EMSCRIPTEN_BINDINGS(my_module) {
  emscripten::function("new_geojsonmap", &new_geojsonmap);
}
#endif

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
		double sparsity =w->GetSparsity();
		std::cout << "sparsity:" << sparsity << std::endl;
		std::stringstream rst; 
		rst << sparsity;
		return rst.str().c_str();
	}
	return "";
}

void local_moran(const std::string& map_uid, 
	const std::string& weight_uid,  
	double* values)
{

}

#ifdef __JSGEODA__
int main() {
    std::cout << "print_json" << std::endl;
	return 0;
}
#endif

void print_json(char* content) {
	std::cout << "print_json" << std::endl;
	std::cout << content << std::endl;
}
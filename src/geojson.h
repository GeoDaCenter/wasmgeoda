#ifndef JSGEODA_GEOJSON
#define JSGEODA_GEOJSON

#include <vector>
#include <map>
#include <string>
#include <rapidjson/document.h>

#include "weights/GeodaWeight.h"
#include "geofeature.h"

class GdaGeojson
{
public:
    // default constructor for std::vector and std::map
    GdaGeojson();

    GdaGeojson(const std::string& file_path);
    
    GdaGeojson(const char* file_name, const char* in_content);

    virtual ~GdaGeojson();

    void Read(const char* file_name, const char* in_content);

    int GetNumObs();

    gda::ShapeType GetMapType();

    const std::vector<gda::Point>& GetCentroids();

    gda::MainMap& GetMainMap();

    // weights related functions:
    GeoDaWeight* CreateQueenWeights(unsigned int order=1, 
        bool include_lower_order = false,
 	    double precision_threshold = 0);

    GeoDaWeight* CreateRookWeights(unsigned int order=1, 
        bool include_lower_order = false,
 	    double precision_threshold = 0);

protected:
    gda::MainMap main_map;

    std::map<std::string, GeoDaWeight*> weights_dict;

     std::vector<gda::Point> centroids;

    // read geojson related functions:
    void init();

    void checkType(const rapidjson::Document& json);

    void readFeatureCollection(const rapidjson::Value& features);

    void createGeometryFeature(size_t i, const rapidjson::Value& geom);

    void addPoint(const rapidjson::Value &coord);

    void addMultiPoints(const rapidjson::Value &coords);

    void addPolygon(const rapidjson::Value &coords);

    void addMultiPolygons(const rapidjson::Value &coords);
};

#endif
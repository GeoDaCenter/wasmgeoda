#ifndef JSGEODA_GEOJSON
#define JSGEODA_GEOJSON

#include <vector>
#include <map>
#include <string>
#include <rapidjson/document.h>

#include "../libgeoda_src/weights/GeodaWeight.h"
#include "../libgeoda_src/geofeature.h"
#include "../libgeoda_src/gda_interface.h"

class GdaGeojson : public AbstractGeoDa
{
public:
    // default constructor for std::vector and std::map
    GdaGeojson();

    GdaGeojson(const std::string& file_path);
    
    GdaGeojson(const char* file_name, const char* in_content);

    virtual ~GdaGeojson();

    void Read(const char* file_name, const char* in_content);

    virtual int GetNumObs() const;

    virtual const std::vector<gda::PointContents*>& GetCentroids();

    virtual int GetMapType();

    virtual std::string GetMapTypeName();

    virtual gda::MainMap& GetMainMap();

    std::vector<double> GetNumericCol(std::string col_name);

    std::vector<std::string> GetStringCol(std::string col_name);

    bool IsNumericCol(std::string col_name);

    // weights related functions:
    GeoDaWeight* CreateQueenWeights(unsigned int order=1, 
        bool include_lower_order = false,
 	    double precision_threshold = 0);

    GeoDaWeight* CreateRookWeights(unsigned int order=1, 
        bool include_lower_order = false,
 	    double precision_threshold = 0);

    GeoDaWeight* CreateKnnWeights(unsigned int k,
        double power = 1.0,
        bool is_inverse = false,
        bool is_arc = false,
        bool is_mile = true);

    GeoDaWeight* CreateDistanceWeights(double dist_thres,
        double power = 1.0,
        bool is_inverse = false,
        bool is_arc = false,
        bool is_mile = true);

    GeoDaWeight* CreateKernelKnnWeights(unsigned int k,
        const std::string& kernel,
        bool adaptive_bandwidth = false,
        bool use_kernel_diagonals = false,
        bool is_arc = false,
        bool is_mile = true);

    GeoDaWeight* CreateKernelWeights(double dist_thres,
        const std::string& kernel,
        bool use_kernel_diagonals = false,
        bool is_arc = false,
        bool is_mile = true);

    GeoDaWeight* GetWeights(const std::string& w_uid) {
        return weights_dict[w_uid];
    }

    double GetMinDistanceThreshold(bool is_arc, bool is_mile);

    std::string GetFilePath() const { return file_path; }

protected:
    std::string file_path;

    gda::MainMap main_map;

    std::map<std::string, std::vector<double> > data_numeric;

    std::map<std::string, std::vector<std::string> > data_string;

    std::map<std::string, GeoDaWeight*> weights_dict;

    std::vector<gda::PointContents*> centroids;

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
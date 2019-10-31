#include <iostream>
#include <sstream> 
#include <algorithm>
#include <limits>
#include <boost/algorithm/string.hpp>

#include "shape/centroid.h"
#include "gda_weights.h"
#include "geojson.h"

using error = std::runtime_error;

GdaGeojson::GdaGeojson()
{

}

GdaGeojson::GdaGeojson(const std::string& file_path)
        : GdaGeojson()
{
#ifndef __JSGEODA__
    std::string filename = file_path.substr(file_path.find_last_of("/") + 1);

    FILE *fp;
    long lSize;
    char *buffer;

    fp = fopen (file_path.c_str(), "rb" );
    if( !fp ) perror("blah.txt"),exit(1);

    fseek( fp , 0L , SEEK_END);
    lSize = ftell( fp );
    rewind( fp );

    /* allocate memory for entire content */
    buffer = (char*)calloc( 1, lSize+1 );
    if( !buffer ) fclose(fp),fputs("memory alloc fails",stderr),exit(1);

    /* copy the file into the buffer */
    if( 1!=fread( buffer , lSize, 1 , fp) )
        fclose(fp),free(buffer),fputs("entire read fails",stderr),exit(1);

    /* do your work here, buffer is a string contains the whole text */
    this->Read(filename.c_str(), buffer);

    fclose(fp);
    free(buffer);
#endif
}

GdaGeojson::GdaGeojson(const char* file_name, const char* in_content)
: GdaGeojson()
{
    this->Read(file_name, in_content);
}

GdaGeojson::~GdaGeojson()
{
    // free memory
    std::map<std::string, GeoDaWeight*>::iterator it;
    for (it = weights_dict.begin(); it != weights_dict.end(); ++it) {
        delete it->second;
    }
    weights_dict.clear();
}

gda::MainMap& GdaGeojson::GetMainMap()
{
    return this->main_map;
}

int GdaGeojson::GetNumObs()
{
    return this->main_map.num_obs;
}

gda::ShapeType GdaGeojson::GetMapType()
{
    return this->main_map.shape_type;
}

const std::vector<gda::Point>& GdaGeojson::GetCentroids()
{
    if (this->centroids.empty()) {
        if (this->main_map.shape_type == gda::POINT_TYP) {
            this->centroids.resize(this->main_map.num_obs);
            for (size_t i=0; i<this->centroids.size(); ++i) {
                gda::PointContents* pt = (gda::PointContents*)this->main_map.records[i];
                this->centroids[i].x = pt->x;
                this->centroids[i].y = pt->y;
            }
        } else if (this->main_map.shape_type == gda::POLYGON) {
            this->centroids.resize(this->main_map.num_obs);
            for (size_t i=0; i<this->centroids.size(); ++i) {
                gda::PolygonContents* poly = (gda::PolygonContents*)this->main_map.records[i];
                Centroid cent(poly);
                cent.getCentroid(this->centroids[i]);
            }
        }
    }
    return this->centroids;
}


double GdaGeojson::GetMinDistanceThreshold(bool is_arc, bool is_mile)
{
    return gda_min_distthreshold(this, is_arc, is_mile);
}

GeoDaWeight* GdaGeojson::CreateQueenWeights(unsigned int order, 
        bool include_lower_order,
 	    double precision_threshold)
{
    std::stringstream w_uid;
    w_uid << "w_queen";
    w_uid << order;
    w_uid << include_lower_order;
    w_uid << precision_threshold;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        w = gda_queen_weights(this, order, include_lower_order, precision_threshold);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

GeoDaWeight* GdaGeojson::CreateRookWeights(unsigned int order, 
        bool include_lower_order,
 	    double precision_threshold)
{
    std::stringstream w_uid;
    w_uid << "w_rook";
    w_uid << order;
    w_uid << include_lower_order;
    w_uid << precision_threshold;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        w = gda_rook_weights(this, order, include_lower_order, precision_threshold);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

GeoDaWeight* GdaGeojson::CreateKnnWeights(unsigned int k,
                              double power,
                              bool is_inverse,
                              bool is_arc,
                              bool is_mile)
{
    std::stringstream w_uid;
    w_uid << "w_knn";
    w_uid << k;
    w_uid << power;
    w_uid << is_inverse;
    w_uid << is_arc;
    w_uid << is_mile;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        gda_knn_weights(this, k, power, is_inverse, is_arc, is_mile);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

GeoDaWeight* GdaGeojson::CreateDistanceWeights(double dist_thres,
                                          double power,
                                          bool is_inverse,
                                          bool is_arc,
                                          bool is_mile)
{
    std::stringstream w_uid;
    w_uid << "w_dist";
    w_uid << dist_thres;
    w_uid << power;
    w_uid << is_inverse;
    w_uid << is_arc;
    w_uid << is_mile;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        gda_distance_weights(this, dist_thres, "", power, is_inverse, is_arc, is_mile);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

GeoDaWeight* GdaGeojson::CreateKernelWeights(double dist_thres,
                                             const std::string& kernel,
                                             bool use_kernel_diagonals,
                                             bool is_arc, bool is_mile)
{
    std::stringstream w_uid;
    w_uid << "w_kernel";
    w_uid << dist_thres;
    w_uid << kernel;
    w_uid << use_kernel_diagonals;
    w_uid << is_arc;
    w_uid << is_mile;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        gda_distance_weights(this, dist_thres, "", 1.0, false, is_arc, is_mile, kernel, use_kernel_diagonals);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

GeoDaWeight* GdaGeojson::CreateKernelWeights(unsigned int k,
                                             const std::string& kernel,
                                             bool adaptive_bandwidth,
                                             bool use_kernel_diagonals,
                                             bool is_arc, bool is_mile)
{
    std::stringstream w_uid;
    w_uid << "w_kernel";
    w_uid << k;
    w_uid << kernel;
    w_uid << adaptive_bandwidth;
    w_uid << use_kernel_diagonals;
    w_uid << is_arc;
    w_uid << is_mile;
    std::string w_uid_str = w_uid.str();
    GeoDaWeight* w = 0;
    if (this->weights_dict.find(w_uid_str) != this->weights_dict.end()) {
        w = this->weights_dict[w_uid.str()];
    } else {
        gda_knn_weights(this, k, 1.0, false, is_arc, is_mile, kernel, 0.0, adaptive_bandwidth, use_kernel_diagonals);
        w->uid = w_uid_str;
        this->weights_dict[w_uid.str()] = w;
    }
    return w;
}

std::vector<double> GdaGeojson::GetNumericCol(std::string col_name)
{
    if (data_numeric.find(col_name) == data_numeric.end()) {
        std::cout << "not found col" <<std::endl;
        return std::vector<double>();
    }
    return data_numeric[col_name];
}


void GdaGeojson::Read(const char* file_name, const char* in_content)
{
    rapidjson::Document json;
    json.Parse(in_content);

    if (!json.IsObject())
        throw error("Geometry must be an object");

    const auto &json_end = json.MemberEnd();

    const auto &features_itr = json.FindMember("features");
    if (features_itr == json_end)
        throw error("Content of features not found");

    const auto &features = features_itr->value;

    this->readFeatureCollection(features);
}

void GdaGeojson::readFeatureCollection(const rapidjson::Value& features)
{
    /*
    {
        "type": "Feature",
        "geometry": {
            "type": "Point",
            "coordinates": [125.6, 10.1]
        },
        "properties": {
            "name": "Dinagat Islands"
        }
    }
    */
    this->main_map.bbox_x_min = std::numeric_limits<double>::max();
    this->main_map.bbox_y_min = std::numeric_limits<double>::max();
    this->main_map.bbox_x_max = std::numeric_limits<double>::lowest();
    this->main_map.bbox_y_max = std::numeric_limits<double>::lowest();
    if (!features.IsArray()) 
        throw error("FeatureCollection is empty");

    const auto &fa = features.GetArray();
    this->main_map.num_obs = fa.Size();
    for (size_t i=0; i<fa.Size(); ++i) {
        const rapidjson::Value &geom = fa[i]["geometry"];
        //std::cout << fa[i]["properties"]["POLY_ID"].GetInt() <<std::endl;
        // get data
        const rapidjson::Value &var_names = fa[i]["properties"];
        for (rapidjson::Value::ConstMemberIterator iter = var_names.MemberBegin();
            iter != var_names.MemberEnd(); ++iter)
        {
            std::string var_name = iter->name.GetString();
            if (iter->value.IsNumber()) {
                data_numeric[var_name].push_back(iter->value.GetDouble());
            } else {
                data_string[var_name].push_back(iter->value.GetString());
            }
        }

        // get geometry
        if (geom.IsNull()) {
            // null geometry
            this->main_map.records.push_back(new gda::NullShapeContents());
        } else {
            this->createGeometryFeature(i, geom);
        }
    }
}

void GdaGeojson::createGeometryFeature(size_t i, const rapidjson::Value& geom)
{
    /*
    "geometry": {
        "type": "Point",
        "coordinates": [125.6, 10.1]
    },
    */
    if (geom["type"].IsNull()) 
        throw error("geometry::type is NULL");

    if (geom["coordinates"].IsNull()) {
        // empty geometry
        this->main_map.records.push_back(new gda::NullShapeContents());
        return;
    }

    const rapidjson::Value &coords = geom["coordinates"];

    std::string geom_type = geom["type"].GetString();
    if (boost::iequals(geom_type, "Point")) {
        this->addPoint(coords);
        this->main_map.shape_type = gda::POINT_TYP;

    } else if (boost::iequals(geom_type, "MultiPoint")) {
        this->addMultiPoints(coords);
        this->main_map.shape_type = gda::POINT_TYP;

    } else if (boost::iequals(geom_type, "Polygon")) {
        this->addPolygon(coords);
        this->main_map.shape_type = gda::POLYGON;

    } else if (boost::iequals(geom_type, "MultiPolygon")) {
        this->addMultiPolygons(coords);
        this->main_map.shape_type = gda::POLYGON;
        
    } else {
        throw error("Geometry::type (Line) is not supported");
    }
}

void GdaGeojson::addPoint(const rapidjson::Value &coord) 
{
    if (coord.Size() < 2) {
        this->main_map.records.push_back(new gda::NullShapeContents());
    } else {
        gda::PointContents* pt = new gda::PointContents();
        pt->x = coord[0].GetDouble();
        pt->y = coord[1].GetDouble();
        this->main_map.set_bbox(pt->x,  pt->y);
        this->main_map.records.push_back(pt);
    }
}

void GdaGeojson::addMultiPoints(const rapidjson::Value &coords) 
{
    // geoda doesn't support multi-points feature, and it is treated by using
    // the first point, and an warning will be raised
    if (coords.IsArray()) {
        if (coords.Size() == 0) {
            this->main_map.records.push_back(new gda::NullShapeContents());
        } else if (coords[0].IsArray()) {
            this->addPoint(coords[0]);
        } else if (coords[0].IsNumber()) {
            this->addPoint(coords);
        }
    }  else {
        this->main_map.records.push_back(new gda::NullShapeContents());
    }
}

void GdaGeojson::addPolygon(const rapidjson::Value &coords) 
{
    // in some cases,
    // [
    //	[
    //		[-80.874755859375,25.8012447357178],[-80.8742065429688,25.9828872680664],[-80.6838607788086,25.9843692779541],[-80.6811676025391,25.960147857666],
    //		[-80.2978363037109,25.9572772979736],[-80.2978515625,25.9739627838135],[-80.1280136108398,25.9771671295166],[-80.1933288574219,25.7596549987793],
    //	]
    //]
    // it could be polygon with holes, e.g. NAT "POLY_ID": 1137, "NAME": "Frederick", "STATE_NAME": "Virginia",
    // [
    //	[
    //		[-78.1545867919922,39.0405921936035],[-78.1680374145508,39.0219383239746],[-78.3112487792969,39.0107421875],[-78.3190536499023,39.021484375],[-78.3138656616211,39.0337600708008],
    //		[-78.3382339477539,39.0413055419922],[-78.3365478515625,39.0490188598633],[-78.3490600585938,39.0552787780762],[-78.3275451660156,39.0907821655273],[-78.3442840576172,39.1029052734375],
    //		[-78.354248046875,39.093318939209],[-78.3922805786133,39.1003036499023],[-78.3968505859375,39.0866737365723],[-78.4345550537109,39.0682640075684],[-78.4559478759766,39.027759552002],
    //		[-78.5369186401367,39.0570259094238],[-78.5018692016602,39.093578338623],[-78.4855194091797,39.1118392944336],[-78.4482498168945,39.1189308166504],[-78.4308395385742,39.1485214233398],
    //		[-78.4026336669922,39.1704902648926],[-78.4243392944336,39.1975250244141],[-78.42333984375,39.2120399475098],[-78.3993988037109,39.2448501586914],[-78.413818359375,39.257438659668],
    //		[-78.3411178588867,39.3413581848145],[-78.3442001342773,39.3508567810059],[-78.3657455444336,39.3615875244141],[-78.3505020141602,39.380729675293],[-78.3478164672852,39.456901550293],
    //		[-78.2771530151367,39.4233665466309],[-78.2297821044922,39.3910140991211],[-78.0336074829102,39.2655372619629],[-78.0997543334961,39.1446685791016],[-78.1089096069336,39.1042823791504],
    //		[-78.1511688232422,39.056022644043],[-78.1545867919922,39.0405921936035]
    //	],
    //	[
    //		[-78.1327972412109,39.1916427612305],[-78.1827621459961,39.2027130126953],[-78.2050399780273,39.1731262207031],
    //		[-78.2054824829102,39.1577110290527],[-78.1625823974609,39.1384582519531],[-78.1396865844727,39.164867401123],[-78.1327972412109,39.1916427612305]
    //	]
    //]
    // in other cases: geom.coordinates is Array[24]: [[1,2],[3,4]]
    // [
    //	 [-80.874755859375,25.8012447357178],[-80.8742065429688,25.9828872680664],[-80.6838607788086,25.9843692779541],[-80.6811676025391,25.960147857666],
    //	 [-80.2978363037109,25.9572772979736],[-80.2978515625,25.9739627838135],[-80.1280136108398,25.9771671295166],[-80.1933288574219,25.7596549987793],
    // ]
    if (coords.IsArray() && coords.Size() > 0) {
        gda::PolygonContents *poly = new gda::PolygonContents();
        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::lowest();
        double maxy = std::numeric_limits<double>::lowest();
        double x, y;

        if (coords[0][0].IsNumber()) {
            // second case
            const rapidjson::Value &xys = coords;
            poly->num_parts = 1;
            poly->num_points = xys.Size();
            poly->parts.push_back(0);
            poly->points.resize(xys.Size());
            poly->holes.push_back(false);

            for (size_t i=0; i<xys.Size(); ++i) {
                const rapidjson::Value &xy = xys[i];
                x = xy[0].GetDouble();
                y = xy[1].GetDouble();
                poly->points[i].x = x;
                poly->points[i].y = y;
                if ( x < minx ) minx = x;
                if ( x >= maxx ) maxx = x;
                if ( y < miny ) miny = y;
                if ( y >= maxy ) maxy = y;
            }


        } else {
            // first case
            const rapidjson::Value &parts = coords;
            size_t n_parts = parts.Size();
            for (size_t i=0; i<n_parts; ++i) {
                const rapidjson::Value &part = parts[i];
                poly->num_parts += 1;
                poly->parts.push_back(poly->num_points);
                bool is_hole = i > 0 ? true : false;
                poly->holes.push_back(is_hole);

                for (size_t j=0; j< part.Size(); ++j) {
                    const rapidjson::Value& xy = part[j];
                    x = xy[0].GetDouble();
                    y = xy[1].GetDouble();

                    poly->points.push_back(gda::Point(x,y));
                    poly->num_points += 1;

                    if ( x < minx ) minx = x;
                    if ( x >= maxx ) maxx = x;
                    if ( y < miny ) miny = y;
                    if ( y >= maxy ) maxy = y;
                }
            }
        }

        poly->box.resize(4);
        poly->box[0] = minx;
        poly->box[1] = miny;
        poly->box[2] = maxx;
        poly->box[3] = maxy;

        this->main_map.set_bbox(minx, miny);
        this->main_map.set_bbox(maxx, maxy);
        this->main_map.records.push_back(poly);
    } else {
        this->main_map.records.push_back(new gda::NullShapeContents());
    }
}

void GdaGeojson::addMultiPolygons(const rapidjson::Value &coords) 
{
    // [
    //	[
    //		[
    //			[ -78.839393615722699, 38.042407989502003 ],
    //			[ -78.883575439453097, 38.0301322937012 ], [ -78.898597717285199, 37.990524291992202 ],
    //		]
    //	],
    //	[
    //		[
    //			[ -78.857666015625, 38.083003997802699 ], [ -78.878311157226605, 38.092292785644503 ],
    //		]
    //	],
    //	[
    //		[
    //			[ -79.040824890136705, 38.145206451416001 ], [ -79.046768188476605, 38.149665832519503 ],
    //		]
    //	]
    //]
    if (coords.IsArray() && coords.Size() > 0) {
        const rapidjson::Value &parts = coords;

        gda::PolygonContents *poly = new gda::PolygonContents();
        poly->num_parts = 0;
        poly->num_points = 0;

        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::lowest();
        double maxy = std::numeric_limits<double>::lowest();
        double x, y;

        for (size_t p=0; p< parts.Size(); ++p) {
            const rapidjson::Value &part = parts[p];
            size_t sub_parts = part.Size();


            for (size_t sp = 0; sp < sub_parts; ++sp) {
                const rapidjson::Value &xys = part[sp];
                size_t n_xy = xys.Size();

                poly->parts.push_back(poly->num_points);
                poly->num_parts += 1;
                bool is_hole = sp > 0 ? true : false;
                poly->holes.push_back(is_hole);

                for (size_t i = 0; i < n_xy; ++i) {
                    const rapidjson::Value &xy = xys[i];
                    x = xy[0].GetDouble();
                    y = xy[1].GetDouble();

                    poly->points.push_back(gda::Point(x, y));
                    poly->num_points += 1;

                    if (x < minx) minx = x;
                    if (x >= maxx) maxx = x;
                    if (y < miny) miny = y;
                    if (y >= maxy) maxy = y;
                }
            }
        }
        poly->box.resize(4);
        poly->box[0] = minx;
        poly->box[1] = miny;
        poly->box[2] = maxx;
        poly->box[3] = maxy;

        this->main_map.set_bbox(minx, miny);
        this->main_map.set_bbox(maxx, maxy);

        this->main_map.records.push_back(poly);
    } else {
        this->main_map.records.push_back(new gda::NullShapeContents());
    }
}

void GdaGeojson::checkType(const rapidjson::Document& json)
{
    const auto &type_itr = json.FindMember("type");
    const auto &json_end = json.MemberEnd();

    if (type_itr == json_end)
        throw error("Geojson must have a type property");

    const auto &type = type_itr->value;
    if (type != "FeatureCollection") {
        throw error("Only FeatureCollection type is supported");
    }
}
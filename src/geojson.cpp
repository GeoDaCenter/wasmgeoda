#include <iostream>
#include <sstream> 
#include <algorithm>
#include <limits>

#include "utils.h"
#include "geojson.h"
#include "gda_weights.h"

using error = std::runtime_error;

GdaGeojson::GdaGeojson()
{

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
    return std::vector<gda::Point>();
}

GeoDaWeight* GdaGeojson::CreateQueenWeights(unsigned int order, 
        bool include_lower_order,
 	    double precision_threshold)
{
    GeoDaWeight* w = gda_queen_weights(this, order, include_lower_order, precision_threshold);
    std::stringstream w_uid;
    w_uid << "w";
    w_uid << order;
    w_uid << include_lower_order;
    w_uid << precision_threshold;
    this->weights_dict[w_uid.str()] = w;
    return w;
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
    this->main_map.bbox_x_max = std::numeric_limits<double>::min();
    this->main_map.bbox_y_max = std::numeric_limits<double>::min();
    if (!features.IsArray()) 
        throw error("FeatureCollection is empty");

    const auto &fa = features.GetArray();
    this->main_map.num_obs = fa.Size();
    for (size_t i=0; i<fa.Size(); ++i) {
        const rapidjson::Value &geom = fa[i]["geometry"];
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
    if (iequals(geom_type, "Point")) {
        this->addPoint(coords);
        this->main_map.shape_type = gda::POINT_TYP;

    } else if (iequals(geom_type, "MultiPoint")) {
        this->addMultiPoints(coords);
        this->main_map.shape_type = gda::POINT_TYP;

    } else if (iequals(geom_type, "Polygon")) {
        this->addPolygon(coords);
        this->main_map.shape_type = gda::POLYGON;

    } else if (iequals(geom_type, "MultiPolygon")) {
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
    // in some cases: geom.coordinates is [Array[24]] [ [[1,2],[3,4]] ]
    // in other cases: geom.coordinates is Array[24]: [[1,2],[3,4]]
    if (coords.IsArray() && coords.Size() > 0) {
        const rapidjson::Value &xys = coords[0][0].IsArray() ? coords[0] : coords;
        gda::PolygonContents *poly = new gda::PolygonContents();
        poly->num_parts = 1;
        poly->num_points = xys.Size();
        poly->parts.push_back(poly->num_points);
        poly->points.resize(xys.Size());

        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::min();
        double maxy = std::numeric_limits<double>::min();
        double x, y;
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
    // [Array[24]] [ [ [[1,2],[3,4]], [[5,6],[7,8]] ] ]
    if (coords.IsArray() && coords.Size() > 0) {
        const rapidjson::Value &parts = coords[0];

        gda::PolygonContents *poly = new gda::PolygonContents();
        poly->num_parts = parts.Size();
        poly->num_points = 0;
        poly->parts.resize(poly->num_parts);

        double minx = std::numeric_limits<double>::max();
        double miny = std::numeric_limits<double>::min();
        double maxx = std::numeric_limits<double>::max();
        double maxy = std::numeric_limits<double>::min();
        double x, y;
        for (size_t part=0; part < poly->num_parts; ++part) {
            const rapidjson::Value &xys = parts[part];
            poly->parts[part] = xys.Size();

            for (size_t i=0; i<xys.Size(); ++i) {
                const rapidjson::Value &xy = xys[i];
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
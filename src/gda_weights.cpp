#include <vector>
#include <set>

#include "weights/PointsToContigWeights.h"
#include "weights/PolysToContigWeights.h"
#include "weights/GalWeight.h"
#include "weights/GeodaWeight.h"
#include "geojson.h"
#include "gda_weights.h"

GeoDaWeight* contiguity_weights(bool is_queen,
                                GdaGeojson* geoda,
                                unsigned int order,
                                bool include_lower_order,
                                double precision_threshold)
{
    if (geoda == 0) return 0;

    int num_obs = geoda->GetNumObs();
    std::cout << "contiguity_weights()" << num_obs << std::endl;
    GalWeight* poW = new GalWeight;
    poW->num_obs = num_obs;
    poW->is_symmetric = true;
    poW->symmetry_checked = true;

    std::cout << "contiguity_weights()maptype" << geoda->GetMapType() << std::endl;
    if (geoda->GetMapType() == gda::POINT_TYP) {
        std::vector<std::set<int> > nbr_map;
        const std::vector<gda::Point>& centroids = geoda->GetCentroids();
        std::vector<double> x(num_obs), y(num_obs);
        for (size_t i=0; i<num_obs; ++i) {
            x[i] = centroids[i].x;
            y[i] = centroids[i].y;
        }
        gda::PointsToContiguity(x, y, is_queen, nbr_map);
        poW->gal = Gda::NeighborMapToGal(nbr_map);

    } else if (geoda->GetMapType() == gda::POLYGON) {
        std::cout << "contiguity_weights()POLYGON" << precision_threshold << std::endl;
        poW->gal = PolysToContigWeights(geoda->GetMainMap(), is_queen, precision_threshold);
        if (order > 1) {
            Gda::MakeHigherOrdContiguity(order, num_obs, poW->gal, include_lower_order);
        }

    } else {
        // line_type not supported yet, should be detected at script side
        delete poW;
        return 0;
    }

    poW->GetNbrStats();
    return (GeoDaWeight*)poW;
}

GeoDaWeight* gda_queen_weights(GdaGeojson* json_map,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = true;
    return contiguity_weights(is_queen, json_map, order, include_lower_order, precision_threshold);
}

GeoDaWeight* gda_rook_weights(GdaGeojson* json_map,
                               unsigned int order,
                               bool include_lower_order,
                               double precision_threshold)
{
    bool is_queen = false;
    return contiguity_weights(is_queen, json_map, order, include_lower_order, precision_threshold);
}
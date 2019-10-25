#ifndef __JSGEODSA_GDA_WEIGHTS__
#define __JSGEODSA_GDA_WEIGHTS__

#include <string>

class GdaGeojson;
class GeoDaWeight;

// APIs of weights creation
/**
 *
 * @param geoda
 * @param polyid
 * @param order
 * @param include_lower_order
 * @param precision_threshold
 * @return
 */
GeoDaWeight* gda_queen_weights(GdaGeojson* json_map, 
                               unsigned int order=1,
                               bool include_lower_order = false,
                               double precision_threshold = 0);

GeoDaWeight* gda_rook_weights(GdaGeojson* json_map, 
                               unsigned int order=1,
                               bool include_lower_order = false,
                               double precision_threshold = 0);
#endif

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

GeoDaWeight* gda_knn_weights(GdaGeojson* geoda, unsigned int k,
                             double power = 1.0,
                             bool is_inverse = false,
                             bool is_arc = false,
                             bool is_mile = true,
                             const std::string& kernel = "",
                             double bandwidth = 0,
                             bool adaptive_bandwidth = false,
                             bool use_kernel_diagnals = false,
                             const std::string& polyid = "");
#endif

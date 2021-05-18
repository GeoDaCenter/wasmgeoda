/**
 * Author: Xun Li <lixun910@gmail.com>
 * Date: 2021-5-13
 *
 * Changes:
 * 2021-5-13
 */

#ifndef JSGEODA_JSGEODA_H
#define JSGEODA_JSGEODA_H



/**
 * CCentroids
 *
 * It is used to return centroids of geometries to js
 */
struct CCentroids {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> get_x() { return x;}
    std::vector<double> get_y() { return y;}
};

/**
 * WeightsResult
 *
 * It is used to return weights creation result to js
 */
struct WeightsResult {
    bool is_valid;
    int weight_type;
    int num_obs;
    bool is_symmetric;
    double sparsity;
    int max_nbrs;
    int min_nbrs;
    double median_nbrs;
    double mean_nbrs;
    std::string uid;

    bool get_is_valid() { return is_valid;}
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

void set_weights_content(GeoDaWeight* w, WeightsResult& rst);

WeightsResult queen_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold);

WeightsResult rook_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold);

WeightsResult knn_weights(std::string map_uid, int k, double power, bool is_inverse, bool is_arc, bool is_mile);

WeightsResult dist_weights(std::string map_uid, double dist_thres, double power, bool is_inverse, bool is_arc,
        bool is_mile);

WeightsResult kernel_weights(std::string map_uid, int k, std::string kernel, bool adaptive_bandwidth,
                             bool use_kernel_diagonals, bool is_arc, bool is_mile);

WeightsResult kernel_bandwidth_weights(std::string map_uid, double dist_thres, std::string kernel,
                                       bool use_kernel_diagonals, bool is_arc, bool is_mile);

double get_min_dist_threshold(std::string map_uid, bool is_arc, bool is_mile);

#endif //JSGEODA_JSGEODA_H

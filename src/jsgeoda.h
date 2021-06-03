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
    std::string map_uid;

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
    std::string get_map_uid() { return map_uid;}
};

/**
 *  Functions of weights
 */


void set_weights_content(GeoDaWeight* w, const string& map_uid, WeightsResult& rst);

WeightsResult queen_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold);

WeightsResult rook_weights(std::string map_uid, int order, int include_lower_order, double precision_threshold);

WeightsResult knn_weights(std::string map_uid, int k, double power, bool is_inverse, bool is_arc, bool is_mile);

WeightsResult dist_weights(std::string map_uid, double dist_thres, double power, bool is_inverse, bool is_arc,
        bool is_mile);

WeightsResult kernel_weights(std::string map_uid, int k, std::string kernel, bool adaptive_bandwidth,
                             bool use_kernel_diagonals, double power, bool is_inverse, bool is_arc, bool is_mile);

WeightsResult kernel_bandwidth_weights(std::string map_uid, double dist_thres, std::string kernel,
                                       bool use_kernel_diagonals, double power, bool is_inverse, bool is_arc,
                                       bool is_mile);

double get_min_dist_threshold(std::string map_uid, bool is_arc, bool is_mile);

/**
 *  Functions of mapping
 */

std::vector<double> natural_breaks(int k, const std::vector<double>& data,
                                    const std::vector<int>& undefs);

std::vector<double> quantile_breaks(int k, const std::vector<double>& data,
                                    const std::vector<int>& undefs);

std::vector<double> percentile_breaks(const std::vector<double>& data,
                                      const std::vector<int>& undefs);

std::vector<double> stddev_breaks(const std::vector<double>& data,
                                  const std::vector<int>& undefs);

std::vector<double> hinge15_breaks(const std::vector<double>& data,
                                   const std::vector<int>& undefs);

std::vector<double> hinge30_breaks(const std::vector<double>& data,
                                   const std::vector<int>& undefs);

std::vector<double> excess_risk(const std::vector<double>& event_data,
                                const std::vector<double>& base_data);

std::vector<double> eb_risk(const std::vector<double>& event_data,
                            const std::vector<double>& base_data);

std::vector<double> spatial_lag(const std::string map_uid,
                                const std::string& weight_uid,
                                const std::vector<double>& data, bool is_binary,
                                bool row_standardize, bool include_diagonal);

std::vector<double> spatial_rate(const std::vector<double>& event_data,
                                 const std::vector<double>& base_data,
                                 const std::string map_uid,
                                 const std::string& weight_uid);

std::vector<double> spatial_eb(const std::vector<double>& event_data,
                               const std::vector<double>& base_data,
                               const std::string map_uid,
                               const std::string& weight_uid);

/**
 * Functions of LISA
 */
// represent lisa results
struct LisaResult {
    bool is_valid;
    std::vector<double> sig_local_vec;
    std::vector<int> sig_cat_vec;
    std::vector<int> cluster_vec;
    std::vector<double> lag_vec;
    std::vector<double> lisa_vec;
    std::vector<int> nn_vec;
    std::vector<std::string> labels;
    std::vector<std::string> colors;

    bool get_is_valid() { return is_valid; }
    std::vector<double> get_sig_local() { return sig_local_vec;}
    std::vector<int> get_sig_cat() { return sig_cat_vec;}
    std::vector<int>  get_cluster() { return cluster_vec;}
    std::vector<double>  get_lag() { return lag_vec;}
    std::vector<double>  get_lisa() { return lisa_vec;}
    std::vector<int>  get_nn() { return nn_vec;}
    std::vector<std::string>  get_labels() { return labels;}
    std::vector<std::string>  get_colors() { return colors;}
};

class LISA;

void set_lisa_content(LISA* lisa, LisaResult& rst);

LisaResult local_moran(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used);

LisaResult local_g(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                   const std::vector<int>& undefs, double significance_cutoff, int permutations,
                   const std::string& permutation_method, int last_seed_used);

LisaResult local_gstar(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                   const std::vector<int>& undefs, double significance_cutoff, int permutations,
                   const std::string& permutation_method, int last_seed_used);

LisaResult local_geary(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used);

LisaResult local_joincount(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                       const std::vector<int>& undefs, double significance_cutoff, int permutations,
                       const std::string& permutation_method, int last_seed_used);

LisaResult quantile_lisa(const std::string map_uid, const std::string weight_uid, int k, int quantile,
                         const std::vector<double>& vals, const std::vector<int>& undefs, double significance_cutoff,
                         int permutations, const std::string& permutation_method, int last_seed_used);

std::vector<std::vector<double> > neighbor_match_test(const std::string map_uid,
                                                      int knn, double power, bool is_inverse, bool is_arc,
                                                      bool is_mile, const std::vector<std::vector<double> >& data,
                                                      const std::string& scale_method, const std::string& dist_type);

LisaResult local_moran_eb(const std::string map_uid, const std::string weight_uid, const std::vector<double>& vals,
                          const std::vector<double>& base_vals, double significance_cutoff, int permutations,
                          const std::string& permutation_method, int last_seed_used);

LisaResult multi_quantile_lisa(const std::string map_uid, const std::string weight_uid,
                               const std::vector<int> &k_s, const std::vector<int> &quantile_s,
                               const std::vector<std::vector<double> > &data,
                               const std::vector<std::vector<int> > &undefs,
                               double significance_cutoff, int permutations,
                               const std::string& permutation_method, int last_seed_used);

LisaResult local_multijoincount(const std::string map_uid, const std::string weight_uid,
                           const std::vector<std::vector<double> > &data,
                           const std::vector<std::vector<int> > &undefs, double significance_cutoff,
                           int permutations, const std::string& permutation_method, int last_seed_used);

LisaResult local_multigeary(const std::string map_uid, const std::string weight_uid,
                            const std::vector<std::vector<double> > &data,
                            const std::vector<std::vector<int> > &undefs, double significance_cutoff,
                            int permutations, const std::string& permutation_method, int last_seed_used);


/**
 * Functions of cartogram
 */

struct CartogramResult {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> r;

    std::vector<double> get_x() { return x;}
    std::vector<double> get_y() { return y;}
    std::vector<double> get_radius() { return r;}
};

#endif //JSGEODA_JSGEODA_H

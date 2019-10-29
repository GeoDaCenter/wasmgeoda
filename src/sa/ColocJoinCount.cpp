//
// Created by Xun Li on 9/27/19.
//

#include "ColocJoinCount.h"

ColocJoinCount::ColocJoinCount(int num_obs, GeoDaWeight *w, const std::vector<std::vector<int> > &_data,
                               const std::vector<std::vector<bool> > &_undefs)
: LISA(num_obs, w), data(_data), undefs(_undefs),
  CLUSTER_NOT_SIG(0),
  CLUSTER_SIG(1),
  CLUSTER_UNDEFINED(2),
  CLUSTER_NEIGHBORLESS(3)
{
    num_vars = data.size();
    z.resize(num_obs, 1);
    z_undefs.resize(num_obs, false);

    Run();
}

ColocJoinCount::~ColocJoinCount() {

}

bool ColocJoinCount::IsUndefsValid() {
    return undefs.size() == num_vars && undefs[0].size() == num_obs;
}

void ColocJoinCount::ComputeLoalSA() {
    // get colocations
    for (size_t i=0; i<num_obs; i++) {
        for (size_t v=0; v<num_vars; ++v) {
            if (IsUndefsValid())
                z_undefs[i] = z_undefs[i] && undefs[v][i];
        }
    }
    for (size_t i=0; i<num_obs; i++) {
        for (size_t v=0; v<num_vars; ++v) {
            if (z_undefs[i])
                z_undefs[i] = z_undefs[i] && undefs[v][i];
        }
    }
    for (size_t i=0; i<num_obs; i++) {
        if (z_undefs[i] == true) {
            lag_vec[i] = 0;
            lisa_vec[i] = 0;
            cluster_vec[i] = CLUSTER_UNDEFINED;
        } else {
            if (weights->GetNbrSize(i) == 0) {
                cluster_vec[i] = CLUSTER_NEIGHBORLESS;
            } else {

            }
        }
    }
}

void
ColocJoinCount::PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double> &permutedSA) {

}

uint64_t ColocJoinCount::CountLargerSA(int cnt, const std::vector<double> &permutedSA) {
    return 0;
}

std::vector<int> ColocJoinCount::GetClusterIndicators() {
    return std::vector<int>();
}

//
// Created by Xun Li on 9/27/19.
//

#ifndef GEODA_COLOCJOINCOUNT_H
#define GEODA_COLOCJOINCOUNT_H

#include <vector>
#include "LISA.h"

class ColocJoinCount : public LISA {
    const unsigned long CLUSTER_NOT_SIG;
    const unsigned long CLUSTER_SIG;
    const unsigned long CLUSTER_UNDEFINED;
    const unsigned long CLUSTER_NEIGHBORLESS;

public:
    ColocJoinCount(int num_obs,
                   GeoDaWeight* w,
                   const std::vector<std::vector<int> >& data,
                   const std::vector<std::vector<bool> >& undefs = std::vector<std::vector<bool> >());

    virtual ~ColocJoinCount();

    virtual void ComputeLoalSA() ;

    virtual void PermLocalSA(int cnt, int perm, const std::vector<int> &permNeighbors, std::vector<double>& permutedSA);

    virtual uint64_t CountLargerSA(int cnt, const std::vector<double>& permutedSA);

    virtual std::vector<int> GetClusterIndicators();

    const std::vector<double> &GetLagValues() const;

    const std::vector<double> &GetLocalGearyValues() const;

protected:
    int num_vars;
    std::vector<int> z;
    std::vector<bool > z_undefs;
    std::vector<std::vector<int> > data;
    std::vector<std::vector<bool> > undefs;

    bool IsUndefsValid();
};


#endif //GEODA_COLOCJOINCOUNT_H

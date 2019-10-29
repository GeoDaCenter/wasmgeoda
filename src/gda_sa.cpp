
#include "weights/GeodaWeight.h"
#include "sa/UniGeary.h"
#include "sa/UniG.h"
#include "sa/UniGstar.h"
#include "sa/UniJoinCount.h"
#include "sa/UniLocalMoran.h"
#include "gda_sa.h"

UniG* gda_localg(GeoDaWeight *w,
                 const std::vector<double> &data,
                 const std::vector<bool> &undefs)
{
    if (w == 0) return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty()) {
        copy_undefs.resize(num_obs, false);
    }
    UniG* localg = new UniG(num_obs, w, data, copy_undefs);
    return localg;
}

UniGstar* gda_localgstar(GeoDaWeight *w,
                         const std::vector<double> &data,
                         const std::vector<bool> &undefs)
{
    if (w == 0) return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty()) {
        copy_undefs.resize(num_obs, false);
    }
    UniGstar* localgstar = new UniGstar(num_obs, w, data, copy_undefs);
    return localgstar;
}

UniLocalMoran* gda_localmoran(GeoDaWeight *w,
                        const std::vector<double> &data,
                        const std::vector<bool> &undefs)
{
    if (w == 0) return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty()) {
        copy_undefs.resize(num_obs, false);
    }
    UniLocalMoran* lisa = new UniLocalMoran(num_obs, w, data, copy_undefs);
    return lisa;
}

UniGeary* gda_geary(GeoDaWeight *w,
                    const std::vector<double> &data,
                    const std::vector<bool> &undefs)
{
    if (w == 0) return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty()) {
        copy_undefs.resize(num_obs, false);
    }
    UniGeary* geary = new UniGeary(num_obs, w, data, copy_undefs);
    return geary;
}


UniJoinCount* gda_joincount(GeoDaWeight *w,
                            const std::vector<double> &data,
                            const std::vector<bool> &undefs)
{
    if (w == 0) return 0;

    int num_obs = w->num_obs;

    std::vector<bool> copy_undefs = undefs;
    if (copy_undefs.empty()) {
        copy_undefs.resize(num_obs, false);
    }
    UniJoinCount* jc= new UniJoinCount(num_obs, w, data, copy_undefs);
    return jc;
}
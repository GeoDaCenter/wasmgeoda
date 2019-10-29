#ifndef __JSGEODSA_GDA_SA__
#define __JSGEODSA_GDA_SA__

#include <string>
#include <vector>

class GeoDaWeight;
class UniLocalMoran;
class UniGeary;
class UniJoinCount;
class UniG;
class UniGstar;

// APIs of local spatial autocorrelation
/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
UniLocalMoran* gda_localmoran(GeoDaWeight *w,
                        const std::vector<double> &data,
                        const std::vector<bool> &undefs = std::vector<bool>());

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
UniGeary* gda_geary(GeoDaWeight *w, const std::vector<double> &data,
                    const std::vector<bool> &undefs = std::vector<bool>());

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
UniJoinCount* gda_joincount(GeoDaWeight *w, const std::vector<double> &data,
                            const std::vector<bool> &undefs = std::vector<bool>());

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
UniG* gda_localg(GeoDaWeight *w, const std::vector<double> &data,
                 const std::vector<bool> &undefs = std::vector<bool>());

/**
 *
 * @param w
 * @param data
 * @param undefs
 * @return
 */
UniGstar* gda_localgstar(GeoDaWeight *w, const std::vector<double> &data,
                         const std::vector<bool> &undefs = std::vector<bool>());

#endif
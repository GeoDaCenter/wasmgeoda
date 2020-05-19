//
// Created by Xun Li on 2019-06-05.
//

#include <math.h>

#include "../GeoDaSet.h"
#include "../GenUtils.h"
#include "../weights/GeodaWeight.h"

#include "LISA.h"

std::map<std::string, bool>  LISA::has_cached_perm;
std::map<std::string, std::vector<std::vector<int> > >  LISA::cached_perm_nbrs;

LISA::LISA(int num_obs, GeoDaWeight* w)
: last_seed_used(123456789),
  permutations(999), reuse_last_seed(true),
  calc_significances(true),row_standardize(true),
  has_undefined(false), has_isolates(w->HasIsolates()),
  user_sig_cutoff(0), weights(w), num_obs(num_obs)
{
    if (has_cached_perm.find(w->GetUID()) == has_cached_perm.end()) {
        has_cached_perm[w->GetUID()] = false;
        cached_perm_nbrs[w->GetUID()] = std::vector<std::vector<int> >(0);
    }
    SetSignificanceFilter(1);
}

LISA::~LISA()
{
}

void LISA::Run()
{
    sig_local_vec.resize(num_obs, 0);
    sig_cat_vec.resize(num_obs, 0);
    cluster_vec.resize(num_obs, 0);
    lag_vec.resize(num_obs, 0);
    lisa_vec.resize(num_obs, 0);
    nn_vec.resize(num_obs, 0);

    for (size_t i=0; i<num_obs; i++) {
        nn_vec[i] = weights->GetNbrSize(i);
    }

    ComputeLoalSA();
    if (calc_significances) {
        CalcPseudoP();
    }
}

void LISA::SetSignificanceFilter(int filter_id)
{
    if (filter_id == -1) {
        // user input cutoff
        significance_filter = filter_id;
        return;
    }
    // 0: >0.05 1: 0.05, 2: 0.01, 3: 0.001, 4: 0.0001
    if (filter_id < 1 || filter_id > 4) return;
    significance_filter = filter_id;
    if (filter_id == 1) significance_cutoff = 0.05;
    if (filter_id == 2) significance_cutoff = 0.01;
    if (filter_id == 3) significance_cutoff = 0.001;
    if (filter_id == 4) significance_cutoff = 0.0001;
}

int LISA::GetSignificanceFilter()
{
    return significance_filter;
}

double LISA::GetSignificanceCutoff()
{
    return significance_cutoff;
}

void LISA::SetSignificanceCutoff(double val)
{
    significance_cutoff = val;
}

double LISA::GetUserCutoff()
{
    return user_sig_cutoff;
}

void LISA::SetUserCutoff(double val)
{
    user_sig_cutoff = val;
}

double LISA::GetFDR()
{
    return fdr;
}
void LISA::SetFDR(double val)
{
    fdr = val;
}

double LISA::GetBO()
{
    return bo;
}
void LISA::SetBO(double val)
{
    bo = val;
}

int LISA::GetNumPermutations()
{
    return permutations;
}
void LISA::SetNumPermutations(int val)
{
    permutations = val;
}

uint64_t LISA::GetLastUsedSeed()
{
    return last_seed_used;
}

bool LISA::IsReuseLastSeed()
{
    return reuse_last_seed;
}

void LISA::SetReuseLastSeed(bool reuse)
{
    reuse_last_seed = reuse;
}

void LISA::SetLastUsedSeed(uint64_t seed)
{
    reuse_last_seed = true;
    last_seed_used = seed;
}

bool LISA::GetHasIsolates()
{
    return has_isolates;
}

bool LISA::GetHasUndefined()
{
    return has_undefined;

}

void LISA::CalcPseudoP()
{
    if (!calc_significances) return;

    CalcPseudoP_range(0, num_obs-1, last_seed_used);

    if (has_cached_perm[weights->GetUID()] == false) {
        has_cached_perm[weights->GetUID()] = true;
    }
}

void LISA::CalcPseudoP_threaded()
{

}

void LISA::CalcPseudoP_range(int obs_start, int obs_end, uint64_t seed_start)
{
    GeoDaSet workPermutation(num_obs);
    int max_rand = num_obs-1;

    std::string wuid = weights->GetUID();
    bool using_cache = has_cached_perm[wuid];
    std::vector<std::vector<int> >& cache = cached_perm_nbrs[wuid];

    for (size_t cnt=obs_start; cnt<=obs_end; cnt++) {

        // get full neighbors even if has undefined value
        int numNeighbors = weights->GetNbrSize(cnt);
        if (numNeighbors == 0) {
            sig_cat_vec[cnt] = 5; // neighborless cat
            // isolate: don't do permutation
            continue;
        }


        std::vector<double> permutedSA(permutations, 0);
        if (using_cache == false) {
            for (size_t perm = 0; perm < permutations; perm++) {
                int rand = 0, newRandom;
                double rng_val;
                while (rand < numNeighbors) {
                    // computing 'perfect' permutation of given size
                    rng_val = Gda::ThomasWangHashDouble(seed_start++) * max_rand;
                    // round is needed to fix issue
                    // https://github.com/GeoDaCenter/geoda/issues/488
                    newRandom = (int) (rng_val < 0.0 ? ceil(rng_val - 0.5) : floor(rng_val + 0.5));

                    if (newRandom != cnt && !workPermutation.Belongs(newRandom) && weights->GetNbrSize(newRandom) > 0) {
                        workPermutation.Push(newRandom);
                        rand++;
                    }
                }
                std::vector<int> permNeighbors(numNeighbors);
                for (int cp = 0; cp < numNeighbors; cp++) {
                    permNeighbors[cp] = workPermutation.Pop();
                }
                cache.push_back(permNeighbors);
                PermLocalSA(cnt, perm, permNeighbors, permutedSA);

            }
        } else {
            for (size_t perm = 0; perm < permutations; perm++) {
                PermLocalSA(cnt, perm, cache[perm], permutedSA);
            }
        }

        uint64_t countLarger = CountLargerSA(cnt, permutedSA);
        double _sigLocal = (countLarger+1.0)/(permutations+1);

        // 'significance' of local Moran
        if (_sigLocal <= 0.0001) sig_cat_vec[cnt] = 4;
        else if (_sigLocal <= 0.001) sig_cat_vec[cnt] = 3;
        else if (_sigLocal <= 0.01) sig_cat_vec[cnt] = 2;
        else if (_sigLocal <= 0.05) sig_cat_vec[cnt] = 1;
        else sig_cat_vec[cnt] = 0;

        sig_local_vec[cnt] = _sigLocal;
        // observations with no neighbors get marked as isolates
        // NOTE: undefined should be marked as well, however, since undefined_cat has covered undefined category, we don't need to handle here
    }
}


std::vector<std::string> LISA::GetDefaultCategories()
{
    std::vector<std::string> cats;
    cats.push_back("p = 0.05");
    cats.push_back("p = 0.01");
    cats.push_back("p = 0.001");
    cats.push_back("p = 0.0001");
    return cats;
}

std::vector<double> LISA::GetDefaultCutoffs()
{
    std::vector<double> cutoffs;
    cutoffs.push_back(0.05);
    cutoffs.push_back(0.01);
    cutoffs.push_back(0.001);
    cutoffs.push_back(0.0001);
    return cutoffs;
}

std::vector<double> LISA::GetLocalSignificanceValues()
{
    return sig_local_vec;
}

std::vector<int> LISA::GetClusterIndicators()
{
    return cluster_vec;
}

std::vector<int> LISA::GetSigCatIndicators()
{
    return sig_cat_vec;
}

std::vector<int> LISA::GetNumNeighbors()
{
    return nn_vec;
}

std::vector<double> LISA::GetSpatialLagValues()
{
    return lag_vec;
}

std::vector<double> LISA::GetLISAValues()
{
    return lisa_vec;
}

bool LISA::IsRowStandardize() const {
    return row_standardize;
}

void LISA::SetRowStandardize(bool rowStandardize) {
    row_standardize = rowStandardize;
}

int LISA::GetNumThreads() const {
    return nCPUs;
}

void LISA::SetNumThreads(int n_threads) {
    nCPUs = n_threads;
}

std::vector<std::string> LISA::GetLabels()
{
    return labels;
}

std::vector<std::string> LISA::GetColors()
{
    return colors;
}

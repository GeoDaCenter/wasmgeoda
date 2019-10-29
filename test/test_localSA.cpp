//
// Created by Xun Li on 2019-06-06.
//

#include <vector>
#include <limits.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../src/gda_weights.h"
#include "../src/gda_sa.h"
#include "../src/geojson.h"
#include "../src/weights/GeodaWeight.h"
#include "../src/sa/UniLocalMoran.h"
#include "../src/sa/UniGeary.h"
#include "../src/sa/UniGstar.h"
#include "../src/sa/UniG.h"
#include "../src/sa/UniJoinCount.h"


using namespace testing;

namespace {

    // Deprecated in >= 0.0.2
    TEST(LOCALSA_TEST, LISA_UNI) {
        std::string file_path = "../data/Guerry.geojson";
        GdaGeojson gda(file_path);
        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<double> data = gda.GetNumericCol("Crm_prp");

        UniLocalMoran* lisa = gda_localmoran(w, data);

        std::vector<int> cvals= lisa->GetClusterIndicators();
        std::vector<double> pvals = lisa->GetLocalSignificanceValues();
        std::vector<double> mvals = lisa->GetLISAValues();
        delete lisa;

        EXPECT_DOUBLE_EQ(mvals[0], 0.015431978309803657);
        EXPECT_DOUBLE_EQ(mvals[1], 0.32706332236560332);
        EXPECT_DOUBLE_EQ(mvals[2], 0.021295296214118884);

        EXPECT_THAT(cvals[0], 0);
        EXPECT_THAT(cvals[1], 0);
        EXPECT_THAT(cvals[2], 1);

        EXPECT_DOUBLE_EQ(pvals[0], 0.41399999999999998);
        EXPECT_DOUBLE_EQ(pvals[1], 0.123);
        EXPECT_DOUBLE_EQ(pvals[2], 0.001);
    }

    TEST(LOCALSA_TEST, GEARY_UNI) {
        std::string file_path = "../data/Guerry.geojson";
        GdaGeojson gda(file_path);

        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<double> data = gda.GetNumericCol("Crm_prp");

        UniGeary* geary = gda_geary(w, data);

        std::vector<int> cvals = geary->GetClusterIndicators();
        std::vector<double> pvals = geary->GetLocalSignificanceValues();
        std::vector<double> gvals = geary->GetLISAValues();
        delete geary;

        EXPECT_DOUBLE_EQ(gvals[0], 7.3980833011783602);
        EXPECT_DOUBLE_EQ(gvals[1], 0.28361195650519017);
        EXPECT_DOUBLE_EQ(gvals[2], 3.6988922226329906);

        EXPECT_THAT(cvals[0], 0);
        EXPECT_THAT(cvals[1], 2);
        EXPECT_THAT(cvals[2], 4);

        EXPECT_DOUBLE_EQ(pvals[0], 0.39800000000000002);
        EXPECT_DOUBLE_EQ(pvals[1], 0.027);
        EXPECT_DOUBLE_EQ(pvals[2], 0.025000000000000001);
    }

    TEST(LOCALSA_TEST, JOINCOUNT_UNI) {
        std::string file_path = "../data/Guerry.geojson";
        GdaGeojson gda(file_path);

        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<double> data = gda.GetNumericCol("nsa");

        UniJoinCount* jc = gda_joincount(w, data);

        std::vector<int> nnvals = jc->GetNumNeighbors();
        std::vector<double> pvals = jc->GetLocalSignificanceValues();
        std::vector<double> jvals = jc->GetLISAValues();
        delete jc;

        EXPECT_DOUBLE_EQ(jvals[0], 2);
        EXPECT_DOUBLE_EQ(jvals[1], 3);
        EXPECT_DOUBLE_EQ(jvals[2], 4);

        EXPECT_THAT(nnvals[0], 2);
        EXPECT_THAT(nnvals[1], 3);
        EXPECT_THAT(nnvals[2], 4);

        EXPECT_DOUBLE_EQ(pvals[0], 0.21299999999999999);
        EXPECT_DOUBLE_EQ(pvals[1], 0.070000000000000007);
        EXPECT_DOUBLE_EQ(pvals[2], 0.017000000000000001);
    }

    TEST(LOCALSA_TEST, LOCALG_UNI) {
        std::string file_path = "../data/Guerry.geojson";
        GdaGeojson gda(file_path);

        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<double> data = gda.GetNumericCol("Crm_prp");

        UniG* localg = gda_localg(w, data);

        std::vector<int> cvals = localg->GetClusterIndicators();
        std::vector<double> pvals = localg->GetLocalSignificanceValues();
        std::vector<double> gvals = localg->GetLISAValues();
        delete localg;

        EXPECT_DOUBLE_EQ(gvals[0], 0.012077920687925825);
        EXPECT_DOUBLE_EQ(gvals[1], 0.0099240961298508561);
        EXPECT_DOUBLE_EQ(gvals[2], 0.018753584525825453);

        EXPECT_THAT(cvals[0], 0);
        EXPECT_THAT(cvals[1], 0);
        EXPECT_THAT(cvals[2], 1);

        EXPECT_DOUBLE_EQ(pvals[0], 0.414);
        EXPECT_DOUBLE_EQ(pvals[1], 0.123);
        EXPECT_DOUBLE_EQ(pvals[2], 0.001);
    }

    TEST(LOCALSA_TEST, LOCALGstar_UNI) {
        std::string file_path = "../data/Guerry.geojson";
        GdaGeojson gda(file_path);

        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<double> data = gda.GetNumericCol("Crm_prp");

        UniGstar* localgstar = gda_localgstar(w, data);

        std::vector<int> cvals = localgstar->GetClusterIndicators();
        std::vector<double> pvals = localgstar->GetLocalSignificanceValues();
        std::vector<double> gvals = localgstar->GetLISAValues();
        delete localgstar;

        EXPECT_DOUBLE_EQ(gvals[0], 0.014177043620524426);
        EXPECT_DOUBLE_EQ(gvals[1], 0.0096136007223101994);
        EXPECT_DOUBLE_EQ(gvals[2], 0.017574324039034434);

        EXPECT_THAT(cvals[0], 0);
        EXPECT_THAT(cvals[1], 0);
        EXPECT_THAT(cvals[2], 1);

        EXPECT_DOUBLE_EQ(pvals[0], 0.414);
        EXPECT_DOUBLE_EQ(pvals[1], 0.123);
        EXPECT_DOUBLE_EQ(pvals[2], 0.001);
    }
}
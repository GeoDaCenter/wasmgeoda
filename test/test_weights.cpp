//
// Created by Xun Li on 2019-06-04.
//

#include <string>
#include <limits.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../src/gda_weights.h"
#include "../src/geojson.h"

using namespace testing;

namespace {

    TEST(WEIGHTS_TEST, QUEEN_CREATE) {
        std::string file_path = "../data/natregimes.geojson";

        GdaGeojson gda(file_path);
        GeoDaWeight* w = gda_queen_weights(&gda);

        EXPECT_THAT(w->num_obs, 3085);
        EXPECT_THAT(w->GetMinNumNbrs(), 1);
        EXPECT_THAT(w->GetMaxNumNbrs(), 14);
        EXPECT_TRUE(w->is_symmetric);
        EXPECT_DOUBLE_EQ(w->GetSparsity(), 0);
        EXPECT_DOUBLE_EQ(w->GetDensity(), 0.19089598070866245);
        EXPECT_DOUBLE_EQ(w->GetMeanNumNbrs(), 5.8891410048622364);

        delete w;
    }

    TEST(WEIGHTS_TEST, ROOK_CREATE) {
        std::string file_path = "../data/natregimes.geojson";

        GdaGeojson gda(file_path);
        GeoDaWeight* w = gda_rook_weights(&gda);

        EXPECT_THAT(w->num_obs, 3085);
        EXPECT_THAT(w->GetMinNumNbrs(), 1);
        EXPECT_THAT(w->GetMaxNumNbrs(), 13);
        EXPECT_TRUE(w->is_symmetric);
        EXPECT_DOUBLE_EQ(w->GetSparsity(), 0);
        EXPECT_DOUBLE_EQ(w->GetDensity(), 0.18059886153789576);
        EXPECT_DOUBLE_EQ(w->GetMeanNumNbrs(), 5.571474878444084);

        delete w;
    }

    TEST(WEIGHTS_TEST, KNN_CREATE) {
        std::string file_path = "../data/natregimes.geojson";

        GdaGeojson gda(file_path);

        GeoDaWeight* w = gda_knn_weights(&gda, 4);

        EXPECT_FALSE(w->is_symmetric);
        EXPECT_THAT(w->num_obs, 3085);
        EXPECT_THAT(w->GetMinNumNbrs(), 4);
        EXPECT_THAT(w->GetMaxNumNbrs(), 4);
        EXPECT_DOUBLE_EQ(w->GetSparsity(), 0);
        EXPECT_DOUBLE_EQ(w->GetDensity(), 0.12965964343598055);
        EXPECT_DOUBLE_EQ(w->GetMeanNumNbrs(), 4);

        delete w;
    }

    /*
    TEST(WEIGHTS_TEST, DIST_CREATE) {
        GeoDa gda("../data/natregimes.shp");
        double min_thres = gda_min_distthreshold(&gda);

        EXPECT_DOUBLE_EQ(min_thres, 1.4657759325950015);

        GeoDaWeight* w = gda_distance_weights(&gda, min_thres);

        EXPECT_FALSE(w->is_symmetric);
        EXPECT_THAT(w->num_obs, 3085);
        EXPECT_THAT(w->GetMinNumNbrs(), 1);
        EXPECT_THAT(w->GetMaxNumNbrs(), 85);
        EXPECT_DOUBLE_EQ(w->GetSparsity(), 0);
        EXPECT_DOUBLE_EQ(w->GetDensity(), 1.1939614751148575);
        EXPECT_DOUBLE_EQ(w->GetMeanNumNbrs(), 36.833711507293351);

        delete w;
    }

    TEST(WEIGHTS_TEST, KERNEL_KNN) {
        GeoDa gda("../data/natregimes.shp");

        double power = 1;
        bool is_inverse = false;
        bool is_arc = false;
        bool is_mile = true;
        std::string kernel = "triangular";
        double bandwidth  = 0;
        bool adaptive_bandwidth = true;
        bool use_kernel_diagonals = false;
        int k = 15;
        GeoDaWeight* w = gda_knn_weights(&gda, k, power, is_inverse,
                is_arc, is_mile, kernel, bandwidth,
                adaptive_bandwidth, use_kernel_diagonals);

        EXPECT_FALSE(w->is_symmetric);
        EXPECT_THAT(w->num_obs, 3085);
        EXPECT_THAT(w->GetMinNumNbrs(), 15);
        EXPECT_THAT(w->GetMaxNumNbrs(), 15);
        EXPECT_DOUBLE_EQ(w->GetSparsity(), 0);
        EXPECT_DOUBLE_EQ(w->GetDensity(), 0.48622366288492708);
        EXPECT_DOUBLE_EQ(w->GetMeanNumNbrs(), 15);

        delete w;
    }
     */
}
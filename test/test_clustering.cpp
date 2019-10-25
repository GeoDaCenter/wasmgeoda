//
// Created by Xun Li on 2019-06-06.
//

#include <vector>
#include <limits.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <../src/gda_weights.h>

using namespace testing;

namespace {

    const char *col_names[6] = {"Crm_prs", "Crm_prp", "Litercy", "Donatns", "Infants", "Suicids"};

    TEST(CLUSTERING_TEST, SKATER) {
        GeoDa gda("../data/Guerry.shp");
        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<std::vector<double> > data;
        for (size_t i=0; i<6; ++i) {
            data.push_back( gda.GetNumericCol(col_names[i]) );
        }
        std::vector<std::vector<int> > clst = gda_skater(4, w, data);
        double totalss = gda_totalsumofsquare(data);
        double withinss = gda_withinsumofsquare(clst, data);
        double ratio = (totalss - withinss) / totalss;
        delete w;

        //EXPECT_THAT(clst, ElementsAre(3,2,3,1,1,1,2,1,2,1,1,1,2,1,1,3,3,3,2,4,3,1,2,1,2,2,4,1,1,1,1,1,4,3,4,1,2,1,4,
        //                              3,3,4,2,1,1,1,4,4,2,2,4,2,2,4,2,3,2,2,4,2,3,1,1,1,2,2,1,2,3,4,2,2,2,2,3,2,1,1,
        //1,1,3,3,3,2,2));
        EXPECT_DOUBLE_EQ(ratio, 0.31564466593112039);
    }

    TEST(CLUSTERING_TEST, REDCAP) {
        GeoDa gda("../data/Guerry.shp");
        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<std::vector<double> > data;
        for (size_t i=0; i<6; ++i) {
            data.push_back( gda.GetNumericCol(col_names[i]) );
        }
        std::string method = "firstorder-singlelinkage";
        std::vector<std::vector<int> > clst = gda_redcap(4, w, data, method);
        double totalss = gda_totalsumofsquare(data);
        double withinss = gda_withinsumofsquare(clst, data);
        double ratio = (totalss - withinss) / totalss;
        delete w;

        // EXPECT_THAT(clst, ElementsAre(3,2,3,1,1,1,2,1,2,1,1,1,2,1,1,3,3,3,2,4,3,1,2,1,2,2,4,1,1,1,1,1,4,3,4,1,2,1,4,
        //        3,3,4,2,1,1,1,4,4,2,2,4,2,2,4,2,3,2,2,4,2,3,1,1,1,2,2,1,2,3,4,2,2,2,2,3,2,1,1,1,1,3,3,3,2,2));
        EXPECT_DOUBLE_EQ(ratio, 0.31564466593112039);
    }

    TEST(CLUSTERING_TEST, MAXP_GREEDY) {
        GeoDa gda("../data/Guerry.shp");
        GeoDaWeight* w = gda_queen_weights(&gda);
        std::vector<std::vector<double> > data;
        for (size_t i=0; i<6; ++i) {
            data.push_back( gda.GetNumericCol(col_names[i]) );
        }
        std::vector<double> bound_vals = gda.GetNumericCol("Pop1831");
        double min_bound = 3236.6700000000001; // 10% of Pop1831

        int initial = 99;
        std::vector<std::vector<int> > clst = gda_maxp(w, data, bound_vals, min_bound, "greedy", initial);
        double totalss = gda_totalsumofsquare(data);
        double withinss = gda_withinsumofsquare(clst, data);
        double ratio = (totalss - withinss) / totalss;

        delete w;
        EXPECT_DOUBLE_EQ(ratio, 0.50701807973320201);

    }
}
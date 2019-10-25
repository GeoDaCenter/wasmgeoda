//
// Created by Xun Li on 2019-06-13.
//

#include <limits.h>
#include <libgeoda.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;

namespace {

    TEST(TableFunctions, intCol) {
        GeoDa gda("../data/natregimes.shp");
        std::vector<long long> vals = gda.GetIntegerCol("POLY_ID");

        EXPECT_THAT(gda.GetNumObs(), 3085);
    }
}
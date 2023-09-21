#include "gtest/gtest.h"
#define EEOS_POLYTROPE
#include "pkd.h"
#include "core/simd.h"
#include "eEOS/eEOS_struct.h"
#include "eEOS/eEOS.h"

template<typename dtype, typename mtype>
class eEOSTestBase : public ::testing::Test {
protected:
    struct eEOSparam params;
    void SetUp() override {
        // These are not meant to be physical, just to check the implementation
        params.dPolyFloorMinOD = 57.0;
        params.dPolyFloorExponent = 1.0;
        params.dPolyFloorDen = 70.0;
        params.dPolyFlooru = 1.0;

        params.dFlooru = 1.0;
        params.dFloorDen = 1.0;
        params.dFloorMinOD = 57.0;
    }
    void TearDown() override {}
    dtype gamma = 5/3.;
    dtype ball = 1.0;
    dtype a_inv3 = 1.0;
};
typedef vec<double,double> dtype;
typedef eEOSTestBase<vec<double,double>,mmask<bool>> eEOSTest;

TEST_F(eEOSTest, LowDensity) {
    // The eEOS should not activate in this case
    dtype dens = 0.1;
    dtype eos = eEOSEnergyFloor(a_inv3, dens, ball, gamma, params);
    EXPECT_EQ(eos, NOT_IN_EEOS);
}
TEST_F(eEOSTest, MidDensity) {
    // The eEOS should not activate in this case
    dtype dens = 10.0;
    dtype eos = eEOSEnergyFloor(a_inv3, dens, ball, gamma, params);
    EXPECT_EQ(eos, NOT_IN_EEOS);
}
TEST_F(eEOSTest, HighDensity) {
    // The eEOS should not activate in this case
    dtype dens = 60.0;
    dtype eos = eEOSEnergyFloor(a_inv3, dens, ball, gamma, params);
    EXPECT_EQ(eos, params.dFlooru);
}
TEST_F(eEOSTest, VeryHighDensity) {
    // The constant energy eEOS should be active
    dtype dens = 80.0;
    dtype eos = eEOSEnergyFloor(a_inv3, dens, ball, gamma, params);
    EXPECT_EQ(eos, pow(dens/params.dPolyFloorDen, params.dPolyFloorExponent) );
}

#include "gtest/gtest.h"

#include <cstdlib>
#include "imf.h"


class IMFTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};


TEST_F(IMFTest, BasicIMFTest) {
    const double dMinMass{0.1};
    const double dMaxMass{100.0};

    // Try to load an unexisting IMF
    EXPECT_DEATH({
        auto wrongIMF = ChooseIMF("other", dMinMass, dMaxMass);
    }, "Undefined IMF type has been given to ChooseIMF: other");

    auto ChabrierIMF = ChooseIMF("chabrier", dMinMass, dMaxMass);
    auto KroupaIMF = ChooseIMF("kroupa", dMinMass, dMaxMass);
    auto SalpeterIMF = ChooseIMF("salpeter", dMinMass, dMaxMass);

    // Check that the IMFs are correctly normalized after initialization
    double norm = ChabrierIMF->MassWeightedIntegration(dMinMass, dMaxMass);
    EXPECT_NEAR(norm, 1.0, 1e-13);

    norm = KroupaIMF->MassWeightedIntegration(dMinMass, dMaxMass);
    EXPECT_NEAR(norm, 1.0, 1e-13);

    norm = SalpeterIMF->MassWeightedIntegration(dMinMass, dMaxMass);
    EXPECT_NEAR(norm, 1.0, 1e-13);

    // Compare with results from Dalla Vecchia & Schaye 2012, Section 2
    double test1 = ChabrierIMF->UnweightedIntegration(6.0, dMaxMass);
    EXPECT_NEAR(test1, 0.01736, 1e-5);
    double test2 = ChabrierIMF->UnweightedIntegration(8.0, dMaxMass);
    EXPECT_NEAR(test2, 0.01180, 1e-5);

    test1 = KroupaIMF->UnweightedIntegration(6.0, dMaxMass);
    EXPECT_NEAR(test1, 0.01736, 3e-3);
    test2 = KroupaIMF->UnweightedIntegration(8.0, dMaxMass);
    EXPECT_NEAR(test2, 0.01180, 3e-3);

    test1 = SalpeterIMF->UnweightedIntegration(6.0, dMaxMass);
    EXPECT_NEAR(test1, 0.01107, 1e-5);
    test2 = SalpeterIMF->UnweightedIntegration(8.0, dMaxMass);
    EXPECT_NEAR(test2, 0.00742, 1e-5);
}

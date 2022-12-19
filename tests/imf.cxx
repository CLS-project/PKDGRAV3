#include "gtest/gtest.h"

#include <cstdlib>
#include "imf.h"


class IMFTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};


TEST_F(IMFTest, ChabrierIMFTest) {
    // Try to load an unexisting IMF
    EXPECT_DEATH({
        auto wrongIMF = ChooseIMF("other", 0.1, 100.0);
    }, "Undefined IMF type has been given to ChooseIMF: other");

    auto IMF = ChooseIMF("chabrier", 0.1, 100.0);

    // Check that the IMF is correctly normalized after initialization
    double norm = IMF->MassWeightedIntegration(0.1, 100.0);
    EXPECT_EQ(norm, 1.);

    // Compare with results from Dalla Vecchia & Schaye 2012, Section 2
    double test1 = IMF->UnweightedIntegration(6.,100.0);
    EXPECT_NEAR(test1, 0.01736, 1e-5);
    double test2 = IMF->UnweightedIntegration(8.,100.0);
    EXPECT_NEAR(test2, 0.01180, 1e-5);
}

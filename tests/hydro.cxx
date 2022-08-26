#include "gtest/gtest.h"
#include <vector>

#include "smooth/smoothfcn.h"
#include "hydro/riemann_own.h"

const double TOL = 5e-6;


class RiemannTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}


    SMF smf;
    std::vector<double> n_unit {1., 0., 0.};
    double rho_R, p_R;
    double rho_L, p_L;
    std::vector<double> v_R;
    std::vector<double> v_L;
    double P_M, S_M, rho_f, p_f, v_f[3];
    double dConstGamma = 1.4;
    void set_R(double rho, double p, std::vector<double> v) {
        rho_R = rho;
        p_R   = p;
        v_R   = v;
    }
    void set_L(double rho, double p, std::vector<double> v) {
        rho_L = rho;
        p_L   = p;
        v_L   = v;
    }

    int solve() {
        return Riemann_solver_exact(dConstGamma, rho_R, p_R, v_R.data(),
                                    rho_L, p_L, v_L.data(),
                                    &P_M, &S_M,
                                    &rho_f, &p_f, &v_f[0],
                                    n_unit.data());

    }


};




TEST_F(RiemannTest, ToroTest1) {
    int iters;

    // Test 1
    set_R(0.125, 0.1, {0.,0.,0.});
    set_L(1., 1., {0.,0.,0.});
    iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 0.30313, 0.30313*TOL);
    EXPECT_NEAR(S_M, 0.92745, 0.92745*TOL);
}

TEST_F(RiemannTest, ToroTest2) {
    int iters;
    // Test 2
    set_R(1., 0.4, {2.,0.,0.});
    set_L(1., 0.4, {-2.,0.,0.});
    iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 0.00189, 0.000005);
    EXPECT_NEAR(S_M, 0.0, 5e-6);
}

TEST_F(RiemannTest, ToroTest3) {
    int iters;
    // Test 3
    set_R(1., 0.01, {0.,0.,0.});
    set_L(1., 1000., {0.,0.,0.});
    iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 460.894, 460.894*TOL);
    EXPECT_NEAR(S_M, 19.5975, 19.5975*TOL);
}

TEST_F(RiemannTest, ToroTest4) {
    int iters;
    // Test 4
    set_R(1., 100., {0.,0.,0.});
    set_L(1., 0.01, {0.,0.,0.});
    iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M,  46.0950, 46.0950*TOL);
    EXPECT_NEAR(S_M, -6.19633, 6.19633*TOL);
}

TEST_F(RiemannTest, ToroTest5) {
    int iters;
    // Test 5
    set_R(5.99242, 46.0950, {-6.19633,0.,0.});
    set_L(5.99924, 460.894, {19.5975,0.,0.});
    iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 1691.64, 1691.64*TOL);
    EXPECT_NEAR(S_M, 8.68975, 8.68975*TOL);
}


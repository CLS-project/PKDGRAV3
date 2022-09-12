#include "gtest/gtest.h"
#define USE_MFM
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
    dvec vn_unit[3] = {1., 0., 0.};
    double rho_R, p_R;
    double rho_L, p_L;
    std::vector<double> v_R;
    std::vector<double> v_L;
    double P_M, S_M, rho_f, p_f, v_f[3];
    double dConstGamma = 1.4;
    dvec vrho_R, vp_R;
    dvec vrho_L, vp_L;
    dvec vv_R[3];
    dvec vv_L[3];
    dvec vP_M, vS_M, vrho_f, vp_f, vv_f[3];
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
    void set_vec_R(double rho, double p, std::vector<double> v) {
        vrho_R = rho;
        vp_R   = p;
        for (auto k=0; k<3; k++)
            vv_R[k]   = v[k];
    }
    void set_vec_L(double rho, double p, std::vector<double> v) {
        vrho_L = rho;
        vp_L   = p;
        for (auto k=0; k<3; k++)
            vv_L[k]   = v[k];
    }

    int solve() {
        return Riemann_solver_exact(dConstGamma, rho_R, p_R, v_R.data(),
                                    rho_L, p_L, v_L.data(),
                                    P_M, S_M,
                                    &rho_f, &p_f, &v_f[0],
                                    n_unit.data());
    }

    dvec solve_vec() {
        dvec gamma = 1.4;
        return Riemann_solver_exact(gamma, vrho_R, vp_R, vv_R,
                                    vrho_L, vp_L, vv_L,
                                    vP_M, vS_M,
                                    &vrho_f, &vp_f, &vv_f[0],
                                    vn_unit);
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
TEST_F(RiemannTest, ToroVec1) {

    set_vec_R(0.125, 0.1, {0.,0.,0.});
    set_vec_L(1., 1., {0.,0.,0.});
    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 0.30313, 0.30313*TOL);
        EXPECT_NEAR(vS_M[i], 0.92745, 0.92745*TOL);
    }
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
TEST_F(RiemannTest, ToroVec2) {

    set_vec_R(1., 0.4, {2.,0.,0.});
    set_vec_L(1., 0.4, {-2.,0.,0.});
    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 0.00189, 0.000005);
        EXPECT_NEAR(vS_M[i], 0.0, 5e-6);
    }
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
TEST_F(RiemannTest, ToroVec3) {

    set_vec_R(1., 0.01, {0.,0.,0.});
    set_vec_L(1., 1000., {0.,0.,0.});
    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 460.894, 460.894*TOL);
        EXPECT_NEAR(vS_M[i], 19.5975, 19.5975*TOL);
    }
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
TEST_F(RiemannTest, ToroVec4) {

    set_vec_R(1., 100., {0.,0.,0.});
    set_vec_L(1., 0.01, {0.,0.,0.});
    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i],  46.0950, 46.0950*TOL);
        EXPECT_NEAR(vS_M[i], -6.19633, 6.19633*TOL);
    }
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

TEST_F(RiemannTest, ToroVec5) {

    set_vec_R(5.99242, 46.0950, {-6.19633,0.,0.});
    set_vec_L(5.99924, 460.894, {19.5975,0.,0.});
    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 1691.64, 1691.64*TOL);
        EXPECT_NEAR(vS_M[i], 8.68975, 8.68975*TOL);
    }
}

TEST_F(RiemannTest, ToroVecAll) {

    double allrho_R[4] = {0.125, 1.,   1.,   5.99242};
    double allrho_L[4] = {1.,    1.,   1.,   5.99924};
    double allp_R[4] =   {0.1,   0.4, 0.01,  46.0950};
    double allp_L[4] =   {1.,    0.4, 1000., 460.894};
    double allv_R[4] =   {0.,    2.,   0.,  -6.19633};
    double allv_L[4] =   {0.,   -2.,   0.,  19.5975};

    vrho_R.load(allrho_R);
    vrho_L.load(allrho_L);
    vp_R.load(allp_R);
    vp_L.load(allp_L);
    vv_R[0].load(allv_R);
    vv_L[0].load(allv_L);
    vv_R[1] = 0.0;
    vv_L[1] = 0.0;
    vv_R[2] = 0.0;
    vv_L[2] = 0.0;

    dvec iters = solve_vec();

    for (auto i=0; i<SIMD_DWIDTH; i++) {
        EXPECT_NE(iters[i], 0);
        //printf("%e %e %e\n", iters[i], vP_M[i], vS_M[i]);
    }
    //Test 1
    EXPECT_NEAR(vP_M[0], 0.30313, 0.30313*TOL);
    EXPECT_NEAR(vS_M[0], 0.92745, 0.92745*TOL);
    //Test 2
    EXPECT_NEAR(vP_M[1], 0.00189, 0.000005);
    EXPECT_NEAR(vS_M[1], 0.0, 5e-6);
    //Test 3
    EXPECT_NEAR(vP_M[2], 460.894, 460.894*TOL);
    EXPECT_NEAR(vS_M[2], 19.5975, 19.5975*TOL);
    //Test 5
    EXPECT_NEAR(vP_M[3], 1691.64, 1691.64*TOL);
    EXPECT_NEAR(vS_M[3], 8.68975, 8.68975*TOL);
}


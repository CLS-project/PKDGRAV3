#include "gtest/gtest.h"
#define USE_MFM
#include <vector>
#include <cmath>

#include "smooth/smoothfcn.h"
#include "hydro/riemann.h"

const double TOL = 5e-6;

template<typename dtype, typename mtype>
class RiemannTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}

    int width() {
        return dtype::width();
    }


    SMF smf;
    std::vector<double> n_unit {1., 0., 0.};
    dtype vn_unit[3] = {1., 0., 0.};
    dtype zero = 0.;
    double rho_R, p_R;
    double rho_L, p_L;
    std::vector<double> v_R;
    std::vector<double> v_L;
    double P_M, S_M, rho_f, p_f, v_f[3];
    double dConstGamma = 1.4;
    dtype vrho_R, vp_R;
    dtype vrho_L, vp_L;
    dtype vv_R[3];
    dtype vv_L[3];
    dtype vP_M, vS_M, vrho_f, vp_f, vv_f[3];
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
    double cs(double g, double rho, double p) {
        return sqrt(g*p/rho);
    }

    dtype solve() {
        mtype mask = static_cast<dtype>(0.) == 0.;
        RiemannSolverExact<dtype,mtype> riemann(dConstGamma,mask);
        return riemann.solve(vrho_R, vp_R, vv_R,
                             vrho_L, vp_L, vv_L,
                             vP_M, vS_M,
                             &vrho_f, &vp_f, &vv_f[0],
                             vn_unit);
    }

};

typedef RiemannTest<vec<double,double>,mmask<bool>> RiemannTestNoVec;
typedef RiemannTest<dvec,dmask> RiemannTestVec;

#if ( defined(HAVE_MM_POW) || defined(HAVE_MM256_POW) || defined(HAVE_MM512_POW) )
TEST_F(RiemannTestVec, ToroVec1) {

    set_vec_R(0.125, 0.1, {0.,0.,0.});
    set_vec_L(1., 1., {0.,0.,0.});
    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 0.30313, 0.30313*TOL);
        EXPECT_NEAR(vS_M[i], 0.92745, 0.92745*TOL);
    }
}

TEST_F(RiemannTestVec, ToroVec2) {

    set_vec_R(1., 0.4, {2.,0.,0.});
    set_vec_L(1., 0.4, {-2.,0.,0.});
    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 0.00189, 0.000005);
        EXPECT_NEAR(vS_M[i], 0.0, 5e-6);
    }
}

TEST_F(RiemannTestVec, ToroVec3) {

    set_vec_R(1., 0.01, {0.,0.,0.});
    set_vec_L(1., 1000., {0.,0.,0.});
    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 460.894, 460.894*TOL);
        EXPECT_NEAR(vS_M[i], 19.5975, 19.5975*TOL);
    }
}

TEST_F(RiemannTestVec, ToroVec4) {

    set_vec_R(1., 100., {0.,0.,0.});
    set_vec_L(1., 0.01, {0.,0.,0.});
    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i],  46.0950, 46.0950*TOL);
        EXPECT_NEAR(vS_M[i], -6.19633, 6.19633*TOL);
    }
}

TEST_F(RiemannTestVec, ToroVec5) {

    set_vec_R(5.99242, 46.0950, {-6.19633,0.,0.});
    set_vec_L(5.99924, 460.894, {19.5975,0.,0.});
    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        EXPECT_NEAR(vP_M[i], 1691.64, 1691.64*TOL);
        EXPECT_NEAR(vS_M[i], 8.68975, 8.68975*TOL);
    }
}

TEST_F(RiemannTestVec, ToroVecAll) {

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

    auto iters = solve();

    for (auto i=0; i<width(); i++) {
        EXPECT_NE(iters[i], 0);
        printf("%e %e %e\n", iters[i], vP_M[i], vS_M[i]);
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

// Generic test, this was extracted from an Evrard collapse, but can be reused
// if a problematic flux is found...
TEST_F(RiemannTestVec, VecTest) {
    double allrho_R[4] = {4.289399e-02,3.237508e-02, 5.528355e-02,6.015821e-02 };
    double allrho_L[4] = {5.161085e-02,3.163507e-02, 9.228625e-02,9.838759e-02 };
    double allp_R[4] =  {8.342137e-04,8.310585e-04, 1.426528e-02,2.818890e-02 };
    double allp_L[4] =  {0.000000e+00,0.000000e+00, 0.000000e+00,0.000000e+00 };
    double allvx_R[4] = {2.872663e-01,3.538408e-03, 1.822213e-02,3.715746e-02 };
    double allvx_L[4] = {-1.996461e-01,-8.412928e-02,-1.822213e-02,-3.715746e-02 };
    double allvy_R[4] = {-6.350583e-02, 1.505490e-02,-1.020411e-02,-4.192078e-02 };
    double allvy_L[4] = {6.350583e-02,1.175719e-02, 1.020411e-02,4.192078e-02 };
    double allvz_R[4] = {-7.044196e-03,2.253765e-02,-3.883064e-03,-4.353969e-02 };
    double allvz_L[4] = {7.044196e-03,5.634412e-03, 3.883064e-03,4.353969e-02 };
    double n_unitx[4] = {-8.133081e-01,-8.340868e-01, 9.763322e-01,8.057771e-03 };
    double n_unity[4] = {-1.497787e-01,-3.649328e-01, -2.116392e-01, -8.610773e-01 };
    double n_unitz[4] = {-5.622244e-01,-4.136705e-01, -4.454557e-02, -5.084103e-01 };

    vrho_R.load(allrho_R);
    vrho_L.load(allrho_L);
    vp_R.load(allp_R);
    vp_L.load(allp_L);
    vv_R[0].load(allvx_R);
    vv_L[0].load(allvx_L);
    vv_R[1].load(allvy_R);
    vv_L[1].load(allvy_L);
    vv_R[2].load(allvz_R);
    vv_L[2].load(allvz_L);
    vn_unit[0].load(n_unitx);
    vn_unit[1].load(n_unity);
    vn_unit[2].load(n_unitz);

    solve();

    EXPECT_EQ(nan_guard(vS_M, zero), 0);
    EXPECT_EQ(nan_guard(vP_M, zero), 0);
}

TEST_F(RiemannTestVec, NaNguard) {
    // This test does not necessarily mean that the implementation is wrong,
    // rather, that the tests in this suite may present false negatives as
    // some of them rely on detecting nan values.
    // In general, using -fno-finite-math-only will allow for the correct use
    // of NaNguard
    double test = nan("");

    dump(test);
    EXPECT_EQ(nan_guard(test, 0.), 1);

    auto vtest = zero/zero;
    dump(vtest);
    EXPECT_EQ(nan_guard(vtest, zero), 1);
}


TEST_F(RiemannTestVec, VecVacuum) {

    double allrho_R[4] = {0., 1.,   1.,   1.0    };
    double allrho_L[4] = {0.,    1.,   1.,   1.0    };
    double allp_R[4] =   {1.,    1.,  2.0,   2.0   };
    double allp_L[4] =   {1.,    1.0, 1.0, 1.0    };
    double allv_R[4] =   {-19.0, 100., 200., -100.  };
    double allv_L[4] =   {-20.0,-100., 100., -200. };

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

    solve();

    //Test 1 input vacuum
    EXPECT_NEAR(vP_M[0], 0., TOL);
    EXPECT_NEAR(vS_M[0], 0., TOL);
    //Test 2 interal vacuum
    EXPECT_NEAR(vP_M[1], 0.0, TOL);
    EXPECT_NEAR(vS_M[1], 0.0, TOL);
    //Test 3
    double S_L = 100. + (2./(dConstGamma-1.0)) * cs(dConstGamma, 1., 1.);
    EXPECT_NEAR(vP_M[2], 0.0, TOL);
    EXPECT_NEAR(vS_M[2], S_L, S_L*TOL);
    //Test 5
    double S_R = -100. - (2./(dConstGamma-1.0)) * cs(dConstGamma, 1., 2.);
    EXPECT_NEAR(vP_M[3], 0.0, TOL);
    EXPECT_NEAR(vS_M[3], S_R, -S_R*TOL);
}
#endif

TEST_F(RiemannTestNoVec, Toro1) {

    set_vec_R(0.125, 0.1, {0.,0.,0.});
    set_vec_L(1., 1., {0.,0.,0.});
    auto iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(vP_M, 0.30313, 0.30313*TOL);
    EXPECT_NEAR(vS_M, 0.92745, 0.92745*TOL);
}

TEST_F(RiemannTestNoVec, Toro2) {

    set_vec_R(1., 0.4, {2.,0.,0.});
    set_vec_L(1., 0.4, {-2.,0.,0.});
    auto iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(vP_M, 0.00189, 0.000005);
    EXPECT_NEAR(vS_M, 0.0, 5e-6);
}

TEST_F(RiemannTestNoVec, Toro3) {

    set_vec_R(1., 0.01, {0.,0.,0.});
    set_vec_L(1., 1000., {0.,0.,0.});
    auto iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(vP_M, 460.894, 460.894*TOL);
    EXPECT_NEAR(vS_M, 19.5975, 19.5975*TOL);
}

TEST_F(RiemannTestNoVec, Toro4) {

    set_vec_R(1., 100., {0.,0.,0.});
    set_vec_L(1., 0.01, {0.,0.,0.});
    auto iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(vP_M,  46.0950, 46.0950*TOL);
    EXPECT_NEAR(vS_M, -6.19633, 6.19633*TOL);
}

TEST_F(RiemannTestNoVec, Toro5) {

    set_vec_R(5.99242, 46.0950, {-6.19633,0.,0.});
    set_vec_L(5.99924, 460.894, {19.5975,0.,0.});
    auto iters = solve();

    EXPECT_NE(iters, 0);
    EXPECT_NEAR(vP_M, 1691.64, 1691.64*TOL);
    EXPECT_NEAR(vS_M, 8.68975, 8.68975*TOL);
}

TEST_F(RiemannTestNoVec, InputVacuumRho) {
    set_R(0., 46.0950, {-6.19633,0.,0.});
    set_L(0., 460.894, {19.5975,0.,0.});
    int iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 0.0, TOL);
    EXPECT_NEAR(S_M, 0.0, TOL);
}

TEST_F(RiemannTestNoVec, InputVacuumP) {
    set_R(1., 0.0, {-6.19633,0.,0.});
    set_L(1., 0.0, {19.5975,0.,0.});
    int iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 0.0, TOL);
    EXPECT_NEAR(S_M, 0.0, TOL);
}

TEST_F(RiemannTestNoVec, InternalVacuum) {
    set_R(1., 1.0, { 100.,0.,0.});
    set_L(1., 1.0, {-100.,0.,0.});
    int iters = solve();
    EXPECT_NE(iters, 0);
    EXPECT_NEAR(P_M, 0.0, TOL);
    EXPECT_NEAR(S_M, 0.0, TOL);
}


#include "gtest/gtest.h"

#include <cstdlib>
#include <array>
#include "core/simd.h"
#include "hydro/hydro.h"
#include "hydro/limiters.h"

const int Ntests = 200;

class RiemannTest : public ::testing::Test {
protected:
    std::array<double, SIMD_DWIDTH> Lst;
    std::array<double, SIMD_DWIDTH> Rst;
    std::array<double, SIMD_DWIDTH> Lstf;
    std::array<double, SIMD_DWIDTH> Rstf;
    void SetUp() override {
        srand(1234);
    }
    void TearDown() override {}

    double get_rand() {
        return (double)rand()/RAND_MAX * 2.0 - 1.;
    }

    void set_inputs() {
        for (auto k=0; k<SIMD_DWIDTH; k++) {
            Lst[k] = get_rand();
            Rst[k] = get_rand();
            Lstf[k] = get_rand();
            Rstf[k] = get_rand();
            // Some chance to have the same value
            if (get_rand()>0.95)
                Rst[k] = Lst[k];
        }
    }

};

TEST_F(RiemannTest, PairwiseLimiterTest) {

    for (auto i=0; i<Ntests; i++) {
        set_inputs();

        // reference
        std::array<double, SIMD_DWIDTH> Lstf0 = Lstf;
        std::array<double, SIMD_DWIDTH> Rstf0 = Rstf;
        for (auto k=0; k<SIMD_DWIDTH; k++)
            genericPairwiseLimiter(Lst[k], Rst[k], &Lstf0[k], &Rstf0[k]);


        // simd version
        dvec Lstv, Rstv, Lstfv, Rstfv;
        Lstv.load( Lst.data() );
        Rstv.load( Rst.data() );
        Lstfv.load( Lstf.data() );
        Rstfv.load( Rstf.data() );
        genericPairwiseLimiter(Lstv, Rstv, &Lstfv, &Rstfv);

        // compare them
        for (auto k=0; k<SIMD_DWIDTH; k++) {
            //printf("%d %e %e %e\n", k, Lstf[k], Lstf0[k], Lstfv[k]);
            EXPECT_NEAR(Lstfv[k], Lstf0[k], 1e-9);
            EXPECT_NEAR(Rstfv[k], Rstf0[k], 1e-9);
        }

    }
    printf("%d %d %d %d %d %d %d\n", equal, pass1, pass2, pass3, pass4, pass5, pass6);
}

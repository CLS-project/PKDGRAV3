//
// Created by tobias on 11.03.24.
//
#include "SPH/sphcooling/SPHCooling.h"
//namespace SPHCooling {
void pkdSPHCoolingInit(pkdContext *pkd, struct coolingVariables *cool_vars) {
    // conversion factors to convert cooling rate
    pkd->Coolingdata.cooling_rate_to_cgs =
        cool_vars->units.dErgUnit / cool_vars->units.dSecUnit * std::pow(cool_vars->units.dKpcUnit * 3.0857e21, -3);
    pkd->Coolingdata.cgs_to_cooling_rate = 1.0 / pkd->Coolingdata.cooling_rate_to_cgs;
    // conversion factors to convert density
    pkd->Coolingdata.density_to_cgs = cool_vars->units.dGmPerCcUnit;
    pkd->Coolingdata.cgs_to_density = 1.0 / pkd->Coolingdata.density_to_cgs;
}
void pkdSPHCoolingUpdate(PKD pkd, double dScaleFactor) {
}
double pkdSPHCooling(PKD pkd, particleStore::ParticleReference &p, float T, SPHOptions *SPHoptions) {
    // double is needed as floats cannot capture the small scale of the cooling rates
    // (will evaluate to 0 and crash the code, because of division by 0)
    using floattype = double;
    bool DoHeating = false;
    // add heating if cooling model is set to 2
    if (SPHoptions->coolingModel == 2) {
        DoHeating = true;
    }
    /*
    if (!(T > 0)){
        std::cout << "T: " << T << std::endl
                  << "u_old: " << u_old << std::endl
                  << "pDelta: " << pDelta << std::endl
                  << "rho: " << p.density() << std::endl;
    }
    assert(T > 0);
    */
    // convert rho from code units to cgs
    floattype rho = p.density();
    floattype rho_cgs = rho * pkd->Coolingdata.density_to_cgs;
    /// Constants
    // gaunt factor for free-free emission from equation 23 in Katz et al. (1996)
    floattype g_ff = 1.1 + 0.34 * std::exp( - std::pow(5.5 - std::log(T), 2) / 3.0);
    // hydrogen mass fraction
    floattype X = 0.76;
    // helium mass fraction
    floattype Y = 0.24;
    // redshift (needed for compton cooling)
    floattype z = 0.0; // 0 as we are not doing cosmological simulations
    /// Precomputations
    floattype sqrtT = std::sqrt(T);
    floattype sqrtTi = 1.0 / sqrtT;
    // powers of temperature with corresponding exponent needed for cooling rates
    floattype T3 = std::pow(T / 1.0e3, -0.2);
    floattype T5i = 1.0 / (1.0 + std::sqrt(T / 1.0e5));
    floattype T6i = 1.0 / (1.0 + std::pow(T / 1.0e6, 0.7));
    /// Recombination and collisional ionization rates
    // (Table 2 from Katz et al. 1996 (TreeSPH))
    floattype a_Hp, a_Hep, a_d, a_Hepp, gamma_eH0, gamma_eHe0, gamma_eHep;
    a_Hp = 8.4e-11 * sqrtTi * T3 * T6i;
    a_Hep = 1.5e-10 * std::pow(T, -0.6353);
    a_d = 1.9e-3 * std::pow(T, -1.5) * std::exp(-470000.0 / T) * (1.0 + 0.3 * std::exp(-94000.0 / T));
    a_Hepp = 3.36e-10 * sqrtTi * T3 * T6i;
    gamma_eH0 = 5.85e-11 * sqrtT * std::exp(-157809.1/T) * T5i;
    gamma_eHe0 = 2.38e-11 * sqrtT * std::exp(-285335.4/T) * T5i;
    gamma_eHep = 5.68e-12 * sqrtT * std::exp(-631515.0/T) * T5i;
    /// Photoionization rates at z = 0 (from Faucher-Giguère 2020)
    // assuming a homogeneous UV background
    // in units ofs^-1
    floattype gamma_gH0, gamma_gHe0, gamma_gHep;
    if (DoHeating) {
        gamma_gH0 = 3.6200e-14;
        gamma_gHe0 = 1.1910e-14;
        gamma_gHep = 2.8420e-16;
    }
    else {
        gamma_gH0 = 0.0;
        gamma_gHe0 = 0.0;
        gamma_gHep = 0.0;
    }
    /// Abundances of H0, H+, He0, He+, and e-
    // with photoionization
    // iterative approach starting with n_e = n_H using a correction based stopping criteria
    floattype y = Y / (4.0 - 4.0 * Y);
    // hydrogen abundance
    floattype n_H = rho_cgs * X / (1.67262158e-24); // rho * X / m_p
    floattype n_e, n_H0, n_Hep, n_He0, n_Hp, n_Hepp;
    // add a small epsilon to the equations that could get divided by zero to improve stability of the program
    // absolute and relative tolerances and maximum number of iterations for stopping criteria
    floattype abstol = 1e-50;
    floattype reltol = 1e-40;
    unsigned maxiter = 100;
    unsigned iter = 0;
    floattype n_e_prev;
    n_e = n_H;
    floattype eps = 1e-40;
    do {
        n_e_prev = n_e;
        n_H0 = n_H * a_Hp / (a_Hp + gamma_eH0 + gamma_gH0 / n_e);
        n_Hp = n_H - n_H0;
        n_Hep = y * n_H /
                (1.0 + (a_Hep + a_d) / (gamma_eHe0 + gamma_gHe0 / n_e + eps) + (gamma_eHep + gamma_gHep / n_e) / a_Hepp);
        n_He0 = n_Hep * (a_Hep + a_d) / (gamma_eHe0 + gamma_gHe0 / n_e + eps);
        n_Hepp = n_Hep * (gamma_eHep + gamma_gHep / n_e) / a_Hepp;
        n_e = n_Hp + n_Hep + 2.0 * n_Hepp;
        iter++;
        // uncomment this part to see which stopping criterion is applied
        /*
        if (std::abs(n_e_prev - n_e) <= abstol) {
            std::cout << "Stopping criterion: Absolute tolerance reached after " << iter << " iterations." << std::endl;
        }
        if (std::abs(n_e_prev - n_e) <= reltol * n_e) {
            std::cout << "Stopping criterion: Relative tolerance reached after " << iter << " iterations." << std::endl;
        }
        if (iter >= maxiter) {
            std::cout << "Stopping criterion: Maximum iterations reached after " << iter << " iterations." << std::endl;
        }
         */
    } while (std::abs(n_e_prev - n_e) > abstol &&
             std::abs(n_e_prev - n_e) > reltol * n_e &&
             iter < maxiter);
    /// Cooling rates
    /// based on Table 1 from Katz et al. (1996)
    // collisional excitation of H0
    floattype ce_H0 = 7.50e-19 * std::exp(-118348.0 / T) * T5i * n_e * n_H0;
    // collisional excitation of He+
    floattype ce_Hep = 5.54e-17 * std::pow(T, -0.397) * std::exp(-473638.0 / T) * T5i * n_e * n_Hep;
    // collisional ionization of H0
    floattype ci_H0 = 1.27e-21 * sqrtT * std::exp(-157809.1 / T) * T5i * n_e * n_H0;
    // collisional ionization of He0
    floattype ci_He0 = 9.38e-22 * sqrtT * std::exp(-285335.4 / T) * T5i * n_e * n_He0;
    // collisional ionization of He+
    floattype ci_Hep = 4.95e-22 * sqrtT * std::exp(-631515.0 / T) * T5i * n_e * n_Hep;
    // recombination of H+
    floattype r_Hp = 8.70e-27 * sqrtT * T3 * T6i * n_e * n_Hp;
    // recombination of He+
    floattype r_Hep = 1.55e-26 * std::pow(T, 0.3647) * n_e * n_Hep;
    // recombination of He++
    floattype r_Hepp = 3.48e-26 * sqrtT * T3 * T6i * n_e * n_Hepp;
    // dielectric recombination of He+
    floattype dr_Hep = 1.24e-13 * std::pow(T, -1.5) * std::exp(-470000.0 / T)
                       * (1.0 + 0.3 * std::exp(-94000.0 / T)) * n_e * n_Hep;
    // free-free emission (all ions)
    floattype ff = 1.42e-27 * g_ff * sqrtT * (n_Hp + n_Hep + 4.0 * n_Hepp) * n_e;
    // compton cooling of the microwave background from Eq. 24 in Katz et al. (1996)
    floattype compton_cool = 5.41e-36 * n_e * T * std::pow(1.0 + z, 4);
    // heating from photoionization (Eq. 39 and 40 in Katz et al. (1996))
    floattype H;
    if (DoHeating) {
        // Photoheating rates at z = 0 (from Faucher-Giguère 2020)
        // assuming a homogeneous UV background
        // in units of erg s^-1
        floattype e_H0 = 1.3784e-25;
        floattype e_He0 = 8.6360e-26;
        floattype e_Hep = 5.9670e-27;
        H = n_H0 * e_H0 + n_He0 * e_He0 + n_Hep * e_Hep;
    }
    else {
        H = 0.0;
    }
    /// Return value
    // return new internal energy u_new with contribution from heating and cooling
    floattype cooling_rate = H -
                             (ce_H0 + ce_Hep + ci_H0 + ci_He0 + ci_Hep + r_Hp + r_Hep + r_Hepp + dr_Hep + ff + compton_cool);
    // uDot from cooling (H-Lambda/rho)
    float uDot_cool = (cooling_rate * pkd->Coolingdata.cgs_to_cooling_rate) / rho;
    // calculate mean mol weight as sum(n_i * m_i) / (sum(n_i) + n_e)
    // SPHoptions->dMeanMolWeight = (n_H0 + n_Hp + 4 * (n_Hep + n_He0 + n_Hepp)) * 1.67262158e-24 / (n_H + n_Hepp + n_Hp + n_Hep + n_He0 + n_e);
    return uDot_cool;
}
//} // namespace SPHCooling

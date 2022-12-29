#include <iostream>
#include <cassert>
#include <memory>
#include <math.h>

/* This universal IMF module uses dynamic polymorphism to select the
 * requested IMF at run time. Usage is abstracted in a single interface
 * consisting in non-virtual, virtual and pure-virtual member functions
 * in a class hierarchy, which cannot, can and must be overriden,
 * respectively, by deriving classes.
 * Addition of a new universal IMF requires defining a new class that
 * derives from UniversalIMFBaseClass, and adding the corresponding
 * 'else if' clause to the ChooseIMF function.
 */
class UniversalIMFBaseClass {
public:
    UniversalIMFBaseClass() = default;
    virtual ~UniversalIMFBaseClass() = default;

    /* Set of non-virtual member functions implementing general algorithms */
    void SampleLogMass(const double dLowerMass, const double dUpperMass,
                       const int nSamples, double *restrict pdMass) const noexcept {
        assert(dUpperMass > dLowerMass);
        const double dDeltaLog = (log10(dUpperMass) - log10(dLowerMass)) /
                                 (nSamples - 1);
        const double dDelta = pow(10.0, dDeltaLog);
        double dMass = dLowerMass;
        for (int i = 0; i < nSamples; ++i) {
            pdMass[i] = dMass;
            dMass *= dDelta;
        }
    }

    double LogTrapezoidalIntegration(const double *restrict pdXVal,
                                     const double *restrict pdYVal,
                                     const int nPoints) const noexcept {
        const double dDeltaLog = log10(pdXVal[1]) - log10(pdXVal[0]);
        double dSum{};
        for (int i = 1; i < nPoints - 1; ++i) dSum += pdXVal[i] * pdYVal[i];
        dSum += 0.5 * (pdXVal[0] * pdYVal[0] + pdXVal[nPoints-1] * pdYVal[nPoints-1]);
        return dSum * M_LN10 * dDeltaLog;
    }

    void UnweightedSample(const double dLowerMass, const double dUpperMass,
                          const int nSamples, double *restrict pdMass,
                          double *restrict pdSample) const noexcept {
        SampleLogMass(dLowerMass, dUpperMass, nSamples, pdMass);
        for (int i = 0; i < nSamples; ++i) {
            pdSample[i] = Evaluate(pdMass[i]);
        }
    }

    void MassWeightedSample(const double dLowerMass, const double dUpperMass,
                            const int nSamples, double *restrict pdMass,
                            double *restrict pdSample) const noexcept {
        SampleLogMass(dLowerMass, dUpperMass, nSamples, pdMass);
        for (int i = 0; i < nSamples; ++i) {
            pdSample[i] = pdMass[i] * Evaluate(pdMass[i]);
        }
    }

protected:
    /* This member function normalizes a given IMF instance. It should be called
     * only once, preferably in the body of deriving classes' constructors */
    void Normalize(const double dMinMass, const double dMaxMass) {
        const double dNormInv = MassWeightedIntegration(dMinMass, dMaxMass);
        assert(dNormInv > 0.0);
        dNorm = dNorm / dNormInv;
    }

public:
    /* Pure-virtual function that encapsulates the mathematical form of a given IMF */
    virtual double Evaluate(const double dMass) const noexcept = 0;

    /* Set of virtual member functions that can be overriden should an implementation
     * different from the standard be preferable; e.g., analytical instead of numerical
     * integration */
    virtual double UnweightedIntegration(const double dLowerMass,
                                         const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);
        constexpr int nSamples{200};
        double adMass[nSamples], adIMF[nSamples];
        UnweightedSample(dLowerMass, dUpperMass, nSamples, adMass, adIMF);
        return LogTrapezoidalIntegration(adMass, adIMF, nSamples);
    }

    virtual double MassWeightedIntegration(const double dLowerMass,
                                           const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);
        constexpr int nSamples{200};
        double adMass[nSamples], adIMF[nSamples];
        MassWeightedSample(dLowerMass, dUpperMass, nSamples, adMass, adIMF);
        return LogTrapezoidalIntegration(adMass, adIMF, nSamples);
    }

protected:
    double dNorm{1.0};
};


/* From Chabrier 2003 (2003PASP..115..763C), Table 1 */
class Chabrier : public UniversalIMFBaseClass {
public:
    explicit Chabrier(const double dMinMass, const double dMaxMass)
        : Chabrier(0.079, 0.69, -2.3, dMinMass, dMaxMass)
    {}

    explicit Chabrier(const double dMc, const double dSigma, const double dPLSlope,
                      const double dMinMass, const double dMaxMass)
        : dLogMc{log10(dMc)},
          dTwoSigmaSq{2.0 * dSigma * dSigma},
          dPLSlope{dPLSlope} {
        Normalize(dMinMass, dMaxMass);
    }

    double Evaluate(const double dMass) const noexcept {
        if (dMass < 1.0) {
            const double dExpNum = log10(dMass) - dLogMc;
            return dNorm * exp(-dExpNum * dExpNum / dTwoSigmaSq) / dMass;
        }
        else {
            const double dPLFac = exp(-dLogMc * dLogMc / dTwoSigmaSq);
            return dNorm * dPLFac * pow(dMass, dPLSlope);
        }
    }

    double UnweightedIntegration(const double dLowerMass,
                                 const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);

        double dTotalInt{};
        if (dLowerMass < 1.0) {
            const double dLNUpperMass = fmin(dUpperMass, 1.0);
            const double dSqrtTwoSigmaSq = sqrt(dTwoSigmaSq);
            const double dLogNormInt =
                sqrt(0.25 * dTwoSigmaSq * M_PI) *
                (erf((log10(dLNUpperMass) - dLogMc) / dSqrtTwoSigmaSq) -
                 erf((log10(dLowerMass) - dLogMc) / dSqrtTwoSigmaSq));
            dTotalInt += dLogNormInt;
        }
        if (dUpperMass > 1.0) {
            const double dPLLowerMass = fmax(1.0, dLowerMass);
            const double dPLFac = exp(-dLogMc * dLogMc / dTwoSigmaSq);
            const double dPowerLawInt = (dPLFac / (dPLSlope + 1.0)) *
                                        (pow(dUpperMass, dPLSlope + 1.0) -
                                         pow(dPLLowerMass, dPLSlope + 1.0));
            dTotalInt += dPowerLawInt;
        }

        return dNorm * dTotalInt;
    }

    double MassWeightedIntegration(const double dLowerMass,
                                   const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);

        const double dFac1 = exp(-dLogMc * dLogMc / dTwoSigmaSq);
        const double dFac2 = dLogMc + 0.5 * dTwoSigmaSq * M_LN10;

        double dTotalInt{};
        if (dLowerMass < 1.0) {
            const double dLNUpperMass = fmin(dUpperMass, 1.0);
            const double dSqrtTwoSigmaSq = sqrt(dTwoSigmaSq);
            const double dLogNormInt =
                sqrt(0.25 * dTwoSigmaSq * M_PI) * M_LN10 * dFac1 *
                exp(dFac2 * dFac2 / dTwoSigmaSq) *
                (erf((log10(dLNUpperMass) - dFac2) / dSqrtTwoSigmaSq) -
                 erf((log10(dLowerMass) - dFac2) / dSqrtTwoSigmaSq));
            dTotalInt += dLogNormInt;
        }
        if (dUpperMass > 1.0) {
            const double dPLLowerMass = fmax(1.0, dLowerMass);
            const double dPowerLawInt = (dFac1 / (dPLSlope + 2.0)) *
                                        (pow(dUpperMass, dPLSlope + 2.0) -
                                         pow(dPLLowerMass, dPLSlope + 2.0));
            dTotalInt += dPowerLawInt;
        }

        return dNorm * dTotalInt;
    }

private:
    double dLogMc;
    double dTwoSigmaSq;
    double dPLSlope;
};


/* From Kroupa 2001 (2001MNRAS.322..231K), Eq. 2 */
class Kroupa : public UniversalIMFBaseClass {
public:
    explicit Kroupa(const double dMinMass, const double dMaxMass)
        : Kroupa(-0.3, -1.3, -2.3, dMinMass, dMaxMass)
    {}

    explicit Kroupa(const double dAlpha0, const double dAlpha1, const double dAlpha2,
                    const double dMinMass, const double dMaxMass)
        : dPLSlope0{dAlpha0},
          dPLSlope1{dAlpha1},
          dPLSlope2{dAlpha2} {
        Normalize(dMinMass, dMaxMass);
    }

    double Evaluate(const double dMass) const noexcept {
        if (dMass < 0.08) {
            const double dPLFac = dNorm * pow(0.08, dPLSlope1 - dPLSlope0) *
                                  pow(0.5, dPLSlope2 - dPLSlope1);
            return dPLFac * pow(dMass, dPLSlope0);
        }
        else if (dMass < 0.5) {
            const double dPLFac = dNorm * pow(0.5, dPLSlope2 - dPLSlope1);
            return dPLFac * pow(dMass, dPLSlope1);
        }
        else {
            return dNorm * pow(dMass, dPLSlope2);
        }
    }

private:
    double dPLSlope0;
    double dPLSlope1;
    double dPLSlope2;
};


/* From Salpeter 1955 (1955ApJ...121..161S), Eq. 5 */
class Salpeter : public UniversalIMFBaseClass {
public:
    explicit Salpeter(const double dMinMass, const double dMaxMass)
        : Salpeter(-2.35, dMinMass, dMaxMass)
    {}

    explicit Salpeter(const double dAlpha, const double dMinMass,
                      const double dMaxMass)
        : dPLSlope{dAlpha} {
        Normalize(dMinMass, dMaxMass);
    }

    double Evaluate(const double dMass) const noexcept {
        return dNorm * pow(dMass, dPLSlope);
    }

private:
    double dPLSlope;
};


/* Function that selects the requested IMF and returns a unique_ptr
 * which can then be used to call the available algorithms.
 * For example:
 * {
 *      ...
 *      auto IMF = ChooseIMF( ... );
 *      IMF->AnAlgorithm( ... );
 *      IMF->AnotherAlgorithm( ... );
 *      ...
 * }
 */
inline std::unique_ptr<UniversalIMFBaseClass>
ChooseIMF(const char *pszIMFType, const double dIMFMinMass, const double dIMFMaxMass) {
    if (strcmp(pszIMFType, "chabrier") == 0) {
        return std::make_unique<Chabrier>(dIMFMinMass, dIMFMaxMass);
    }
    else if (strcmp(pszIMFType, "kroupa") == 0) {
        return std::make_unique<Kroupa>(dIMFMinMass, dIMFMaxMass);
    }
    else if (strcmp(pszIMFType, "salpeter") == 0) {
        return std::make_unique<Salpeter>(dIMFMinMass, dIMFMaxMass);
    }
    else {
        std::cerr << "ERROR: Undefined IMF type has been given to ChooseIMF: " <<
                  pszIMFType << std::endl;
        exit(1);
    }
}


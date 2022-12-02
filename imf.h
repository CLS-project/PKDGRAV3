#include <iostream>
#include <cassert>
#include <memory>
#include <math.h>

class UniversalIMFBaseClass {
public:
    UniversalIMFBaseClass() = default;
    virtual ~UniversalIMFBaseClass() = default;

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
    void Normalize(const double dMinMass, const double dMaxMass) {
        const double dNormInv = MassWeightedIntegration(dMinMass, dMaxMass);
        assert(dNormInv > 0.0);
        dNorm = 1.0 / dNormInv;
    }

public:
    virtual double Evaluate(const double dMass) const noexcept = 0;

    virtual double UnweightedIntegration(const double dLowerMass,
                                         const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);
        constexpr int nSamples = 200;
        double adMass[nSamples], adIMF[nSamples];
        UnweightedSample(dLowerMass, dUpperMass, nSamples, adMass, adIMF);
        return LogTrapezoidalIntegration(adMass, adIMF, nSamples);
    }

    virtual double MassWeightedIntegration(const double dLowerMass,
                                           const double dUpperMass) const noexcept {
        assert(dUpperMass > dLowerMass);
        constexpr int nSamples = 200;
        double adMass[nSamples], adIMF[nSamples];
        MassWeightedSample(dLowerMass, dUpperMass, nSamples, adMass, adIMF);
        return LogTrapezoidalIntegration(adMass, adIMF, nSamples);
    }

protected:
    double dNorm{1.0};
};


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


class Kroupa : public UniversalIMFBaseClass {
public:
    explicit Kroupa(const double dMinMass, const double dMaxMass) {
        Normalize(dMinMass, dMaxMass);
    }

    double Evaluate(const double dMass) const noexcept {
        return 1.0;
    }

    double UnweightedIntegration(const double dLowerMass, const double dUpperMass) {
        return 1.0;
    }

    double MassWeightedIntegration(const double dLowerMass, const double dUpperMass) {
        return 1.0;
    }
};


class Salpeter : public UniversalIMFBaseClass {
public:
    explicit Salpeter(const double dMinMass, const double dMaxMass) {
        Normalize(dMinMass, dMaxMass);
    }

    double Evaluate(const double dMass) const noexcept {
        return 1.0;
    }

    double UnweightedIntegration(const double dLowerMass, const double dUpperMass) {
        return 1.0;
    }

    double MassWeightedIntegration(const double dLowerMass, const double dUpperMass) {
        return 1.0;
    }
};


inline std::unique_ptr<UniversalIMFBaseClass>
ChooseIMF(const char *pszIMFType, const double dIMFMinMass, const double dIMFMaxMass) {
    using IMFBase = UniversalIMFBaseClass;

    if (strcmp(pszIMFType, "chabrier") == 0) {
        return std::unique_ptr<IMFBase>(new Chabrier{dIMFMinMass, dIMFMaxMass});
    }
    else if (strcmp(pszIMFType, "kroupa") == 0) {
        return std::unique_ptr<IMFBase>(new Kroupa{dIMFMinMass, dIMFMaxMass});
    }
    else if (strcmp(pszIMFType, "salpeter") == 0) {
        return std::unique_ptr<IMFBase>(new Salpeter{dIMFMinMass, dIMFMaxMass});
    }
    else {
        std::cerr << "ERROR: Undefined IMF type has been given to ChooseIMF: " <<
                  pszIMFType << std::endl;
        exit(1);
    }
}


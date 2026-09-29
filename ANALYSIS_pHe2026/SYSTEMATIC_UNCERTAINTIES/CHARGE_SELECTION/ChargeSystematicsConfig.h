#ifndef DAMPE_CHARGE_SYSTEMATICS_CONFIG_H
#define DAMPE_CHARGE_SYSTEMATICS_CONFIG_H

#include <string>
#include <stdexcept>
#include <cmath>

// Positive widthFraction = wider charge window; negative = narrower.
// Strength = 1 preserves the existing affine MC correction exactly.
// Strength = 0 disables it. Values 0.9/1.1 are ONLY stress tests,
// NOT 1-sigma uncertainties without external calibration information.
struct ChargeSysConfig {
    std::string tag;
    double psdWidthFraction = 0.0;
    double stkWidthFraction = 0.0;
    double psdCorrectionStrength = 1.0;
    double stkCorrectionStrength = 1.0;
};

inline ChargeSysConfig MakeChargeSysConfig(const std::string& tag) {
    ChargeSysConfig c;
    c.tag = tag;

    if      (tag == "nominal") {}
    else if (tag == "psd_tight05")  c.psdWidthFraction = -0.05;
    else if (tag == "psd_loose05")  c.psdWidthFraction = +0.05;
    else if (tag == "stk_tight05")  c.stkWidthFraction = -0.05;
    else if (tag == "stk_loose05")  c.stkWidthFraction = +0.05;
    else if (tag == "both_tight05") { c.psdWidthFraction = -0.05; c.stkWidthFraction = -0.05; }
    else if (tag == "both_loose05") { c.psdWidthFraction = +0.05; c.stkWidthFraction = +0.05; }
    else if (tag == "psd_tight10")  c.psdWidthFraction = -0.10;
    else if (tag == "psd_loose10")  c.psdWidthFraction = +0.10;
    else if (tag == "stk_tight10")  c.stkWidthFraction = -0.10;
    else if (tag == "stk_loose10")  c.stkWidthFraction = +0.10;
    else if (tag == "both_tight10") { c.psdWidthFraction = -0.10; c.stkWidthFraction = -0.10; }
    else if (tag == "both_loose10") { c.psdWidthFraction = +0.10; c.stkWidthFraction = +0.10; }
    // MC-only correction-strength variations (diagnostics, not 1 sigma).
    else if (tag == "psd_corr90")  c.psdCorrectionStrength = 0.9;
    else if (tag == "psd_corr110") c.psdCorrectionStrength = 1.1;
    else if (tag == "psd_corroff") c.psdCorrectionStrength = 0.0;
    else if (tag == "stk_corr90")  c.stkCorrectionStrength = 0.9;
    else if (tag == "stk_corr110") c.stkCorrectionStrength = 1.1;
    else if (tag == "stk_corroff") c.stkCorrectionStrength = 0.0;
    else if (tag == "both_corroff") { c.psdCorrectionStrength = 0.0; c.stkCorrectionStrength = 0.0; }
    else throw std::invalid_argument("Unknown charge-systematics scenario: " + tag);

    return c;
}

inline void VaryChargeWindow(double &low, double &high, double widthFraction) {
    if (!std::isfinite(low) || !std::isfinite(high) || !(high > low) || widthFraction <= -1.0)
        throw std::invalid_argument("Invalid charge selection window/width variation");
    if (widthFraction == 0.0) return; // preserve nominal limits exactly
    const double center = 0.5 * (low + high);
    const double halfWidth = 0.5 * (high - low) * (1.0 + widthFraction);
    low  = center - halfWidth;
    high = center + halfWidth;
}

inline double VaryCorrectionStrength(double rawCharge, double nominalCorrectedCharge,
                                     double strength) {
    if (strength == 1.0) return nominalCorrectedCharge; // preserve original MC nominal exactly
    if (strength == 0.0) return rawCharge;
    return rawCharge + strength * (nominalCorrectedCharge - rawCharge);
}

#endif


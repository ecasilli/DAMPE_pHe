#ifndef DAMPE_PHE_RATIO_H
#define DAMPE_PHE_RATIO_H
#include <math.h>
#include <float.h>
#include <stdexcept>

// ROOT 5 / CINT variant. Identical numerical model to the first version.
// Energy argument: kinetic energy per nucleus in GeV.
// Log-log interpolation of the two supplied SEP2026 fluxes.
namespace DAMPEpHe {
inline int size() { return 42; }
inline double minEnergyGeV() { return 49.8588; }
inline double maxEnergyGeV() { return 789808.0; }
inline bool outsideSupport(double T) {
    return T < minEnergyGeV() || T > maxEnergyGeV();
}
inline double ratio(double T) {
    // No isfinite (C++11); this also rejects NaN and infinity.
    if (!(T > 0.) || T > DBL_MAX)
        throw std::invalid_argument("MC energy must be positive and finite");

    // Simple arrays, explicitly assigned once, for old CINT compatibility.
    static double logT[42];
    static double logR[42];
    static int initialized = 0;
    if (!initialized) {
        logT[0] = 3.9091950104170992; logR[0] = 0.49247090848393316;
        logT[1] = 4.3692047849538236; logR[1] = 0.41096336894354346;
        logT[2] = 4.3697085779688445; logR[2] = 0.41091669222222249;
        logT[3] = 4.8297207469997563; logR[3] = 0.34851818756181263;
        logT[4] = 4.8302207178585164; logR[4] = 0.34847668036038781;
        logT[5] = 5.2902398309186189; logR[5] = 0.309186023987988;
        logT[6] = 5.290733682974615; logR[6] = 0.30915782290843374;
        logT[7] = 5.7507551151137495; logR[7] = 0.29524044575819186;
        logT[8] = 5.7512466807665801; logR[8] = 0.29521231002073911;
        logT[9] = 6.2112725416326899; logR[9] = 0.28457728002267224;
        logT[10] = 6.2117594448586084; logR[10] = 0.28454982042905907;
        logT[11] = 6.6717898779478686; logR[11] = 0.27702109764002891;
        logT[12] = 6.6722721573539499; logR[12] = 0.2769937577488264;
        logT[13] = 7.1323058399938022; logR[13] = 0.2658107158401517;
        logT[14] = 7.132784970415214; logR[14] = 0.26579444852231049;
        logT[15] = 7.5928249239126648; logR[15] = 0.26131577498722258;
        logT[16] = 7.5932978174087706; logR[16] = 0.26128501001731408;
        logT[17] = 8.0533402081077945; logR[17] = 0.24392161923903544;
        logT[18] = 8.0538106344150648; logR[18] = 0.24388192241608131;
        logT[19] = 8.5138576346267349; logR[19] = 0.21046289419007813;
        logT[20] = 8.5143234376089261; logR[20] = 0.21041132839186538;
        logT[21] = 8.9743749709419145; logR[21] = 0.1509661514251146;
        logT[22] = 8.9748362551504464; logR[22] = 0.15088537769010024;
        logT[23] = 9.4348909329878481; logR[23] = 0.053777596285849327;
        logT[24] = 9.4353490628357353; logR[24] = 0.053673730686554677;
        logT[25] = 9.8954100169067107; logR[25] = -0.091368098251798813;
        logT[26] = 9.895861870791272; logR[26] = -0.091467304106462421;
        logT[27] = 10.355925301101841; logR[27] = -0.2724914417872597;
        logT[28] = 10.356374679806208; logR[28] = -0.27259380985070791;
        logT[29] = 10.816442727620782; logR[29] = -0.4457973774464854;
        logT[30] = 10.816887490087714; logR[30] = -0.44593380327294696;
        logT[31] = 11.27696006393596; logR[31] = -0.54698699409686213;
        logT[32] = 11.277400299352061; logR[32] = -0.54702807531286624;
        logT[33] = 11.737476025981893; logR[33] = -0.50029770328655232;
        logT[34] = 11.737913109299237; logR[34] = -0.50026622590110392;
        logT[35] = 12.197995109900756; logR[35] = -0.39274583740798619;
        logT[36] = 12.198425919196559; logR[36] = -0.39265549805485733;
        logT[37] = 12.658510394095886; logR[37] = -0.20594912662695464;
        logT[38] = 12.658938728887062; logR[38] = -0.20585448707595333;
        logT[39] = 13.119027820614827; logR[39] = -0.059839700655317074;
        logT[40] = 13.119451538914458; logR[40] = -0.059731231158785647;
        logT[41] = 13.579545156930006; logR[41] = -0.08966919069686341;
        initialized = 1;
    }
    const double x = log(T);
    if (x <= logT[0]) return exp(logR[0]);
    if (x >= logT[41]) return exp(logR[41]);

    // Binary search without STL comparator/lambda.
    int lo = 0;
    int hi = 41;
    while (hi-lo > 1) {
        const int mid = (lo+hi)/2;
        if (x < logT[mid]) hi = mid;
        else lo = mid;
    }
    const double u = (x-logT[lo])/(logT[hi]-logT[lo]);
    return exp(logR[lo] + u*(logR[hi]-logR[lo]));
}
inline double compositionWeight(double T, bool isHelium) {
    const double r = ratio(T);
    return isHelium ? 2./(1.+r) : 2.*r/(1.+r);
}
}
#endif


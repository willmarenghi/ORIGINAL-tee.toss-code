#pragma once
#include <algorithm>
#include <array>
#include <cmath>

namespace kit {

// ------------------------------------------------------------
// Metal-only coefficients (MVP)
// ------------------------------------------------------------
struct Coeff {
    double e;            // collision efficiency (ball-bat COR-like term)
    double k_potential;  // effective bat-efficiency factor for "potential" line
    double c;            // intercept (m/s); keep 0 for MVP
};

// Default metal bat coefficients 
inline constexpr Coeff kMetal{ /*e*/0.28, /*k_potential*/0.92, /*c*/0.0 };
// Default wood bat coefficients
inline constexpr Coeff kWood{ /*e*/0.23, /*k_potential*/0.92, /*c*/0.0 };

// ACTUAL MODEL NUMBERS

// Tee Mode:
//     Wood e_value = 0.
//     Wood K_potential = 1.21
//     Metal e_value = 0.26
//     Metal K_potential = 1.10
//         plateVelo = 0mph or 0m/s

// Toss Mode:
//     Wood e_value = 0.22
//     Wood k_potential = 1.18
//     Metal e_value = 0.26
//     Metal k_potential = 1.10
//         plateVelo = 5.36 mps (12 mph)

// ------------------------------------------------------------
// Spin × Distance drag (plate/release speed ratio, unitless)
//
// Spin bins: [1000,1200,1400,1600,1800,2000,2200,2400,2600)
// Distances snap to nearest of {12.19, 13.41, 14.02, 14.63, 15.24, 16.46} m
// ------------------------------------------------------------
struct MultRow { double dist_m; std::array<double,8> m; };

inline constexpr std::array<int,9> kSpinBins = {
    1000,1200,1400,1600,1800,2000,2200,2400,2600
};

// Converted from feet to meters (1 ft = 0.3048 m)
inline constexpr std::array<MultRow,6> kMultiplier = {{
    {12.192, {0.969,0.970,0.972,0.974,0.975,0.976,0.978,0.980}}, // 40 ft
    {13.411, {0.964,0.966,0.968,0.970,0.972,0.973,0.975,0.977}}, // 44 ft
    {14.021, {0.956,0.959,0.961,0.962,0.964,0.965,0.967,0.969}}, // 46 ft
    {14.630, {0.942,0.945,0.948,0.951,0.953,0.955,0.958,0.961}}, // 48 ft
    {15.240, {0.927,0.931,0.935,0.939,0.942,0.945,0.948,0.952}}, // 50 ft
    {16.459, {0.901,0.906,0.910,0.914,0.918,0.922,0.926,0.931}}, // 54 ft
}};

inline int spinIndex(int rpm) {
    int clamped = std::max(1000, std::min(rpm, 2599));
    for (int i = 0; i < 8; ++i) {
        if (clamped >= kSpinBins[i] && clamped < kSpinBins[i+1]) return i;
    }
    return 7; // highest bin if out-of-range
}

inline double nearestDistanceKey(double dist_m) {
    double best = kMultiplier[0].dist_m;
    double bestDiff = std::abs(dist_m - best);
    for (const auto& row : kMultiplier) {
        double d = std::abs(dist_m - row.dist_m);
        if (d < bestDiff) {
            best = row.dist_m;
            bestDiff = d;
        }
    }
    return best;
}

inline double plateMultiplier(double dist_m, int spin_rpm) {
    const double dkey = nearestDistanceKey(dist_m);
    const int idx  = spinIndex(spin_rpm);
    for (const auto& row : kMultiplier) {
        if (std::abs(row.dist_m - dkey) < 0.001) return row.m[idx];
    }
    // Fallback: roughly 16.46m (54'), mid spin
    return 0.914;
}

// ------------------------------------------------------------
// Public API
// ------------------------------------------------------------
struct Inputs {
    double batSpeed_mps;      // BS at impact (m/s)
    double pitchRelease_mps;  // V_release (m/s)
    double distance_m;        // distance in meters (snaps to nearest key)
    int    spin_rpm;          // pitch spin (rpm)
    double EV_measured_mps;   // measured EV (m/s)
};

struct Outputs {
    // Telemetry for debugging/tuning
    double distance_key_m;
    int    spin_bin_index;
    double multiplier_used;     // unitless drag multiplier
    double plateVelocity_mps;   // V_plate = V_release * multiplier

    // Model outputs
    double potentialEV_mps;     // EV_pot = k(1+e)*BS + e*V_plate + c
    double potentialSmash;      // EV_pot / BS

    // Derived from measured vs model
    double smash_measured;      // EV_measured / BS
    double squaredUp_pct_raw;   // 100 * EV_measured / EV_pot
    double squaredUp_pct_ui;    // UI cap only: min(raw, 100)
};

// Main computation for metal/wood bats
inline Outputs computePotential(const Inputs& in, const Coeff batCoeff) {
    Outputs out{};

    // 1) Drag-adjusted plate velocity
    out.distance_key_m    = nearestDistanceKey(in.distance_m);
    out.spin_bin_index    = spinIndex(in.spin_rpm);
    out.multiplier_used   = plateMultiplier(in.distance_m, in.spin_rpm);
    out.plateVelocity_mps = in.pitchRelease_mps * out.multiplier_used;

    // 2) Potential EV from collision model:
    //    EV_pot = k(1+e)*BS + e*V_plate + c
    const double BS  = std::max(1e-6, in.batSpeed_mps);
    out.potentialEV_mps = batCoeff.k_potential * (1.0 + batCoeff.e) * BS
                          + batCoeff.e * out.plateVelocity_mps
                          + batCoeff.c;

    // 3) Derived metrics
    out.potentialSmash    = out.potentialEV_mps / BS;
    out.smash_measured    = in.EV_measured_mps / BS;

    const double denom    = std::max(1e-6, out.potentialEV_mps);
    out.squaredUp_pct_raw = (in.EV_measured_mps / denom) * 100.0;
    out.squaredUp_pct_ui  = std::min(out.squaredUp_pct_raw, 100.0); // UI cap only

    return out;
}

#if 0 // Enable for cage debugging
#include <cstdio>
inline void debugPrint(const Inputs& in, const Outputs& out) {
  std::printf("[KIT] d=%.2fm→%.2fm spin=%drpm bin=%d mult=%.3f | "
              "Vrel=%.1f Vplate=%.1f | EV_meas=%.1f EV_pot=%.1f | "
              "Sm_meas=%.3f Sm_pot=%.3f SqUp_raw=%.1f%%\n",
              in.distance_m, out.distance_key_m, in.spin_rpm, out.spin_bin_index,
              out.multiplier_used,
              in.pitchRelease_mps, out.plateVelocity_mps,
              in.EV_measured_mps, out.potentialEV_mps,
              out.smash_measured, out.potentialSmash, out.squaredUp_pct_raw);
}
#endif

} // namespace kit

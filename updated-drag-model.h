#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

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
//     Wood e_value = 0.22
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
// Physics-based drag model
//
// Replaces the old spin-bin × distance lookup table with a
// continuous formula using exact spin RPM and exact distance.
//
// Formula:
//   omega     = spin_rpm * 2*pi/60                    (rad/s)
//   S         = r * omega / v_release                 (spin parameter, dimensionless)
//   CD0       = cd0ForVelocityAndDistance(v, dist_m)  (2D baseline drag)
//   CD_eff    = CD0 - kCDspin * S                    (effective drag coeff)
//   k_eff     = kAeroFactor * CD_eff                 (drag decay constant, 1/m)
//   multiplier = exp(-k_eff * dist_m)
//
// CD0 is a 2D lookup: velocity × distance.
//   Slow anchor (40mph): calibrated from 45 Stalker+KIT pitches (June 2026).
//                        RMS error = 1.37 mph across 25-44ft.
//   Fast anchor (90mph): back-calculated so 90mph/2200rpm/60.5ft = 10% loss
//                        (MLB standard). CD0_fast = 0.5199 (flat across distances
//                        — no high-speed distance data yet).
//   Between anchors:     linearly interpolated by release velocity.
//   Below 40mph:         clamped to slow anchor values.
//   Above 90mph:         clamped to fast anchor values.
//
// Why two anchors? The drag coefficient drops significantly as velocity rises
// (Reynolds number effect / drag crisis). A baseball at 40mph experiences
// roughly 2.6x more drag than at 90mph.
//
// Supported distance range: 25-54 ft (7.62-16.46 m)
// Velocity range: calibrated 30-47mph (slow), anchored at 90mph (fast)
// Spin range: any RPM
// ------------------------------------------------------------
namespace drag {
    constexpr double kBallRadius_m = 0.03683;
    constexpr double kAirDensity   = 1.225;
    constexpr double kBallArea_m2  = 3.14159265358979 * kBallRadius_m * kBallRadius_m;
    constexpr double kBallMass_kg  = 0.1417;
    constexpr double kAeroFactor   = kAirDensity * kBallArea_m2 / (2.0 * kBallMass_kg);

    // Spin correction — higher spin → smaller CD (backspin lift effect).
    // Calibrated via least-squares regression on 45 pitches (Stalker + KIT radar, 2026-06).
    // CD_eff = CD0(v, dist) - kCDspin * S,  where S = r*omega/v
    constexpr double kCDspin = 0.9946;

    // 2D CD0 table: two velocity anchors × 8 distances.
    // CD0 interpolates linearly between kVelSlow and kVelFast by release velocity.
    constexpr double kVelSlow_mps = 17.882;  // 40 mph — slow-toss anchor
    constexpr double kVelFast_mps = 40.234;  // 90 mph — fast-pitch anchor

    struct CD0Point { double dist_m; double cd0_slow; double cd0_fast; };
    inline constexpr CD0Point kCD0Table[] = {
        //  dist_m    slow(40mph)  fast(90mph)
        {  7.620,    0.7638,      0.5199 },  // 25ft — slow: calibrated; fast: anchored to MLB 10%
        {  9.144,    0.7934,      0.5199 },  // 30ft — slow: calibrated
        { 10.668,    0.7654,      0.5199 },  // 35ft — slow: calibrated
        { 12.192,    0.8580,      0.5199 },  // 40ft — slow: calibrated
        { 13.411,    0.8673,      0.5199 },  // 44ft — slow: calibrated
        { 14.630,    0.8673,      0.5199 },  // 48ft — slow: PLACEHOLDER (needs data)
        { 15.240,    0.8673,      0.5199 },  // 50ft — slow: PLACEHOLDER (needs data)
        { 16.459,    0.8673,      0.5199 },  // 54ft — slow: PLACEHOLDER (needs data)
    };
    inline constexpr int kCD0TableSize = 8;
}

// Returns the 2D CD0 baseline drag coefficient.
// Interpolates linearly along the distance axis, then blends between the
// slow (40mph) and fast (90mph) velocity anchors.
inline double cd0ForVelocityAndDistance(double v_mps, double dist_m) {
    const auto* t  = drag::kCD0Table;
    const int   n  = drag::kCD0TableSize;

    // 1) Interpolate distance axis to get cd0_slow and cd0_fast at this dist
    double cd0_slow, cd0_fast;
    if (dist_m <= t[0].dist_m) {
        cd0_slow = t[0].cd0_slow;
        cd0_fast = t[0].cd0_fast;
    } else if (dist_m >= t[n-1].dist_m) {
        cd0_slow = t[n-1].cd0_slow;
        cd0_fast = t[n-1].cd0_fast;
    } else {
        cd0_slow = t[n-1].cd0_slow;
        cd0_fast = t[n-1].cd0_fast;
        for (int i = 0; i + 1 < n; ++i) {
            if (dist_m >= t[i].dist_m && dist_m <= t[i+1].dist_m) {
                double frac = (dist_m - t[i].dist_m) / (t[i+1].dist_m - t[i].dist_m);
                cd0_slow = t[i].cd0_slow + frac * (t[i+1].cd0_slow - t[i].cd0_slow);
                cd0_fast = t[i].cd0_fast + frac * (t[i+1].cd0_fast - t[i].cd0_fast);
                break;
            }
        }
    }

    // 2) Blend between slow and fast anchors based on release velocity
    const double v_clamped = std::max(drag::kVelSlow_mps,
                                      std::min(v_mps, drag::kVelFast_mps));
    const double vfrac = (v_clamped - drag::kVelSlow_mps)
                       / (drag::kVelFast_mps - drag::kVelSlow_mps);
    return cd0_slow + vfrac * (cd0_fast - cd0_slow);
}

// Debug telemetry from the drag calculation
struct DragDebug {
    double omega_rads;   // angular velocity (rad/s)
    double spin_param;   // S = r*omega/v (dimensionless)
    double CD0_used;     // distance-dependent baseline before spin correction
    double CD_eff;       // effective drag coefficient after spin correction
    double k_eff;        // drag decay constant (1/m)
};

// plateMultiplierAdjusted — continuous physics-based drag multiplier.
// CD0 is looked up from the distance-dependent calibration table.
// Pass a non-null dbg pointer to capture intermediate values for logging.
inline double plateMultiplierAdjusted(double dist_m, int spin_rpm,
                                       double pitchRelease_mps,
                                       DragDebug* dbg = nullptr) {
    const double v     = std::max(pitchRelease_mps, 0.1);
    const double omega = spin_rpm * (2.0 * M_PI / 60.0);
    const double S     = drag::kBallRadius_m * omega / v;
    const double CD0   = cd0ForVelocityAndDistance(v, dist_m);
    const double CD    = std::max(0.05, CD0 - drag::kCDspin * S);
    const double k     = drag::kAeroFactor * CD;
    if (dbg) {
        dbg->omega_rads = omega;
        dbg->spin_param = S;
        dbg->CD0_used   = CD0;
        dbg->CD_eff     = CD;
        dbg->k_eff      = k;
    }
    return std::exp(-k * dist_m);
}

// ------------------------------------------------------------
// Public API
// ------------------------------------------------------------
struct Inputs {
    double batSpeed_mps;      // BS at impact (m/s)
    double pitchRelease_mps;  // V_release (m/s)
    double distance_m;        // pitch distance in meters (25-54 ft / 7.62-16.46 m)
    int    spin_rpm;          // pitch spin, exact RPM from radar
    double EV_measured_mps;   // measured EV (m/s)
};

struct Outputs {
    // Drag model telemetry
    DragDebug drag;              // omega, spin_param, CD_eff, k_eff
    double multiplier_used;      // unitless drag multiplier
    double plateVelocity_mps;    // V_plate = V_release * multiplier

    // Model outputs
    double potentialEV_mps;      // EV_pot = k(1+e)*BS + e*V_plate + c
    double potentialSmash;       // EV_pot / BS

    // Derived from measured vs model
    double smash_measured;       // EV_measured / BS
    double squaredUp_pct_raw;    // 100 * EV_measured / EV_pot
    double squaredUp_pct_ui;     // UI cap only: min(raw, 100)
};

// Main computation for metal/wood bats
inline Outputs computePotential(const Inputs& in, const Coeff batCoeff) {
    Outputs out{};

    // 1) Drag-adjusted plate velocity
    out.multiplier_used   = plateMultiplierAdjusted(in.distance_m, in.spin_rpm,
                                                     in.pitchRelease_mps, &out.drag);
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
    out.squaredUp_pct_ui  = std::min(out.squaredUp_pct_raw, 100.0);

    return out;
}

#if 0 // Enable for cage debugging
#include <cstdio>
inline void debugPrint(const Inputs& in, const Outputs& out) {
  std::printf("[KIT] d=%.2fm spin=%drpm(%.1frad/s) S=%.4f CD0=%.4f CD_eff=%.4f k=%.5f/m | "
              "mult=%.4f Vrel=%.1f Vplate=%.1f | "
              "EV_meas=%.1f EV_pot=%.1f | "
              "Sm_meas=%.3f Sm_pot=%.3f SqUp_raw=%.1f%%\n",
              in.distance_m, in.spin_rpm,
              out.drag.omega_rads, out.drag.spin_param,
              out.drag.CD0_used, out.drag.CD_eff, out.drag.k_eff,
              out.multiplier_used,
              in.pitchRelease_mps, out.plateVelocity_mps,
              in.EV_measured_mps, out.potentialEV_mps,
              out.smash_measured, out.potentialSmash, out.squaredUp_pct_raw);
}
#endif

// ------------------------------------------------------------
// Unit tests  (compile with: #define KIT_RUN_TESTS 1)
// ------------------------------------------------------------
#ifdef KIT_RUN_TESTS

// Distances to test (ft → m): 25,30,35,40,44,48,50,54
inline constexpr double kTestDists_m[] = {
    7.620, 9.144, 10.668, 12.192, 13.411, 14.630, 15.240, 16.459
};
inline constexpr int kTestSpins[] = { 1000, 1700, 2600 };          // low, mid, high
inline constexpr double kTestVelos[] = { 11.18, 17.88, 26.82 };   // 25, 40, 60 mph

static void kitRunTests() {
    int pass = 0, fail = 0;
    const char* spinLabel[] = { "low ", "mid ", "high" };
    const char* veloLabel[] = { "25mph", "40mph", "60mph" };

    std::printf("=== KIT Drag Model Unit Tests ===\n");
    std::printf("%-6s %-5s %-5s  %-8s  %s\n",
                "dist_m", "spin", "velo", "mult", "checks");

    for (double d : kTestDists_m) {
        for (int si = 0; si < 3; ++si) {
            for (int vi = 0; vi < 3; ++vi) {
                DragDebug dbg{};
                double mult = plateMultiplierAdjusted(d, kTestSpins[si],
                                                      kTestVelos[vi], &dbg);
                bool ok = true;
                // multiplier must be (0, 1]
                if (mult <= 0.0 || mult > 1.0) ok = false;
                // must be closer to 1 at shorter distance than at 54ft
                double mult54 = plateMultiplierAdjusted(16.459, kTestSpins[si],
                                                         kTestVelos[vi]);
                if (d < 16.459 && mult <= mult54) ok = false;
                // CD must be positive
                if (dbg.CD_eff <= 0.0) ok = false;

                std::printf("%.3fm  %-4s  %-5s  %.4f  %s\n",
                            d, spinLabel[si], veloLabel[vi], mult,
                            ok ? "PASS" : "FAIL");
                ok ? ++pass : ++fail;
            }
        }
    }

    std::printf("\n%d passed, %d failed\n", pass, fail);
}
#endif // KIT_RUN_TESTS

// ------------------------------------------------------------
// Legacy lookup table (kept for reference / calibration only)
// Remove once physics model is validated in production.
// ------------------------------------------------------------
#if 0
struct MultRow { double dist_m; std::array<double,8> m; };

inline constexpr std::array<int,9> kSpinBins = {
    1000,1200,1400,1600,1800,2000,2200,2400,2600
};

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
    return 7;
}

inline double nearestDistanceKey(double dist_m) {
    double best = kMultiplier[0].dist_m;
    double bestDiff = std::abs(dist_m - best);
    for (const auto& row : kMultiplier) {
        double d = std::abs(dist_m - row.dist_m);
        if (d < bestDiff) { best = row.dist_m; bestDiff = d; }
    }
    return best;
}

inline double plateMultiplier(double dist_m, int spin_rpm) {
    const double dkey = nearestDistanceKey(dist_m);
    const int idx  = spinIndex(spin_rpm);
    for (const auto& row : kMultiplier) {
        if (std::abs(row.dist_m - dkey) < 0.001) return row.m[idx];
    }
    return 0.914;
}
#endif // legacy table

} // namespace kit

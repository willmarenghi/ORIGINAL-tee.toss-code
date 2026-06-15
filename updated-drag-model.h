#pragma once
#include <algorithm>
#include <cmath>

namespace kit {

// ------------------------------------------------------------
// Bat coefficients
// ------------------------------------------------------------
struct Coeff {
    double e;            // collision efficiency (ball-bat COR-like term)
    double k_potential;  // effective bat-efficiency factor for "potential" line
    double c;            // intercept (m/s); keep 0 for MVP
};

inline constexpr Coeff kMetal{ /*e*/0.28, /*k_potential*/0.92, /*c*/0.0 };
inline constexpr Coeff kWood { /*e*/0.23, /*k_potential*/0.92, /*c*/0.0 };

// ACTUAL MODEL NUMBERS
//
// Tee Mode:
//     Wood e_value = 0.22,  Wood K_potential = 1.21
//     Metal e_value = 0.26, Metal K_potential = 1.10
//     plateVelo = 0 m/s
//
// Toss Mode:
//     Wood e_value = 0.22,  Wood k_potential = 1.18
//     Metal e_value = 0.26, Metal k_potential = 1.10
//     plateVelo = 5.36 m/s (12 mph)

// ------------------------------------------------------------
// Physics-based drag model
//
// Formula:
//   omega      = spin_rpm * 2*pi/60                    (rad/s)
//   S          = r * omega / v_release                 (spin parameter, dimensionless)
//   CD0        = cd0ForVelocityAndDistance(v, dist_m)  (2D baseline drag coefficient)
//   CD_eff     = CD0 - kCDspin * S                    (effective drag after spin correction)
//   k_eff      = kAeroFactor * CD_eff                 (drag decay constant, 1/m)
//   multiplier = exp(-k_eff * dist_m)                 (fraction of release velocity reaching plate)
//
// CD0 is a 2D lookup: velocity x distance.
//   Slow anchor (40 mph): calibrated from 82 Stalker+KIT pitches (June 2026).
//                         RMS error = 1.69 mph across 25-50 ft.
//   Fast anchor (90 mph): back-calculated so 90mph/2200rpm/60.5ft = 10% loss
//                         (MLB standard). CD0_fast = 0.5199 flat across distances.
//   Between anchors:      linearly interpolated by release velocity.
//   Below 40 mph:         clamped to slow anchor values.
//   Above 90 mph:         clamped to fast anchor values.
//
// Why two anchors? Drag coefficient drops as velocity rises (Reynolds number /
// drag crisis effect). A baseball at 40 mph experiences ~2.6x more drag than at 90 mph.
//
// Supported distance range: 25-50 ft (7.62-15.24 m)
// Spin range: any RPM
// ------------------------------------------------------------
namespace drag {
    constexpr double kBallRadius_m = 0.03683;
    constexpr double kAirDensity   = 1.225;
    constexpr double kBallArea_m2  = 3.14159265358979 * kBallRadius_m * kBallRadius_m;
    constexpr double kBallMass_kg  = 0.1417;
    constexpr double kAeroFactor   = kAirDensity * kBallArea_m2 / (2.0 * kBallMass_kg);

    // Spin correction coefficient.
    // CD_eff = CD0(v, dist) - kCDspin * S,  where S = r*omega/v
    // Calibrated via least-squares regression on 82 Stalker+KIT pitches (June 2026).
    constexpr double kCDspin = 0.9946;

    // Velocity anchors for 2D CD0 interpolation
    constexpr double kVelSlow_mps = 17.882;  // 40 mph — slow-toss anchor
    constexpr double kVelFast_mps = 40.234;  // 90 mph — fast-pitch anchor

    // 2D CD0 calibration table: 7 distances x 2 velocity anchors.
    // cd0_slow: calibrated from real radar data (Stalker + KIT, June 2026).
    // cd0_fast: 0.5199 at all distances — derived from MLB 10% velocity loss standard.
    struct CD0Point { double dist_m; double cd0_slow; double cd0_fast; };
    inline constexpr CD0Point kCD0Table[] = {
        {  7.620,  0.7366,  0.5199 },  // 25ft — 14 pitches
        {  9.144,  0.8101,  0.5199 },  // 30ft — 15 pitches
        { 10.668,  0.7994,  0.5199 },  // 35ft — 15 pitches
        { 12.192,  0.8339,  0.5199 },  // 40ft — 10 pitches
        { 13.411,  0.7965,  0.5199 },  // 44ft — 10 pitches
        { 14.630,  0.7557,  0.5199 },  // 48ft —  8 pitches
        { 15.240,  0.6644,  0.5199 },  // 50ft — 10 pitches
    };
    inline constexpr int kCD0TableSize = 7;
}

// Returns the 2D CD0 baseline drag coefficient.
// Step 1: interpolates along the distance axis.
// Step 2: blends between slow (40mph) and fast (90mph) velocity anchors.
inline double cd0ForVelocityAndDistance(double v_mps, double dist_m) {
    const auto* t = drag::kCD0Table;
    const int   n = drag::kCD0TableSize;

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

    const double v_clamped = std::max(drag::kVelSlow_mps,
                                      std::min(v_mps, drag::kVelFast_mps));
    const double vfrac = (v_clamped - drag::kVelSlow_mps)
                       / (drag::kVelFast_mps - drag::kVelSlow_mps);
    return cd0_slow + vfrac * (cd0_fast - cd0_slow);
}

// Debug telemetry from the drag calculation
struct DragDebug {
    double omega_rads;  // angular velocity (rad/s)
    double spin_param;  // S = r*omega/v (dimensionless)
    double CD0_used;    // 2D baseline drag before spin correction
    double CD_eff;      // effective drag coefficient after spin correction
    double k_eff;       // drag decay constant (1/m)
};

// Returns the plate velocity multiplier: fraction of release velocity reaching the plate.
// Pass a non-null dbg pointer to capture intermediate values for diagnostics.
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
    double batSpeed_mps;      // bat speed at impact (m/s)
    double pitchRelease_mps;  // release velocity (m/s)
    double distance_m;        // pitch distance (25-50 ft / 7.62-15.24 m)
    int    spin_rpm;          // pitch spin from radar (exact RPM)
    double EV_measured_mps;   // measured exit velocity (m/s)
};

struct Outputs {
    DragDebug drag;              // intermediate drag values for diagnostics
    double multiplier_used;      // drag multiplier (unitless, 0-1)
    double plateVelocity_mps;    // V_plate = V_release * multiplier

    double potentialEV_mps;      // EV_pot = k*(1+e)*BS + e*V_plate + c
    double potentialSmash;       // EV_pot / bat speed

    double smash_measured;       // measured EV / bat speed
    double squaredUp_pct_raw;    // 100 * measured EV / potential EV
    double squaredUp_pct_ui;     // capped at 100% for display
};

inline Outputs computePotential(const Inputs& in, const Coeff batCoeff) {
    Outputs out{};

    out.multiplier_used   = plateMultiplierAdjusted(in.distance_m, in.spin_rpm,
                                                     in.pitchRelease_mps, &out.drag);
    out.plateVelocity_mps = in.pitchRelease_mps * out.multiplier_used;

    const double BS     = std::max(1e-6, in.batSpeed_mps);
    out.potentialEV_mps = batCoeff.k_potential * (1.0 + batCoeff.e) * BS
                        + batCoeff.e * out.plateVelocity_mps
                        + batCoeff.c;

    out.potentialSmash    = out.potentialEV_mps / BS;
    out.smash_measured    = in.EV_measured_mps / BS;

    const double denom    = std::max(1e-6, out.potentialEV_mps);
    out.squaredUp_pct_raw = (in.EV_measured_mps / denom) * 100.0;
    out.squaredUp_pct_ui  = std::min(out.squaredUp_pct_raw, 100.0);

    return out;
}

} // namespace kit

// =============================================================================
// test_drag.cpp  —  KIT Drag Model Offline Test
// =============================================================================
// Compile:
//   g++ -std=c++17 -DKIT_RUN_TESTS test_drag.cpp -o test_drag
// Run:
//   ./test_drag        (Mac/Linux)
//   test_drag.exe      (Windows)
// =============================================================================

#include "updated-drag-model.h"
#include <cstdio>
#include <cmath>

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static double mps_to_mph(double mps) { return mps * 2.23694; }
static double ft_to_m(double ft)     { return ft  * 0.3048;  }

static void separator(char c = '-', int n = 78) {
    for (int i = 0; i < n; ++i) std::putchar(c);
    std::putchar('\n');
}

// -----------------------------------------------------------------------------
// TEST 1 — Full multiplier table: distance × spin, at 3 pitch velocities
//          Shows the big picture of how the formula behaves across the range.
// -----------------------------------------------------------------------------
static void test1_multiplier_grid() {
    separator('=');
    std::printf("TEST 1 — Multiplier Grid  (distance × spin RPM)\n");

    const double velocities_mps[] = { 11.18, 17.88, 26.82 }; // 25, 40, 60 mph
    const char*  velo_labels[]    = { "25 mph (slow lob)",
                                      "40 mph (medium toss)",
                                      "60 mph (fast overhand)" };

    const double dists_ft[]  = { 25, 30, 35, 40, 44, 48, 50 };
    const int    spins[]     = { 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600 };
    const int    nDist       = 7;
    const int    nSpin       = 9;

    for (int vi = 0; vi < 3; ++vi) {
        double v = velocities_mps[vi];
        separator();
        std::printf("  Pitch velocity: %s (%.2f m/s)\n", velo_labels[vi], v);
        separator();

        // Header row
        std::printf("  dist  |");
        for (int si = 0; si < nSpin; ++si)
            std::printf(" %4d  |", spins[si]);
        std::printf("\n");
        separator();

        for (int di = 0; di < nDist; ++di) {
            double d = ft_to_m(dists_ft[di]);
            std::printf("  %2.0fft  |", dists_ft[di]);
            for (int si = 0; si < nSpin; ++si) {
                double mult = kit::plateMultiplierAdjusted(d, spins[si], v);
                std::printf(" %.4f |", mult);
            }
            std::printf("\n");
        }
    }
    std::printf("\n");
}

// -----------------------------------------------------------------------------
// TEST 2 — Delta vs old lookup table
//          Side-by-side: old bin-snapped value vs new formula value.
//          Old table only covered 40-54 ft at bin-center spins.
// -----------------------------------------------------------------------------
static void test2_delta_vs_old_table() {
    separator('=');
    std::printf("TEST 2 — New Formula vs Old Lookup Table (delta)\n");
    std::printf("  Positive delta = new formula gives MORE velocity reaching plate\n");
    std::printf("  Negative delta = new formula gives LESS velocity reaching plate\n");
    separator();

    // Old table values exactly as they were
    struct OldRow {
        double dist_ft;
        double dist_m;
        // bin centers: 1100,1300,1500,1700,1900,2100,2300,2500
        double old_mult[8];
    };

    const OldRow oldTable[] = {
        { 40, 12.192, {0.969,0.970,0.972,0.974,0.975,0.976,0.978,0.980} },
        { 44, 13.411, {0.964,0.966,0.968,0.970,0.972,0.973,0.975,0.977} },
        { 46, 14.021, {0.956,0.959,0.961,0.962,0.964,0.965,0.967,0.969} },
        { 48, 14.630, {0.942,0.945,0.948,0.951,0.953,0.955,0.958,0.961} },
        { 50, 15.240, {0.927,0.931,0.935,0.939,0.942,0.945,0.948,0.952} },
        { 54, 16.459, {0.901,0.906,0.910,0.914,0.918,0.922,0.926,0.931} },
    };

    const int   binCenterSpins[] = {1100,1300,1500,1700,1900,2100,2300,2500};
    const double v_ref_mps       = 17.88; // 40 mph reference velocity

    std::printf("  %-6s | %-6s |  OLD    NEW    DELTA\n", "dist", "spin");
    separator();

    double sumSqErr = 0.0;
    int    nPts     = 0;

    for (const auto& row : oldTable) {
        for (int si = 0; si < 8; ++si) {
            double newMult = kit::plateMultiplierAdjusted(row.dist_m,
                                                           binCenterSpins[si],
                                                           v_ref_mps);
            double delta = newMult - row.old_mult[si];
            sumSqErr += delta * delta;
            ++nPts;

            // Flag large deltas
            const char* flag = (std::abs(delta) > 0.010) ? "  <<" : "";
            std::printf("  %2.0fft   | %4drpm | %.4f  %.4f  %+.4f%s\n",
                        row.dist_ft, binCenterSpins[si],
                        row.old_mult[si], newMult, delta, flag);
        }
        separator();
    }

    double rms = std::sqrt(sumSqErr / nPts);
    std::printf("  RMS error vs old table: %.4f (%.2f%%)\n\n", rms, rms * 100.0);
}

// -----------------------------------------------------------------------------
// TEST 3 — Velocity sensitivity
//          Same distance + spin, sweep pitch velocity from 10-70 mph.
//          Shows how much the formula changes with pitch speed.
// -----------------------------------------------------------------------------
static void test3_velocity_sensitivity() {
    separator('=');
    std::printf("TEST 3 — Velocity Sensitivity  (fixed dist=46ft, spin=1800rpm)\n");
    std::printf("  Shows multiplier change as pitch speed varies\n");
    separator();

    const double dist_m  = ft_to_m(46.0);
    const int    spin    = 1800;
    const double velos[] = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70 };

    std::printf("  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n",
                "mph", "m/s", "S", "CD_eff", "k_eff", "mult");
    separator();

    for (double mph : velos) {
        double v = mph / 2.23694;
        kit::DragDebug dbg{};
        double mult = kit::plateMultiplierAdjusted(dist_m, spin, v, &dbg);
        std::printf("  %-8.1f  %-8.2f  %-8.4f  %-8.4f  %-8.5f  %.4f\n",
                    mph, v, dbg.spin_param, dbg.CD_eff, dbg.k_eff, mult);
    }
    std::printf("\n");
}

// -----------------------------------------------------------------------------
// TEST 4 — Short distance extrapolation (25-39 ft)
//          The old model had NO data here. Shows what the new formula gives.
// -----------------------------------------------------------------------------
static void test4_short_distance() {
    separator('=');
    std::printf("TEST 4 — Short Distance Extrapolation  (25-39 ft, new territory)\n");
    std::printf("  Old model had no data here. These are pure physics formula outputs.\n");
    separator();

    const double shortDists_ft[] = { 25, 27, 29, 31, 33, 35, 37, 39 };
    const int    spins[]         = { 1000, 1500, 2000, 2500 };
    const double v_mps           = 17.88; // 40 mph

    std::printf("  %-6s |", "dist");
    for (int sp : spins)
        std::printf(" %4drpm |", sp);
    std::printf("\n");
    separator();

    for (double ft : shortDists_ft) {
        double d = ft_to_m(ft);
        std::printf("  %2.0fft   |", ft);
        for (int sp : spins) {
            double mult = kit::plateMultiplierAdjusted(d, sp, v_mps);
            std::printf("  %.4f |", mult);
        }
        std::printf("\n");
    }
    std::printf("\n");
}

// -----------------------------------------------------------------------------
// TEST 5 — Edge cases
//          Things that could break the formula or produce bad values.
// -----------------------------------------------------------------------------
static void test5_edge_cases() {
    separator('=');
    std::printf("TEST 5 — Edge Cases\n");
    separator();

    struct Case {
        const char* label;
        double dist_m;
        int    spin;
        double v_mps;
    };

    const Case cases[] = {
        { "Minimum distance (25ft)",        ft_to_m(25),  1800, 17.88 },
        { "Maximum distance (54ft)",        ft_to_m(54),  1800, 17.88 },
        { "Very slow pitch (10mph)",        ft_to_m(46),  1800, 10.0/2.23694 },
        { "Very fast pitch (80mph)",        ft_to_m(46),  1800, 80.0/2.23694 },
        { "Minimum spin (500rpm)",          ft_to_m(46),   500, 17.88 },
        { "Maximum spin (3500rpm)",         ft_to_m(46),  3500, 17.88 },
        { "Zero spin (0rpm)",               ft_to_m(46),     0, 17.88 },
        { "Very short + very slow (toss)",  ft_to_m(25),   800,  8.94 },
        { "Full overhand 54ft high spin",   ft_to_m(54),  2600, 26.82 },
        { "Near-zero velocity (1mph)",      ft_to_m(46),  1800,  1.0/2.23694 },
    };

    std::printf("  %-38s  %-8s  %-8s  %s\n", "case", "mult", "CD_eff", "valid?");
    separator();

    for (const auto& c : cases) {
        kit::DragDebug dbg{};
        double mult = kit::plateMultiplierAdjusted(c.dist_m, c.spin, c.v_mps, &dbg);
        bool valid = (mult > 0.0 && mult <= 1.0 && dbg.CD_eff > 0.0);
        std::printf("  %-38s  %.4f    %.4f    %s\n",
                    c.label, mult, dbg.CD_eff, valid ? "PASS" : "FAIL <<<");
    }
    std::printf("\n");
}

// -----------------------------------------------------------------------------
// TEST 6 — Plate velocity through full computePotential()
//          Runs the whole model end-to-end with realistic inputs.
// -----------------------------------------------------------------------------
static void test6_full_model() {
    separator('=');
    std::printf("TEST 6 — Full Model (computePotential) End-to-End\n");
    separator();

    struct Scenario {
        const char* label;
        double batSpeed_mph;
        double pitchRelease_mph;
        double dist_ft;
        int    spin_rpm;
        double EV_measured_mph;
    };

    const Scenario scenarios[] = {
        { "Youth tee toss  25ft slow",   55,  25,  25, 1200,  68 },
        { "Travel ball 46ft med",        65,  45,  46, 1800,  82 },
        { "HS overhand 54ft fast",       70,  60,  54, 2400,  91 },
        { "Soft toss 30ft",              60,  30,  30, 1000,  75 },
        { "Overhand 48ft high spin",     68,  55,  48, 2600,  88 },
    };

    for (const auto& s : scenarios) {
        kit::Inputs in{};
        in.batSpeed_mps     = s.batSpeed_mph     / 2.23694;
        in.pitchRelease_mps = s.pitchRelease_mph  / 2.23694;
        in.distance_m       = ft_to_m(s.dist_ft);
        in.spin_rpm         = s.spin_rpm;
        in.EV_measured_mps  = s.EV_measured_mph  / 2.23694;

        kit::Outputs out = kit::computePotential(in, kit::kMetal);

        std::printf("  %s\n", s.label);
        std::printf("    Input:   BS=%.1fmph  Vrel=%.1fmph  dist=%.0fft"
                    "  spin=%drpm  EV_meas=%.1fmph\n",
                    s.batSpeed_mph, s.pitchRelease_mph,
                    s.dist_ft, s.spin_rpm, s.EV_measured_mph);
        std::printf("    Drag:    S=%.4f  CD=%.4f  k=%.5f/m  mult=%.4f\n",
                    out.drag.spin_param, out.drag.CD_eff,
                    out.drag.k_eff, out.multiplier_used);
        std::printf("    Output:  Vplate=%.1fmph  EV_pot=%.1fmph"
                    "  SqUp=%.1f%%  SmashPot=%.3f\n\n",
                    out.plateVelocity_mps * 2.23694,
                    out.potentialEV_mps   * 2.23694,
                    out.squaredUp_pct_ui,
                    out.potentialSmash);
        separator();
    }
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main() {
    separator('=');
    std::printf("  KIT Baseball — Drag Model Offline Test Suite\n");
    separator('=');
    std::printf("\n");

    test1_multiplier_grid();
    test2_delta_vs_old_table();
    test3_velocity_sensitivity();
    test4_short_distance();
    test5_edge_cases();
    test6_full_model();

    return 0;
}

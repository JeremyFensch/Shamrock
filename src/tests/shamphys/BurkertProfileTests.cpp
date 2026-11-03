// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2026 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

/**
 * @file BurkertProfileTests.cpp
 * @author Jeremy Fensch (jeremy.fensch@ens-lyon.fr)
 * @brief Unit tests of the Burkert halo profile
 */

#include "shambase/aliases_float.hpp"
#include "shambase/constants.hpp"
#include "shambackends/vec.hpp"
#include "shamphys/BurkertProfile.hpp"
#include "shamtest/shamtest.hpp"
#include <cmath>

namespace {

    /// closed form of f(x)/x^3, without the small x series
    f64 closed_form_mass_shape_over_x3(f64 x) {
        return (std::log1p(x * x) + 2 * std::log1p(x) - 2 * std::atan(x)) / (x * x * x);
    }

    /// enclosed mass by Simpson integration of 4 pi r^2 rho(r) over [0, r]
    f64 integrated_enclosed_mass(f64 r, f64 rho0, f64 r0, u32 n) {
        auto integrand = [&](f64 s) {
            return 4 * shambase::constants::pi<f64>
                   * s * s * shamphys::burkert_density(s, rho0, r0);
        };
        f64 h   = r / n;
        f64 sum = integrand(0) + integrand(r);
        for (u32 i = 1; i < n; i++) {
            sum += ((i % 2 == 1) ? 4 : 2) * integrand(i * h);
        }
        return sum * h / 3;
    }

} // namespace

NEW_TEST(Unittest, "shamphys/BurkertProfile", 1) {

    const f64 rho0 = 2.5;
    const f64 r0   = 0.7;
    const f64 G    = 1.3;
    const f64 pi   = shambase::constants::pi<f64>;

    { // the small x series matches the closed form where the closed form is still accurate
        for (f64 x : {0.01, 0.03, 0.0499999, 0.05}) {
            REQUIRE_FLOAT_EQUAL_NAMED(
                "series vs closed form",
                shamphys::burkert_mass_shape_over_x3(x) / closed_form_mass_shape_over_x3(x),
                1.,
                1e-11);
        }

        REQUIRE_FLOAT_EQUAL(shamphys::burkert_mass_shape_over_x3(0.), 4. / 3., 1e-16);
    }

    { // enclosed mass vs numerical integration of the density
        for (f64 r : {1e-3, 0.02, 0.3, 0.7, 2.0, 15.0}) {
            f64 m_ref = integrated_enclosed_mass(r, rho0, r0, 20000);
            REQUIRE_FLOAT_EQUAL_NAMED(
                "M(r) vs integral of rho",
                shamphys::burkert_enclosed_mass(r, rho0, r0) / m_ref,
                1.,
                1e-10);
        }
    }

    { // acceleration
        // regular at the center
        f64_3 a0 = shamphys::burkert_acceleration(f64_3{0, 0, 0}, rho0, r0, G);
        REQUIRE_EQUAL(a0.x(), 0.);
        REQUIRE_EQUAL(a0.y(), 0.);
        REQUIRE_EQUAL(a0.z(), 0.);

        // homogeneous core near the center: a = -4/3 pi G rho0 r
        f64_3 dr_small = {1e-9, -2e-9, 3e-9};
        f64_3 a_small  = shamphys::burkert_acceleration(dr_small, rho0, r0, G);
        f64_3 a_core   = (-4. / 3. * pi * G * rho0) * dr_small;
        REQUIRE_FLOAT_EQUAL(sycl::length(a_small - a_core) / sycl::length(a_core), 0., 1e-8);

        // -G M(r) / r^2 pointing to the center
        for (f64 r : {0.01, 0.5, 3.0}) {
            f64_3 dir = f64_3{1, 2, -2} / 3.;
            f64_3 acc = shamphys::burkert_acceleration(r * dir, rho0, r0, G);
            f64 g_ref = G * shamphys::burkert_enclosed_mass(r, rho0, r0) / (r * r);
            REQUIRE_FLOAT_EQUAL(sycl::length(acc + g_ref * dir) / g_ref, 0., 1e-14);
        }
    }
}

// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2026 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

#pragma once

/**
 * @file BurkertProfile.hpp
 * @author Jeremy Fensch (jeremy.fensch@ens-lyon.fr)
 * @brief Burkert (1995) dark matter halo density profile, enclosed mass and gravitational
 * acceleration
 *
 * The density profile is
 * \f[
 *  \rho(r) = \frac{\rho_0 r_0^3}{(r + r_0)(r^2 + r_0^2)}
 * \f]
 * and, with \f$x = r / r_0\f$, the mass enclosed within the radius \f$r\f$ is
 * \f[
 *  M(r) = \pi \rho_0 r_0^3 f(x), \quad f(x) = \ln(1 + x^2) + 2 \ln(1 + x) - 2 \arctan(x)
 * \f]
 * so that the gravitational acceleration at \f$\mathbf{r}\f$ from the halo center is
 * \f[
 *  \mathbf{a} = - \frac{G M(r)}{r^3} \mathbf{r} = - \pi G \rho_0 \frac{f(x)}{x^3} \mathbf{r}
 * \f]
 * which is regular at the center (\f$f(x)/x^3 \to 4/3\f$).
 */

#include "shambase/constants.hpp"
#include "shambackends/math.hpp"
#include "shambackends/sycl.hpp"

namespace shamphys {

    /**
     * @brief \f$f(x)/x^3\f$ for the Burkert profile, with
     * \f$f(x) = \ln(1 + x^2) + 2 \ln(1 + x) - 2 \arctan(x)\f$, \f$x = r/r_0 \geq 0\f$
     *
     * The terms of \f$f\f$ cancel at small \f$x\f$ (\f$f(x) \sim 4 x^3 / 3\f$ while each term is
     * \f$\mathcal{O}(x)\f$), so below `x_series` the Taylor series
     * \f$f(x) = \sum_{k \geq 0} \left[\frac{4 x^{4k+3}}{4k+3} - \frac{x^{4k+4}}{k+1}\right]\f$
     * is used instead. With the terms up to \f$k = 2\f$ its relative truncation error is
     * \f$\sim x^{12}/5\f$, below \f$10^{-16}\f$ for \f$x < 0.05\f$, where the closed form loses
     * \f$\sim 1.5 \epsilon / x^2\f$ (\f$< 10^{-13}\f$ in double precision).
     */
    template<class T>
    inline T burkert_mass_shape_over_x3(T x) {
        constexpr T x_series = T(0.05);

        if (x < x_series) {
            T x4 = x * x * x * x;
            // 4/3 - x + x^4 (4/7 - x/2) + x^8 (4/11 - x/3)
            return T(4) / T(3) - x + x4 * (T(4) / T(7) - x / T(2))
                   + x4 * x4 * (T(4) / T(11) - x / T(3));
        }

        T f = sycl::log1p(x * x) + T(2) * sycl::log1p(x) - T(2) * sycl::atan(x);
        return f / (x * x * x);
    }

    /**
     * @brief Density of the Burkert profile
     *
     * @param r distance to the halo center
     * @param rho0 central density
     * @param r0 core radius
     */
    template<class T>
    inline T burkert_density(T r, T rho0, T r0) {
        return rho0 * r0 * r0 * r0 / ((r + r0) * (r * r + r0 * r0));
    }

    /**
     * @brief Mass of the Burkert profile enclosed within the radius r
     *
     * @param r distance to the halo center
     * @param rho0 central density
     * @param r0 core radius
     */
    template<class T>
    inline T burkert_enclosed_mass(T r, T rho0, T r0) {
        T x = r / r0;
        return shambase::constants::pi<T> * rho0 * r * r * r * burkert_mass_shape_over_x3(x);
    }

    /**
     * @brief Gravitational acceleration of the Burkert profile
     *
     * @param dr position relative to the halo center
     * @param rho0 central density
     * @param r0 core radius
     * @param G gravitational constant
     */
    template<class Tvec, class T>
    inline Tvec burkert_acceleration(Tvec dr, T rho0, T r0, T G) {
        T x = sycl::sqrt(sham::length2(dr)) / r0;
        return (-shambase::constants::pi<T> * G * rho0 * burkert_mass_shape_over_x3(x)) * dr;
    }

} // namespace shamphys

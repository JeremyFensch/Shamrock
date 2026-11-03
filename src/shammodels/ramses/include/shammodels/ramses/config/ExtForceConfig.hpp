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
 * @file ExtForceConfig.hpp
 * @author Léodasce Sewanou (leodasce.sewanou@ens-lyon.fr)
 * @author Timothée David--Cléris (tim.shamrock@proton.me)
 * @brief External forces config of the Ramses solver + json serialization/deserialization
 */

#include "shambase/exception.hpp"
#include "shambase/overloaded.hpp"
#include "shambackends/type_convert.hpp"
#include "shambackends/vec.hpp"
#include "shammodels/ramses/config/enum_ExtForceCoupling.hpp"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace shammodels::basegodunov {

    /**
     * @brief List of the external forces applied to the gas
     *
     * The accelerations of all the forces are summed (in the order of the list) before being
     * applied to the conservative state, with the time quadrature selected by `coupling`.
     */
    template<class Tvec>
    struct ExtForceConfig {
        using Tscal = shambase::VecComponent<Tvec>;

        /// Constant gravitational acceleration g
        struct ConstantGravity {
            Tvec g;
        };

        /**
         * @brief Gravitational acceleration of a static Burkert (1995) dark matter halo
         *
         * Density profile \f$\rho(r) = \rho_0 r_0^3 / ((r + r_0)(r^2 + r_0^2))\f$ around
         * `center`, see shamphys/BurkertProfile.hpp. The gravitational constant is the one of the
         * unit system of the solver.
         */
        struct Burkert {
            Tvec center;
            Tscal rho0;
            Tscal r0;
        };

        using Variant = std::variant<ConstantGravity, Burkert>;

        std::vector<Variant> ext_forces;

        /// Time quadrature of the external forces source term
        ExtForceCoupling coupling = ExtForceCoupling::Midpoint;

        inline void set_coupling(ExtForceCoupling mode) { coupling = mode; }

        inline void add_constant_gravity(Tvec g) { ext_forces.push_back(ConstantGravity{g}); }

        inline void add_burkert(Tvec center, Tscal rho0, Tscal r0) {
            if (!(rho0 >= 0)) {
                shambase::throw_with_loc<std::invalid_argument>(
                    "Burkert external force: rho0 must be >= 0");
            }
            if (!(r0 > 0)) {
                shambase::throw_with_loc<std::invalid_argument>(
                    "Burkert external force: r0 must be > 0");
            }
            ext_forces.push_back(Burkert{center, rho0, r0});
        }

        inline bool has_ext_forces() const { return !ext_forces.empty(); }

        /// Whether one of the forces depends on the cell positions
        inline bool requires_coordinates() const {
            for (const auto &var_force : ext_forces) {
                if (std::holds_alternative<Burkert>(var_force)) {
                    return true;
                }
            }
            return false;
        }
    };

    template<class Tvec>
    inline void to_json(nlohmann::json &j, const ExtForceConfig<Tvec> &p) {
        using EF = ExtForceConfig<Tvec>;

        nlohmann::json force_list = nlohmann::json::array();
        for (const auto &var_force : p.ext_forces) {
            std::visit(
                shambase::overloaded{
                    [&](const typename EF::ConstantGravity &cfg) {
                        force_list.push_back(
                            nlohmann::json{{"force_type", "constant_gravity"}, {"g", cfg.g}});
                    },
                    [&](const typename EF::Burkert &cfg) {
                        force_list.push_back(
                            nlohmann::json{
                                {"force_type", "burkert"},
                                {"center", cfg.center},
                                {"rho0", cfg.rho0},
                                {"r0", cfg.r0}});
                    },
                },
                var_force);
        }

        j = nlohmann::json{{"force_list", force_list}, {"coupling", p.coupling}};
    }

    template<class Tvec>
    inline void from_json(const nlohmann::json &j, ExtForceConfig<Tvec> &p) {
        p.ext_forces.clear();

        // configs dumped before the coupling was configurable used the midpoint rule
        p.coupling = ExtForceCoupling::Midpoint;
        if (j.contains("coupling")) {
            j.at("coupling").get_to(p.coupling);
        }

        for (const auto &jforce : j.at("force_list")) {
            const std::string force_type = jforce.at("force_type").get<std::string>();
            if (force_type == "constant_gravity") {
                p.add_constant_gravity(jforce.at("g").get<Tvec>());
            } else if (force_type == "burkert") {
                using Tscal = typename ExtForceConfig<Tvec>::Tscal;
                p.add_burkert(
                    jforce.at("center").get<Tvec>(),
                    jforce.at("rho0").get<Tscal>(),
                    jforce.at("r0").get<Tscal>());
            } else {
                shambase::throw_with_loc<std::runtime_error>(
                    "Invalid external force type: " + force_type);
            }
        }
    }

} // namespace shammodels::basegodunov

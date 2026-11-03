// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2026 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

/**
 * @file ExtForceConfigTests.cpp
 * @author Timothée David--Cléris (tim.shamrock@proton.me)
 * @brief Unit tests for the RAMSES external forces config JSON (de)serialization
 */

#include "shammodels/ramses/SolverConfig.hpp"
#include "shamtest/shamtest.hpp"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <variant>

namespace {

    using Config           = shammodels::basegodunov::SolverConfig<f64_3, i64_3>;
    using ExtForceCoupling = shammodels::basegodunov::ExtForceCoupling;
    using ConstantGravity  = Config::ExtForceConfig::ConstantGravity;
    using Burkert          = Config::ExtForceConfig::Burkert;

    /// REQUIRE_EQUAL needs a formattable type, so modes are compared through their json name
    std::string coupling_name(ExtForceCoupling mode) {
        return nlohmann::json(mode).get<std::string>();
    }

    Config json_roundtrip(const Config &in) {
        nlohmann::json j = in;
        return nlohmann::json::parse(j.dump(4)).get<Config>();
    }

} // namespace

NEW_TEST(Unittest, "shammodels/ramses/SolverConfig::ext_force_config", 1) {

    { // json names
        REQUIRE_EQUAL(coupling_name(ExtForceCoupling::Midpoint), std::string("midpoint"));
        REQUIRE_EQUAL(coupling_name(ExtForceCoupling::Trapezoidal), std::string("trapezoidal"));

        nlohmann::json j = "simpson";
        REQUIRE_EXCEPTION_THROW(j.get<ExtForceCoupling>(), std::runtime_error);
    }

    { // default is the midpoint rule without forces
        Config cfg{};
        REQUIRE(!cfg.has_ext_forces());
        REQUIRE_EQUAL(
            coupling_name(cfg.ext_force_config.coupling),
            coupling_name(ExtForceCoupling::Midpoint));
    }

    for (ExtForceCoupling mode : {ExtForceCoupling::Midpoint, ExtForceCoupling::Trapezoidal}) {
        // forces and coupling json roundtrip
        Config in{};
        in.ext_force_config.add_constant_gravity({0., -0.1, 0.});
        in.ext_force_config.set_coupling(mode);

        Config out = json_roundtrip(in);
        REQUIRE_EQUAL(coupling_name(out.ext_force_config.coupling), coupling_name(mode));
        REQUIRE_EQUAL(out.ext_force_config.ext_forces.size(), 1_u64);

        const ConstantGravity *g
            = std::get_if<ConstantGravity>(&out.ext_force_config.ext_forces[0]);
        REQUIRE(g != nullptr);
        if (g != nullptr) {
            REQUIRE_EQUAL(g->g.y(), -0.1);
        }
    }

    { // Burkert halo json roundtrip, it requires the cell coordinates
        Config in{};
        REQUIRE(!in.is_coordinate_field_required());
        in.ext_force_config.add_constant_gravity({0., -0.1, 0.});
        in.ext_force_config.add_burkert({0.5, 0.25, -1.}, 3., 0.2);
        REQUIRE(in.is_coordinate_field_required());

        Config out = json_roundtrip(in);
        REQUIRE(out.is_coordinate_field_required());
        REQUIRE_EQUAL(out.ext_force_config.ext_forces.size(), 2_u64);

        const Burkert *b = std::get_if<Burkert>(&out.ext_force_config.ext_forces[1]);
        REQUIRE(b != nullptr);
        if (b != nullptr) {
            REQUIRE_EQUAL(b->center.x(), 0.5);
            REQUIRE_EQUAL(b->center.y(), 0.25);
            REQUIRE_EQUAL(b->center.z(), -1.);
            REQUIRE_EQUAL(b->rho0, 3.);
            REQUIRE_EQUAL(b->r0, 0.2);
        }
    }

    { // invalid Burkert parameters
        Config cfg{};
        REQUIRE_EXCEPTION_THROW(
            cfg.ext_force_config.add_burkert({0., 0., 0.}, 1., 0.), std::invalid_argument);
        REQUIRE_EXCEPTION_THROW(
            cfg.ext_force_config.add_burkert({0., 0., 0.}, -1., 1.), std::invalid_argument);
    }

    { // configs dumped before the coupling was configurable fall back on the midpoint rule
        Config in{};
        in.ext_force_config.add_constant_gravity({0., -0.1, 0.});
        in.ext_force_config.set_coupling(ExtForceCoupling::Trapezoidal);

        nlohmann::json j = in;
        j.at("ext_force_config").erase("coupling");

        Config out = j.get<Config>();
        REQUIRE_EQUAL(
            coupling_name(out.ext_force_config.coupling),
            coupling_name(ExtForceCoupling::Midpoint));
        REQUIRE_EQUAL(out.ext_force_config.ext_forces.size(), 1_u64);
    }
}

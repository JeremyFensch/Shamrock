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
 * @file AddExtForceBurkert.hpp
 * @author Jeremy Fensch (jeremy.fensch@ens-lyon.fr)
 * @brief Node adding the gravitational acceleration of a static Burkert dark matter halo to the
 * external forces acceleration
 */

#include "shambackends/vec.hpp"
#include "shamrock/solvergraph/IFieldSpan.hpp"
#include "shamrock/solvergraph/Indexes.hpp"
#include "shamsolvergraph/node/INode.hpp"

#define NODE_EDGES(X_RO, X_RW)                                                                     \
    /* ------------------- inputs ------------------- */                                           \
    X_RO(shamrock::solvergraph::Indexes<u32>, sizes)                                               \
    X_RO(shamrock::solvergraph::IFieldSpan<Tvec>, spans_coordinates)                               \
                                                                                                   \
    /* ------------------- outputs ------------------- */                                          \
    X_RW(shamrock::solvergraph::IFieldSpan<Tvec>, spans_a_ext)

namespace shammodels::basegodunov::modules {

    /**
     * @brief Add the gravitational acceleration of a static Burkert (1995) halo, evaluated at the
     * cell centers, to the external forces acceleration
     *
     * See shamphys/BurkertProfile.hpp for the profile.
     */
    template<class Tvec>
    class AddExtForceBurkert : public shamrock::solvergraph::INode {

        using Tscal = shambase::VecComponent<Tvec>;

        u32 block_size;
        Tvec center;
        Tscal rho0;
        Tscal r0;
        Tscal G;

        public:
        AddExtForceBurkert(u32 block_size, Tvec center, Tscal rho0, Tscal r0, Tscal G)
            : block_size(block_size), center(center), rho0(rho0), r0(r0), G(G) {}

        EXPAND_NODE_EDGES(NODE_EDGES)

        void _impl_evaluate_internal();

        inline virtual std::string _impl_get_label() const { return "AddExtForceBurkert"; }

        virtual std::string _impl_get_tex() const;
    };

} // namespace shammodels::basegodunov::modules

#undef NODE_EDGES

// -------------------------------------------------------//
//
// SHAMROCK code for hydrodynamics
// Copyright (c) 2021-2026 Timothée David--Cléris <tim.shamrock@proton.me>
// SPDX-License-Identifier: CeCILL Free Software License Agreement v2.1
// Shamrock is licensed under the CeCILL 2.1 License, see LICENSE for more information
//
// -------------------------------------------------------//

/**
 * @file AddExtForceBurkert.cpp
 * @author Jeremy Fensch (jeremy.fensch@ens-lyon.fr)
 * @brief
 */

#include "shambase/string.hpp"
#include "sham/format/format.hpp"
#include "shambackends/kernel_call_distrib.hpp"
#include "shammodels/ramses/modules/AddExtForceBurkert.hpp"
#include "shamphys/BurkertProfile.hpp"
#include "shamrock/patch/PatchDataField.hpp"
#include "shamsys/NodeInstance.hpp"

namespace shammodels::basegodunov::modules {

    template<class Tvec>
    void AddExtForceBurkert<Tvec>::_impl_evaluate_internal() {
        auto edges = get_edges();

        edges.spans_coordinates.check_sizes(edges.sizes.indexes);
        edges.spans_a_ext.check_sizes(edges.sizes.indexes);

        shambase::DistributedData<u32> cell_counts
            = edges.sizes.indexes.template map<u32>([&](u64 id, u32 block_count) {
                  return block_count * block_size;
              });

        sham::distributed_data_kernel_call(
            shamsys::instance::get_compute_scheduler_ptr(),
            sham::DDMultiRef{edges.spans_coordinates.get_spans()},
            sham::DDMultiRef{edges.spans_a_ext.get_spans()},
            cell_counts,
            [center = center, rho0 = rho0, r0 = r0, G = G](
                u32 i, const Tvec *__restrict xyz, Tvec *__restrict a_ext) {
                a_ext[i] += shamphys::burkert_acceleration(xyz[i] - center, rho0, r0, G);
            });
    }

    template<class Tvec>
    std::string AddExtForceBurkert<Tvec>::_impl_get_tex() const {
        auto xyz   = get_ro_edge_base(1).get_tex_symbol();
        auto a_ext = get_rw_edge_base(0).get_tex_symbol();

        std::string tex = R"tex(
            Add the gravitational acceleration of a Burkert halo to the external forces acceleration

            \begin{align}
            \mathbf{r}_i &= {xyz}_i - \mathbf{r}_c, \quad x_i = \vert \mathbf{r}_i \vert / r_0 \\
            {a_ext}_i &= {a_ext}_i - \pi G \rho_0
                \frac{\ln(1 + x_i^2) + 2 \ln(1 + x_i) - 2 \arctan(x_i)}{x_i^3} \mathbf{r}_i \\
            \mathbf{r}_c &= ({cx}, {cy}, {cz}), \quad \rho_0 = {rho0}, \quad r_0 = {r0},
                \quad G = {G}
            \end{align}
        )tex";

        shambase::replace_all(tex, "{a_ext}", a_ext);
        shambase::replace_all(tex, "{xyz}", xyz);
        shambase::replace_all(tex, "{cx}", sham::format("{}", center.x()));
        shambase::replace_all(tex, "{cy}", sham::format("{}", center.y()));
        shambase::replace_all(tex, "{cz}", sham::format("{}", center.z()));
        shambase::replace_all(tex, "{rho0}", sham::format("{}", rho0));
        shambase::replace_all(tex, "{r0}", sham::format("{}", r0));
        shambase::replace_all(tex, "{G}", sham::format("{}", G));

        return tex;
    }

} // namespace shammodels::basegodunov::modules

template class shammodels::basegodunov::modules::AddExtForceBurkert<f64_3>;

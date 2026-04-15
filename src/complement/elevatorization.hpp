// Copyright (C) 2022  The Kofola Authors
//
// Kofola is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Kofola is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// kofola
#include "../util/helpers.hpp"
#include "../algorithms/abstract_complement_alg.hpp"

// spot
#include <spot/twa/twa.hh>
#include <spot/twaalgos/isdet.hh>


namespace kofola
{ // {{{

    // limit determinism macrostates
    struct RBL
    {
        std::set<unsigned> R{};
        std::set<unsigned> B{};
        unsigned l = 0;

        bool operator==(const RBL& other) const {
            return R == other.R && B == other.B && l == other.l;
        }

        bool operator<(const RBL& other) const {
            if (R != other.R) return R < other.R;
            if (B != other.B) return B < other.B;
            return l < other.l;
        }
    };

    inline std::ostream& operator<<(std::ostream& os, const RBL& x)
    {
        os << "R={";
        for (auto it = x.R.begin(); it != x.R.end(); ++it) {
            if (it != x.R.begin()) os << ",";
            os << *it;
        }
        os << "}, B={";
        for (auto it = x.B.begin(); it != x.B.end(); ++it) {
            if (it != x.B.begin()) os << ",";
            os << *it;
        }
        os << "}, l=" << x.l;
        return os;
    };

    class Elevatorization {
    public:
        Elevatorization(const spot::twa_graph_ptr& aut);

        spot::twa_graph_ptr elevatorize(bool only_non_buchi=false);

        void limit_deter(size_t part_index, size_t scc_idx);

        void create_deter_part(size_t scc_idx, kofola::AccClause disjunct);

        void remove_all_cols_within_scc(size_t scc_idx);

        void add_cols_within_scc(size_t part_index, spot::acc_cond::mark_t fins);
    
        std::set<unsigned> get_succ_excluding_colors(
            const std::set<unsigned>& states,
            const bdd& bdd,
            const spot::acc_cond::mark_t& col
        );

        std::set<unsigned> get_succ_including_colors(
            const std::set<unsigned>& states,
            const bdd& bdd,
            const spot::acc_cond::mark_t& col,
            const spot::acc_cond::mark_t& without_cols
        );

    private:
        std::unique_ptr<kofola::cmpl_info> info_;
        const spot::twa_graph_ptr& aut_;
        unsigned old_aut_num_states_;
        spot::scc_info si_;

        // Support for each state of the source automaton.
        std::vector<bdd> support_;

        // Propositions compatible with all transitions of a state.
        std::vector<bdd> compat_;

        // Whether a SCC is deterministic or not
        std::string scc_types_;

        std::tuple<size_t,
                kofola::PartitionToTypeMap,
                kofola::StateToPartitionMap,
                kofola::SCCToPartitionMap,
                kofola::PartitionToAccMap
        > partitions_;

        unsigned new_inf_col_;
        unsigned new_fin_col_;
    };

} // kofola }}}

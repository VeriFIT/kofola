// Copyright (C) 2017-2019 Laboratoire de Recherche et Développement
// de l'Epita.
// Copyright (C) 2022  The COLA Authors
//
// COLA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// COLA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

// kofola
#include "kofola.hpp"
#include "types.hpp"
#include "complement_tela.hpp"
#include "util.hpp"

#include "abstract_complement_alg.hpp"
#include "complement_alg_mh.hpp"
#include "complement_alg_ncsb.hpp"
#include "complement_alg_ncsb_delay.hpp"
#include "complement_alg_safra.hpp"
// #include "complement_alg_rank.hpp"
#include "complement_alg_rank2.hpp"
#include "complement_alg_init_det.hpp"
#include "complement_alg_subs_tuple.hpp"

#include <deque>
#include <map>
#include <set>
#include <stack>
#include <queue>

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/dot.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/cleanacc.hh>
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/misc/bddlt.hh>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/misc/version.hh>
#include <spot/twa/acc.hh>

namespace cola
{
    // complementation Buchi automata
    class tnba_complement
    {
    public:
        /// target of a transition (including colours)
        using state_col = std::pair<unsigned, std::set<unsigned>>;
        /// target of a transition (including colours tagged by partition index)
        using state_taggedcol = std::pair<unsigned, std::set<std::pair<unsigned, unsigned>>>;

        /// return type of get_succ (set of pairs <state, set of colours>)
        using vec_state_col = std::vector<state_col>;
        /// vec_state_col where colours are tagged by their partition
        using vec_state_taggedcol = std::vector<state_taggedcol>;

        // new interface
        using abs_cmpl_alg_p = std::unique_ptr<kofola::abstract_complement_alg>;
        using vec_algorithms = std::vector<abs_cmpl_alg_p>;

        using abs_cmpl_ms_p = std::shared_ptr<kofola::abstract_complement_alg::mstate>;
        using vec_macrostates = std::vector<abs_cmpl_ms_p>;

    private:
        // The source automaton.
        spot::const_twa_graph_ptr aut_;

        // Direct simulation on source automaton.
        kofola::Simulation dir_sim_;

        // vector of reachable states for every state
        kofola::ReachableVector reachable_vector_;

        // SCCs information of the source automaton.
        spot::scc_info si_;

        // general info for complementation
        std::unique_ptr<kofola::cmpl_info> info_;

        // vector of algorithms for complementation
        vec_algorithms alg_vec_;

        //
        size_t num_colours_;

        //
        spot::acc_cond::acc_code final_code_;

        //
        std::vector<spot::acc_cond> vec_acc_code_;

        //
        std::map<unsigned, unsigned> part_col_offset_;

        //
        spot::acc_cond::acc_code sink_acc_code_;

        std::tuple<size_t,
                kofola::PartitionToTypeMap,
                kofola::StateToPartitionMap,
                kofola::SCCToPartitionMap
        > partitions_;

        bool is_sink_created_ = false;
        unsigned sink_state_ = static_cast<unsigned>(-1);

        // Number of states in the input automaton.
        unsigned nb_states_;

        // state_simulator
        // state_simulator simulator_;

        // delayed simulation
        // delayed_simulation delayed_simulator_;

        // The parity automata being built.
        spot::twa_graph_ptr res_;

        // the number of indices
        unsigned sets_ = 0;

        unsigned num_colors_;

        // Association between labelling states and state numbers of the
        // DPA.
        // std::unordered_map<complement_mstate, unsigned, complement_mstate_hash> rank2n_;

        // States to process.
        // std::deque<std::pair<complement_mstate, unsigned>> todo_;

        // Support for each state of the source automaton.
        std::vector<bdd> support_;

        // Propositions compatible with all transitions of a state.
        std::vector<bdd> compat_;

        // is accepting for states
        std::vector<bool> is_accepting_;

        // Whether a SCC is deterministic or not
        std::string scc_types_;

        // State names for graphviz display
        std::vector<std::string>* names_;

        // the index of each weak SCCs
        std::vector<unsigned> weaksccs_;
        // the index of each deterministic accepting SCCs
        std::vector<unsigned> acc_detsccs_;
        // the index of each deterministic accepting SCCs
        std::vector<unsigned> acc_nondetsccs_;

        // Show Rank states in state name to help debug
        bool show_names_;

        std::map<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned> rank_bounds_; // TODO

        std::string get_det_string(const std::vector<state_rank> &states);


    public:
        tnba_complement(const spot::twa_graph_ptr &aut, spot::scc_info& si);

        unsigned get_num_states();

        spot::scc_info & get_scc_info();

        void reduce_and_compute_simulation();

        std::set<int> reachable_vertices(std::vector<std::vector<int>> list, std::set<int> from);

        std::vector<std::set<int>> get_reachable_vector();

        // ######################################################################
        // NEW INTERFACE
        // ######################################################################

        /*using abs_cmpl_alg_p = std::unique_ptr<kofola::abstract_complement_alg>;
        using vec_algorithms = std::vector<abs_cmpl_alg_p>;

        using abs_cmpl_ms_p = std::shared_ptr<kofola::abstract_complement_alg::mstate>;
        using vec_macrostates = std::vector<abs_cmpl_ms_p>;*/

        /// the uberstate - combination of all partial macrostates
        class uberstate
        { // {{{
        private:  // DATA MEMBERS

            /// all reached states
            std::set<unsigned> reached_states_;
            /// vector of partial macrostates
            vec_macrostates part_macrostates_;
            /// index of the active component for round robin
            int active_scc_;

            /// shared breakpoint
            std::set<unsigned> shared_breakpoint_;

        public:  // METHODS

            /// constructor
            uberstate(const std::set<unsigned>& reached_states,
                      const vec_macrostates& part_macrostates,
                      int active_scc,
                      const std::set<unsigned>& shared_breakpoint);

            /// move constructor
            uberstate(uberstate&& us) = default;

            /// deleted copy constructor and assignment operator
            uberstate(const uberstate& us);

            uberstate& operator=(const uberstate& us) = delete;

            /// converts to string
            std::string to_string() const;

            /// output stream operator
            friend std::ostream& operator<<(std::ostream& os, const uberstate& us)
            {
                os << us.to_string();
                return os;
            }

            /// returns the set of all reached states
            const std::set<unsigned>& get_reach_set() const;

            /// returns the partial macrostates
            const vec_macrostates& get_part_macrostates() const;

            /// returns the index of the active SCC (INACTIVE_SCC if no active)
            const int get_active_scc() const;

            /// get shared breakpoint
            const std::set<unsigned>& get_shared_breakpoint() const;

            /// total ordering operator to allow use in std::set and std::map
            bool operator<(const uberstate& rhs) const;

            bool operator==(const uberstate& rhs) const;

        }; // uberstate }}}

        /// functor for comparison of uberstate pointers
        struct uberstate_ptr_less_ftor
        {
            bool operator()(const uberstate* lhs, const uberstate* rhs) const
            {
                assert(lhs && rhs);
                return *lhs < *rhs;
            }
        };


        // Here we have a bidirectional map between uberstates and state
        // identifiers (unsigned).  The uberstates are physically stored only at
        // 'num_to_uberstate_map_', the reason being that they contain vectors of
        // unique_ptr (no copy is therefore allowed).

        /// maps uberstates to state numbers
        std::map<const uberstate*, unsigned, uberstate_ptr_less_ftor> uberstate_to_num_map_;
        /// maps state numbers to uberstates
        std::vector<std::shared_ptr<uberstate>> num_to_uberstate_map_;
        /// counter of states (to be assigned to uberstates) - 0 is reserved for sink
        unsigned cnt_state_ = 0;

        // reserved colours
        // static const unsigned SINK_COLOUR = 0;
        enum {SINK_COLOUR = 0};
        // static const unsigned RR_COLOUR = 1;      // colour for round robin ### we will choose it dynamically
        static const size_t RESERVED_COLOURS = 1; // how many colours are reserved

        /// index of active SCC if no SCC is active
        static const int INACTIVE_SCC = -1;

        /// accessor into the uberstate table
        unsigned uberstate_to_num(const uberstate& us) const;

        /// translates state number to uberstate
        const uberstate& num_to_uberstate(unsigned num) const;

        /// inserts an uberstate (by moving) and returns its assigned number (if
        /// not present), or just returns the number of an equal uberstate (if
        /// present)
        unsigned insert_uberstate(const uberstate& us);


        /// computes the Cartesian product of a vector of sets (no repetitions
        /// assumed in the inputs)
        template<class A>
        std::vector<std::vector<A>> compute_cartesian_prod(
                const std::vector<std::vector<A>> vec_of_sets){ // {{{
            const size_t length = vec_of_sets.size();
            std::vector<std::vector<A>> result;

            // this vector will iterate over all possible tuples of indices
            std::vector<size_t> indices(length, 0);

            while (true) {
                std::vector<A> vec;
                for (size_t i = 0; i < length; ++i) {
                    assert(indices[i] < vec_of_sets[i].size());
                    vec.push_back(vec_of_sets[i][indices[i]]);
                }

                assert(vec.size() == length);
                result.push_back(std::move(vec));

                // generate the next vector of indices, if possible
                bool generated = false;
                for (size_t j = 0; j < length; ++j) {
                    ++(indices[j]);
                    if (indices[j] < vec_of_sets[j].size()) { // indices is set
                        generated = true;
                        break;
                    } else { // we need to move into the next index
                        indices[j] = 0;
                    }
                }

                if (!generated) { break; }
            }

            return result;
        } // compute_cartesian_prod() }}}

        /// removes duplicit values (warning: can change order!)
        template <class T>
        static void remove_duplicit(T& t){ // {{{
            std::sort(t.begin(), t.end());
            t.erase(std::unique(t.begin(), t.end()), t.end());
        } // remove_duplicit }}}

        static int get_next_active_scc(const vec_algorithms& alg_vec, int prev);


        static kofola::PartitionToSCCMap create_part_to_scc_map(
                const kofola::SCCToPartitionMap& scc_to_part_map);

        static kofola::SCCToSCCSetMap create_scc_to_pred_sccs_map(
                const spot::scc_info&           si,
                const kofola::ReachableVector&  reach_vec);


        /// gets all successors of an uberstate wrt a vector of algorithms and a
        /// symbol
        vec_state_taggedcol get_succ_uberstates(
                const uberstate&       src,
                const bdd&             symbol);

        /// gets all initial uberstates wrt a vector of algorithms
        std::vector<unsigned> get_initial_uberstates();


        /// partitions the SCCs of the input automaton according to decomposition
        /// options, returns a triple (num_partitions, partition_types,
        /// state_to_partition_map, scc_to_partition_map)
        static std::tuple<size_t,
                kofola::PartitionToTypeMap,
                kofola::StateToPartitionMap,
                kofola::SCCToPartitionMap
        >
        create_partitions(
                const spot::scc_info&   scc_inf,
                const kofola::options&  options);


        /// selects the algorithms to run on the SCCs
        void select_algorithms();

        ///
        unsigned int get_cnt_state_();

        void inc_cnt_state_();

        ///
        bdd get_support_at(unsigned s);

        ///
        bdd get_compat_at(unsigned s);

        ///
        void handle_sink_state();

        ///
        bool get_is_sink_created();

        ///
        unsigned get_sink_state();

        /// new modular complementation procedure
        spot::twa_graph_ptr run_new();

        ///
        void prep_for_compl();

        ///
        std::vector<spot::acc_cond> get_vec_acc_cond();

        ///
        spot::acc_cond::acc_code get_final_acc_code();

        ///
        unsigned int get_alg_vec_mincolour_at_i(unsigned i);

        ///
        std::map<unsigned int, unsigned int> get_part_col_offset();

        ///
        size_t  set_acc_cond();
    };
}

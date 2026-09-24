#pragma once

// spot
#include <spot/twa/twa.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/sccinfo.hh>

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "../util/helpers.hpp"
#include "abstract_determinize_alg.hpp"

namespace kofola {
    using scc_partitions_t = std::tuple<size_t,
            kofola::PartitionToTypeMap,
            kofola::StateToPartitionMap,
            kofola::SCCToPartitionMap,
            kofola::PartitionToAccMap>;

    /// Determinize a transition-based Emerson-Lei automaton (TELA).
    ///
    /// This is currently a thin wrapper around `tela_determinize`.
    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut);

    /// Modular determinization of transition-based Emerson-Lei automata (TELA).
    ///
    /// The construction mirrors the modular complementation
    /// (@p helpers::tnba_complement) but is considerably simpler: the SCCs of
    /// the input automaton are split into partitions, a (deterministic) partial
    /// determinization algorithm is run on each of them synchronously over a
    /// shared subset construction, and the acceptance condition of the result is
    /// the *disjunction* of the (colour-shifted) partial conditions.  Since
    /// L(A) is the union of the languages of the runs eventually staying in one
    /// partition, this yields a deterministic automaton for L(A).
    class tela_determinize
    {
    public: // TYPES

        using alg_p = std::unique_ptr<kofola::abstract_determinize_alg>;
        using vec_algorithms = std::vector<alg_p>;

        using ms_p = std::shared_ptr<kofola::abstract_determinize_alg::mstate>;
        using vec_macrostates = std::vector<ms_p>;

        /// target of a transition (including colours tagged by partition index)
        using state_taggedcol = std::pair<unsigned, std::set<std::pair<unsigned, unsigned>>>;

        /// the uberstate - combination of all partial macrostates
        class uberstate
        { // {{{
        private:  // DATA MEMBERS

            /// all reached states
            std::set<unsigned> reached_states_;
            /// vector of partial macrostates
            vec_macrostates part_macrostates_;

        public:  // METHODS

            /// constructor
            uberstate(const std::set<unsigned>& reached_states,
                      const vec_macrostates& part_macrostates);

            /// copy and move constructors
            uberstate(const uberstate& us) = default;
            uberstate(uberstate&& us) = default;

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

    private: // DATA MEMBERS

        spot::twa_graph_ptr aut_;
        std::shared_ptr<spot::scc_info> scc_;
        kofola::scc_partitions_t partitions_;

        /// general info about the automaton (shared with complementation)
        std::unique_ptr<kofola::cmpl_info> info_;

        /// vector of partial determinization algorithms (one per partition)
        vec_algorithms alg_vec_;

        /// map of partitions to sets of SCCs they contain (referenced by info_)
        kofola::PartitionToSCCMap part_to_scc_map_;
        /// empty auxiliary structures referenced by info_
        kofola::ReachableVector reachable_vector_;
        kofola::SCCToSCCSetMap scc_to_pred_sccs_map_;
        kofola::Simulation dir_sim_;

        /// support of each state of the input automaton
        std::vector<bdd> support_;
        /// propositions compatible with all transitions of a state
        std::vector<bdd> compat_;
        /// is accepting for states
        std::vector<bool> is_accepting_;

        /// number of colours of the result
        size_t num_colours_;
        /// acceptance condition of the result
        spot::acc_cond::acc_code final_code_;
        /// acceptance conditions of the partial algorithms
        std::vector<spot::acc_cond> vec_acc_code_;
        /// offset of the colour block of each partition
        std::map<unsigned, unsigned> part_col_offset_;

        // bidirectional map between uberstates and state identifiers
        /// maps uberstates to state numbers
        std::map<const uberstate*, unsigned, uberstate_ptr_less_ftor> uberstate_to_num_map_;
        /// maps state numbers to uberstates
        std::vector<std::shared_ptr<uberstate>> num_to_uberstate_map_;
        /// counter of states (to be assigned to uberstates)
        unsigned cnt_state_;

        /// show macrostate labels in the output
        bool show_names_;

    public: // METHODS

        tela_determinize(const spot::twa_graph_ptr& aut,
                         std::shared_ptr<spot::scc_info> scc,
                         kofola::scc_partitions_t partitions);

        /// modular determinization procedure
        spot::twa_graph_ptr run_new();

    private: // METHODS

        /// selects the algorithms to run on the partitions
        void select_algorithms();

        /// creates algorithm for an inherently weak partition
        alg_p create_inherently_weak_algorithm(size_t partition_index);
        /// creates algorithm for a deterministic partition
        alg_p create_deterministic_algorithm(size_t partition_index);
        /// creates algorithm for a strongly deterministic partition
        alg_p create_strongly_deterministic_algorithm(size_t partition_index);
        /// creates algorithm for a nondeterministic partition
        alg_p create_nondeterministic_algorithm(size_t partition_index);
        /// creates algorithm for an initial almost deterministic partition
        alg_p create_initial_almost_deterministic_algorithm(size_t partition_index);

        /// accessor into the uberstate table
        const uberstate& num_to_uberstate(unsigned num) const;

        /// inserts an uberstate and returns its assigned number (if not present),
        /// or just returns the number of an equal uberstate (if present)
        unsigned insert_uberstate(const uberstate& us);

        /// gets the initial uberstate
        unsigned get_initial_uberstate();

        /// gets the (unique) successor of an uberstate over a symbol
        state_taggedcol get_succ_uberstate(const uberstate& src, const bdd& symbol);

        /// computes the acceptance condition of the result
        void set_acc_cond();

        /// returns an automaton with the empty language
        spot::twa_graph_ptr make_empty_aut() const;
    };
}

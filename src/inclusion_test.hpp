#pragma once

// kofola
#include "kofola.hpp"
#include "abstract_complement_alg.hpp"
#include "complement_sync.hpp"

// spot
#include <spot/twa/twa.hh>

spot::formula formula();

namespace kofola
{
    class inclusionTest {
    public:
        inclusionTest();

        using intersect_mstate = std::pair<unsigned, unsigned>; // (A,B)
        /// target of a transition (including colours tagged by partition index)
        using state_col = std::pair<intersect_mstate, std::set<unsigned>>;
        /// vec_state_col where colours are tagged by their partition
        using vec_state_col = std::vector<state_col>;

        /// true if A ⊆ B holds, false otherwise (calls tarjan_is_empty)
        bool test(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B);

    private:
        const int UNDEFINED = -1;

        spot::twa_graph_ptr preprocessed_orig_aut_B_;

        /// states of the inspected resulting automaton
        std::map<intersect_mstate , vec_state_col> intersect_states_;
        /// tarjan variables
        std::map<intersect_mstate, signed> indices_;
        std::map<intersect_mstate, signed> lowlinks_;
        std::map<intersect_mstate, bool> on_stack_;
        signed index_ = 0;
        std::stack<intersect_mstate> tarjan_stack_;
        /// end of tarjan variables

        /// stores which colors can be cycled from resulting macrostate
        std::map<intersect_mstate, std::set<unsigned>> cols_visited_;
        /// stores position of macrostate in the current dfs traversal
        std::map<intersect_mstate, unsigned> dfs_state_pos_;
        /// stores colors of strongly connected components
        std::vector<std::set<unsigned>> scc_cols_;

        /// stores accepting condition
        std::set<unsigned> inf_acc_cols_conj_;
        spot::acc_cond::acc_code acceptance_cond_;

        /// to stop searching when counter-example
        bool decided_ = false;
        bool empty_ = false;

        /// to preprocess autB
        /// TODO just copypaste from complement_sync.cpp
        spot::twa_graph_ptr preprocess(const spot::twa_graph_ptr &aut);

        /// returns "alphabet" of automaton aut_A
        std::pair<bdd, bdd> symbols_from_A(const spot::twa_graph_ptr &aut_A);

        /// returns set of all successors for given intersect_mstate
        std::set<state_col>
        get_all_succs(const spot::twa_graph_ptr &aut_A, cola::tnba_complement &aut_B, intersect_mstate &mstate,
                      std::map<intersect_mstate, vec_state_col> &intersect_states);

        ///
        bool is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol);

        /// returns cross product: states_A x states_B
        std::set<state_col>
        get_cross_prod(const spot::twa_graph_ptr &aut_A, unsigned aut_A_src, cola::tnba_complement &aut_B,
                       std::set<unsigned> &states_A, cola::tnba_complement::vec_state_taggedcol &states_B,
                       std::map<intersect_mstate, vec_state_col> &intersect_states, const bdd &letter);

        /// performs emptiness check: whether aut_A⊆aut_B using tarjan's algo
        void tarjan_is_empty(intersect_mstate &src_mstate, const spot::twa_graph_ptr &aut_A, cola::tnba_complement &aut_B, unsigned ith_on_path, std::vector<std::set<unsigned>> cols_path);
    };
}


#pragma once

// kofola
#include "emptiness_check.hpp"
#include "kofola.hpp"
#include "abstract_successor.hpp"

namespace kofola {
    class  inclusion_check;
    using intersect_mstate = std::pair<unsigned, unsigned>;

    class  inclusion_mstate : public abstr_succ::abstract_successor::mstate {
    private:
        intersect_mstate state_; /// systems, last state of aut
    public:
        inclusion_mstate() {
            ;
        }

        bool eq(const inclusion_mstate& rhs) {
            return (state_ == rhs.state_ && acc_ == rhs.acc_);
        }

        bool neq(const inclusion_mstate& rh) {
            return !eq(rh);
        }

        bool lt(const inclusion_mstate& rhs) {
            if(state_ != rhs.state_) {return state_ < rhs.state_;}
        }

        friend class inclusion_check;
    };

    class inclusion_check : public abstr_succ::abstract_successor {
        /// target of a transition (including colours tagged by partition index)
        using state_col = std::pair<intersect_mstate, std::set<unsigned>>;
        /// vec_state_col where colours are tagged by their partition
        using vec_state_col = std::vector<state_col>;
    private:
        spot::twa_graph_ptr preprocessed_orig_aut_B_;
        std::vector<std::unique_ptr<kofola::inclusion_mstate>> init_states_;
        std::map<intersect_mstate , vec_state_col> intersect_states_;

        spot::twa_graph_ptr aut_A_;
        cola::tnba_complement aut_B_compl_;
        spot::acc_cond::acc_code acc_cond_;

    public:
        inclusion_check(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B);

        cola::tnba_complement init_compl_aut_b(const spot::twa_graph_ptr &aut_B);

        bool inclusion();

        /// to preprocess autB
        /// TODO just copypaste from complement_sync.cpp
        spot::twa_graph_ptr preprocess(const spot::twa_graph_ptr &aut);

        /// returns "alphabet" of automaton aut_A
        std::pair<bdd, bdd> symbols_from_A(const spot::twa_graph_ptr &aut_A);

        std::vector<std::unique_ptr<abstract_successor::mstate>> get_initial_states() override;

        /// returns set of all successors for given intersect_mstate
        std::vector<std::unique_ptr<abstract_successor::mstate>> get_succs(const std::unique_ptr<abstract_successor::mstate> &src) override;

        bool is_accepting(spot::acc_cond::mark_t inf_cond) override;

        ///
        bool is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol);

        /// returns cross product: states_A x states_B
        std::vector<std::unique_ptr<inclusion_mstate>>
        get_cartesian_prod(unsigned aut_A_src, std::set<unsigned> &states_A,
                           cola::tnba_complement::vec_state_taggedcol &states_B,
                           const bdd &letter);

    };
}

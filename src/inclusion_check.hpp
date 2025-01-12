/**
 * @file inclusion_check.hpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief declarations for on the fly inclusion check
 * @version 0.1
 * @date 2024-05-03
 * 
 * 
 */

#pragma once

// kofola
#include "emptiness_check.hpp"
#include "kofola.hpp"

namespace kofola {
    class  inclusion_check;
    /// (aut_A,aut_B)
    using intersect_mstate = std::pair<unsigned, unsigned>;

    /// class for macrostate used within inclusion procedure
    class  inclusion_mstate {
    private:
        intersect_mstate state_;
        spot::acc_cond::mark_t acc_; /// acc. marks on transition incoming to the state
        bdd trans_cond_; /// transition condition to come to this state
        bool encountered_ = false;
    public:
        inclusion_mstate() {
            ;
        }

        spot::acc_cond::mark_t get_acc() {return acc_; }
        void set_acc(spot::acc_cond::mark_t new_acc) {acc_ = new_acc; }

        void set_encountered(bool val) {encountered_ = val;}
        bool get_encountered() {return encountered_;}

        /// equality of inclusion macrostates
        bool eq(const inclusion_mstate& rhs) const {
            const auto *rhs_incl_ms = dynamic_cast<const inclusion_mstate*>(&rhs);
            return (state_ == rhs_incl_ms->state_);
        }

        /// ordering of inclusion macrostate
        bool lt(const inclusion_mstate& rhs) const {
            const auto *rhs_incl_ms = dynamic_cast<const inclusion_mstate*>(&rhs);
            if(state_ != rhs_incl_ms->state_) {return state_ < rhs_incl_ms->state_;}

            return false;
        }

        ~inclusion_mstate() {}

        const intersect_mstate& get_intersect_state() const { return state_; }

        friend class inclusion_check;
    };

    /// main class for on the fly inclussion procedure
    // class inclusion_check : public abstract_successor {
    class inclusion_check {
        /// target of a transition (including colours tagged by partition index)
        using state_col = std::pair<intersect_mstate, std::set<unsigned>>;
        /// vec_state_col where colours are tagged by their partition
        using vec_state_col = std::vector<state_col>;
    private:
        /// can be omitted
        spot::twa_graph_ptr preprocessed_orig_aut_B_;
        /// is not used, should be deleted 
        std::map<intersect_mstate , vec_state_col> intersect_states_;
        /// to not generate redundantly new instances
        std::map<unsigned, std::vector<std::pair<cola::tnba_complement::vec_state_taggedcol, bdd>>> compl_state_storage_;

        std::vector<std::shared_ptr<kofola::inclusion_mstate>> init_states_;
        spot::twa_graph_ptr aut_A_;
        bdd msupport_;
        bdd n_s_compat_;
        cola::tnba_complement aut_B_compl_;
        /// acc_cond that should be satisfied for inclusion to not hold
        spot::acc_cond::acc_code acc_cond_;
        /// acc mark for aut_A
        unsigned int first_col_to_use_;
        std::set<unsigned> infs_from_compl_;

        /// maps state of aut_A to states of aut_B that direct simulate it
        std::unordered_map<unsigned, std::vector<unsigned>> dir_simul_;
        /// offset computed for the use of computing direct simulation
        unsigned offset_ = 0;

    public:
        /// constructor that stores automata such that aut_A subset aut_B can be decided (calling inclusion() method should follow)
        inclusion_check(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B);

        /// decide which states from aut_B simulate states in aut_A
        void compute_simulation(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B);

        /// for debugging purposes only
        void print_mstate(const std::shared_ptr<inclusion_mstate> a);

        /// preprocess automaton A if required 
        spot::twa_graph_ptr init_aut_A(const spot::twa_graph_ptr &aut_A);

        /// to obtain cola::tnba_complement instance in the constructor and preprocess if required
        cola::tnba_complement init_compl_aut_b(const spot::twa_graph_ptr &aut_B);

        /// returns union of aut_A and aut_B for the purpose of simulations, firstly states from aut_A are inserted, then states 
        /// frin aut_B and finally, initial state of union automaton with transitions to initial state of aut_A and aut_B respectively
        /// TODO should be done properly (so far suffices)
        spot::twa_graph_ptr aut_union(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B);

        /// method to be called after instantiation of this class to decide inclusion, calls emptiness checker 
        bool inclusion();

        /// to preprocess autB
        /// TODO just copypaste from complement_sync.cpp
        spot::twa_graph_ptr preprocess(const spot::twa_graph_ptr &aut);

        /// returns "alphabet" of automaton aut_A
        std::pair<bdd, bdd> symbols_from_A(const spot::twa_graph_ptr &aut_A);

        /// implements getter fot initial states, so the emptiness check can obtain them
        std::vector<std::shared_ptr<inclusion_mstate>> get_initial_states();

        /// get successors for complement from compl_state over letter
        cola::tnba_complement::vec_state_taggedcol get_successors_compl(unsigned compl_state, const bdd& letter);

        /// returns set of all successors (inclusion macrostates) for given inclusion macrostate, for the need of emptiness check
        std::vector<std::shared_ptr<inclusion_mstate>> get_succs(const std::shared_ptr<inclusion_mstate> &src);

        /// returns true if a is early simul. less than b
        bool subsum_less_early(const std::shared_ptr<inclusion_mstate>& a, const std::shared_ptr<inclusion_mstate>& b);

        /// returns true if a is early+1 simul. less than b
        bool subsum_less_early_plus(const std::shared_ptr<inclusion_mstate>& a, const std::shared_ptr<inclusion_mstate>& b);

        /// returns true when the provided mark_t satisfies acc_cond_, for the need of emptiness check
        bool is_accepting(spot::acc_cond::mark_t inf_cond);

        /// decides whether the transition from src to dst over symbol is accepting in aut_A
        bool is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol);

        /// returns product: states_A x states_B
        std::vector<std::shared_ptr<inclusion_mstate>>
        get_cartesian_prod(unsigned aut_A_src, std::set<unsigned> &states_A,
                           cola::tnba_complement::vec_state_taggedcol &states_B,
                           const bdd &letter);

    };

    bool operator<(const inclusion_mstate& lhs,
                   const inclusion_mstate& rhs);
}

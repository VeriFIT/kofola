/**
 * @file hyperltl_mc.hpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief Declarations for hyperltl model checking
 * @version 0.1
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include "kofola.hpp"
#include "hyperltl_formula_processor.hpp"
#include "abstract_successor.hpp"

namespace kofola {
    class  hyperltl_mc_mstate;

    // main class implementing hyperltl model checking
    class hyperltl_mc : public abstract_successor {
        using mc_macrostate = std::vector<unsigned>; // first-system, following-formula

        const parsed_hyperltl_form_ptr& parsed_hyperltl_f_; /// info regarding formula
        spot::twa_graph_ptr built_aut_; /// inductively built automaton
        std::vector<spot::kripke_graph_ptr> kripke_structs_; /// kripke structs representing behaviors of systems
        std::vector<spot::kripke_graph_ptr> curr_kripke_structs_; /// kripke structure that is existentially quantified over
        bool one_system_only_; /// predicate to show only one system is used for all quantifications
    public:
        /// prints SAT or UNSAT based on whether parsed_hyperltl_f is satisfieed by kripke_structs 
        /// (maybe should be moved to separate method)
        hyperltl_mc(const parsed_hyperltl_form_ptr& parsed_hyperltl_f, std::vector<spot::kripke_graph_ptr> kripke_structs);

        /// moves n kripke structs from kripke_structs_ to curr_kripke_structs_
        void use_last_n_kripke_structs(unsigned n);

        /// returns vector of conditions corresponding to vector of states in sys_state, 
        /// where ith state corresponds to ith kripke struct in curr_kripke_structs_
        std::vector<bdd> get_sys_state_conds(std::vector<unsigned> sys_state);

        /// compute successors of states for each system in curr_kripke_structs_
        std::vector<std::vector<unsigned>> system_successors(std::vector<unsigned> src);

        /// removes variable vars from condition bdd (used for existential projection)
        bdd remove_bdd_vars(bdd cond, std::vector<int> vars);

        /// for each current system returns vector of initial states
        std::vector<std::vector<unsigned>> get_systems_init();

        /// performs existential projection with respect to trace variables in exist_trac_vars and builds corresponding 
        /// automaton
        spot::twa_graph_ptr existential_projection(const std::vector<std::string>& exist_trac_vars);

        /// computes product of N vectors
        std::vector<std::vector<unsigned>> prod(const std::vector<std::vector<unsigned>>& sets);

        /// returns bdd representing cond of system state only consisting of APs that are about to be removed and replacing 
        /// the bdd with system variables to bdd with automaton variables 
        bdd get_bdd_pair_system_to_aut(const std::string& trace_vars, spot::twa_graph_ptr &composed, bdd cond,  unsigned ith_sys);

        /// performs nfold self composition regarding trac_vars
        spot::twa_graph_ptr n_fold_self_composition(std::vector<std::string> trac_vars);

        /// for the use of emptiness check, returns if cond satisfies the condition of built_aut_
        bool is_accepting(spot::acc_cond::mark_t cond) override;

        /// returns initial states, for emptiness check 
        std::vector<std::shared_ptr<abstract_successor::mstate>> get_initial_states() override;

        /// interface for emptiness check to obtain successors of src
        std::vector<std::shared_ptr<abstract_successor::mstate>> get_succs(const std::shared_ptr<abstract_successor::mstate> &src) override;

        /// internal computation of successors for src over exist_trac_vars
        std::vector<std::shared_ptr<hyperltl_mc_mstate>> get_succs_internal(std::vector<unsigned> src, std::vector<std::string> exist_trac_vars);

        /// for debugging
        void print_mstate(const std::shared_ptr<abstract_successor::mstate> a) override {
            return;
        };
    };

/// class used as a macrostate for emptiness check
class  hyperltl_mc_mstate : public abstract_successor::mstate {
    private:
        std::vector<unsigned> state_; /// systems, last state of aut
    public:
        hyperltl_mc_mstate() {
            ;
        }

        bool eq(const abstract_successor::mstate& rhs) const override {
            const auto *rhs_hypetltlmc_ms = dynamic_cast<const hyperltl_mc_mstate*>(&rhs);
            return (state_ == rhs_hypetltlmc_ms->state_);
        }

        bool lt(const abstract_successor::mstate& rhs) const override {
            const auto *rhs_incl_ms = dynamic_cast<const hyperltl_mc_mstate*>(&rhs);
            if(state_ != rhs_incl_ms->state_) {return state_ < rhs_incl_ms->state_;}
            return false;
        }

        friend class hyperltl_mc;
    };

} // kofola

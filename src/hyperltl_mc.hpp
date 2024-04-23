#pragma once

#include "kofola.hpp"
#include "hyperltl_formula_processor.hpp"
#include "abstract_successor.hpp"

namespace kofola {
    class  hyperltl_mc_mstate;

    class hyperltl_mc : public abstract_successor {
        using mc_macrostate = std::vector<unsigned>; // first-system, second-formula

        const parsed_hyperltl_form_ptr& parsed_hyperltl_f_; /// info regarding formula
        spot::twa_graph_ptr built_aut_; /// inductively built automaton
        std::vector<spot::kripke_graph_ptr> kripke_structs_; /// kripke struct representing system behavior
        std::vector<spot::kripke_graph_ptr> curr_kripke_structs_;
        bool one_system_only_;
    public:
        hyperltl_mc(const parsed_hyperltl_form_ptr& parsed_hyperltl_f, std::vector<spot::kripke_graph_ptr> kripke_structs);

        void use_last_n_kripke_structs(unsigned n);

        std::vector<bdd> get_sys_state_conds(std::vector<unsigned> sys_state);

        std::vector<std::vector<unsigned>> system_successors(std::vector<unsigned> src);

        bdd remove_bdd_vars(bdd cond, std::vector<int> vars);

        std::vector<std::vector<unsigned>> get_systems_init();

        spot::twa_graph_ptr existential_projection(const std::vector<std::string>& exist_trac_vars);

        std::vector<std::vector<unsigned>> prod(const std::vector<std::vector<unsigned>>& sets);

        bdd get_bdd_pair_system_to_aut(const std::string& trace_vars, spot::twa_graph_ptr &composed, bdd cond,  unsigned ith_sys);

        spot::twa_graph_ptr n_fold_self_composition(std::vector<std::string> trac_vars);

        bool is_accepting(spot::acc_cond::mark_t cond) override;

        std::vector<std::shared_ptr<abstract_successor::mstate>> get_initial_states() override;

        std::vector<std::shared_ptr<abstract_successor::mstate>> get_succs(const std::shared_ptr<abstract_successor::mstate> &src) override;

        std::vector<std::shared_ptr<hyperltl_mc_mstate>> get_succs_internal(std::vector<unsigned> src, std::vector<std::string> exist_trac_vars);

        void print_mstate(const std::shared_ptr<abstract_successor::mstate> a) override {
            return;
        };
    };

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

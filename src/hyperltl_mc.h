#pragma once

#include "kofola.hpp"
#include "hyperltl_formula_processor.hpp"

namespace kofola {
    class hyperltl_mc {
        using mc_macrostate = std::pair<unsigned, unsigned>; // first-system, second-formula

        const parsed_hyperltl_form_ptr& parsed_hyperltl_f_; /// info regarding formula
        spot::twa_graph_ptr built_aut_; /// inductively built automaton
        spot::kripke_graph_ptr system_; /// kripke struct representing system behavior
        bddPair *aut_to_system_;
        bdd aps_to_remove_;
        bool negate_; /// indicates if formula preprocessing caused the top level negation
        bdd aps_not_in_system_;
    public:
        hyperltl_mc(const parsed_hyperltl_form_ptr& parsed_hyperltl_f, spot::kripke_graph_ptr  system);

        void mark_redundant_aps_formula();

        std::vector<unsigned> system_successors(spot::kripke_graph_state* s);

        bddPair *get_bdd_pair_aut_to_system();

        bdd get_relevant_aut_aps(const std::vector<std::string>& exist_trac_vars, spot::twa_graph_ptr &projected);

        spot::twa_graph_ptr existential_projection(const std::vector<std::string>& exist_trac_vars);

        std::vector<std::vector<unsigned>> prod(const std::vector<std::vector<unsigned>>& sets);

        bdd get_bdd_pair_system_to_aut(const std::string& trace_vars, spot::twa_graph_ptr &composed, bdd cond);

        spot::twa_graph_ptr n_fold_self_composition(std::vector<std::string> trac_vars);
    };

} // kofola

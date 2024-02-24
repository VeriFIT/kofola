#include "hyperltl_mc.h"
#include "inclusion_test.hpp"
#include <spot/twa/twagraph.hh>
#include <utility>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/genem.hh>

#include <map>

namespace kofola {
    hyperltl_mc::hyperltl_mc(const parsed_hyperltl_form_ptr& parsed_hyperltl_f, spot::kripke_graph_ptr  system):
    parsed_hyperltl_f_(parsed_hyperltl_f),
    built_aut_(parsed_hyperltl_f->aut),
    system_(std::move(system))
    {
        Quantification q;
        while(!(parsed_hyperltl_f->q_list.empty())) {
            q = parsed_hyperltl_f->q_list.front();
            parsed_hyperltl_f->q_list.pop_front();

            if(q.type == static_cast<unsigned int>(QuantificationType::Exists)) {
                existential_projection(q.trace_vars);
            }
            else if(parsed_hyperltl_f->q_list.empty()) {
                // n-fold self composition and inclusion,
                // although when there is only seq. of existential, we can
                // perform emptiness check on the fly
                unsigned n = q.trace_vars.size();
                auto aut_A = n_fold_self_composition(n, q.trace_vars);
                kofola::inclusionTest().test(aut_A, built_aut_)
                ? std::cout << "SAT" << std::endl
                : std::cout << "UNSAT" << std::endl;

                return;
            }
            else {
                built_aut_ = kofola::complement_tela(built_aut_);
                built_aut_ = existential_projection(q.trace_vars);
                built_aut_ = kofola::complement_tela(built_aut_);
            }
        }
    }

    std::list<unsigned> hyperltl_mc::system_successors(spot::kripke_graph_state* s) {
        std::list<unsigned> res;

        for (auto i: system_->succ(s))
        {
            auto tmp = i->dst();
            res.emplace_back(system_->state_number(tmp));
        }

        return res;
    }

    bddPair *hyperltl_mc::get_bdd_pair_aut_to_system() {
        auto dict = built_aut_->get_dict();
        auto system_dict = system_->get_dict();

        std::vector<int> aut_part;
        std::vector<int> sys_part;

        for(const auto& ap: built_aut_->ap()) {
            for (const auto &system_ap: system_->ap()) {
                if (system_ap.ap_name() == parsed_hyperltl_f_->aps_map[ap.ap_name()].atomic_prop) {
                    aut_part.emplace_back(dict->var_map[ap]);
                    sys_part.emplace_back(system_dict->var_map[system_ap]);
                    break;
                }
            }
        }

        bddPair *res = bdd_newpair();
        bdd_setpairs(res, aut_part.data(), sys_part.data(), aut_part.size());
        return res;
    }

    bdd hyperltl_mc::get_relevant_aut_aps(const std::list<std::string>& exist_trac_vars) {
        auto dict = built_aut_->get_dict();
        auto system_dict = system_->get_dict();
        // obtain all aps from the automaton
        std::vector<int> aps_to_remove;
        for(const auto& ap: built_aut_->ap()) {
            bool remove = true;

            for(const auto& trace_var: exist_trac_vars) {
                if(parsed_hyperltl_f_->aps_map[ap.ap_name()].trace_var == trace_var) {
                    remove = false;
                    break;
                }
            }

            if(remove)
                aps_to_remove.emplace_back(dict->var_map[ap]);
        }

        return bdd_makeset(aps_to_remove.data(), aps_to_remove.size());
    }

    spot::twa_graph_ptr hyperltl_mc::existential_projection(const std::list<std::string>& exist_trac_vars) {
        auto init_formula = parsed_hyperltl_f_->aut->get_init_state_number();
        auto init_system = system_->get_init_state_number();

        mc_macrostate initial_macrostate = std::make_pair(init_formula, init_system);

        auto projected = spot::make_twa_graph(built_aut_->get_dict());
        unsigned spot_state = projected->new_state();

        std::queue<std::pair<mc_macrostate, unsigned>> macrostates_spotstate;
        macrostates_spotstate.emplace(std::make_pair(initial_macrostate, spot_state));

        std::pair<mc_macrostate, unsigned> s;
        std::map<mc_macrostate, bool> used_states;
        used_states[initial_macrostate] = true;

        while(!macrostates_spotstate.empty()) {
            s = macrostates_spotstate.back();
            macrostates_spotstate.pop();

            unsigned system_src = s.first.first;
            unsigned aut_src = s.first.second;
            unsigned new_aut_src = s.second;

            auto kripke_state = system_->state_from_number(system_src);
            auto system_cond = system_->state_condition(kripke_state);

            // TODO consider only relevant cond
            auto system_succs = system_successors(kripke_state);

            auto aut_aps = built_aut_->ap();


            auto aps_to_remove = get_relevant_aut_aps(exist_trac_vars); // aps to remove

            for (const auto &t : built_aut_->out(aut_src)) {
                auto relevant_aut_aps = bdd_exist(t.cond, aps_to_remove);
                auto to_compare = bdd_replace(relevant_aut_aps,get_bdd_pair_aut_to_system());

                if (bdd_implies(system_cond, to_compare)) {
                    for(auto s_succ: system_succs) {
                        mc_macrostate mstate = std::make_pair(s_succ, t.dst);
                        spot_state = projected->new_state();
                        projected->new_edge(new_aut_src, spot_state, system_cond, t.acc);
                        if(used_states.count(mstate) == 0) {
                            macrostates_spotstate.emplace(std::make_pair(mstate, spot_state));
                            used_states[mstate] = true;
                        }
                    }
                }
            }
        }

        return projected;
    }

    spot::twa_graph_ptr hyperltl_mc::n_fold_self_composition(unsigned n, const std::list<std::string>& exist_trac_vars) {

    }
} // kofola
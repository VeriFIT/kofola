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
        mark_redundant_aps_formula();
        Quantification q;
        while(!(parsed_hyperltl_f->q_list.empty())) {
            q = parsed_hyperltl_f->q_list.back();
            parsed_hyperltl_f->q_list.pop_back();

            if(q.type == static_cast<unsigned int>(QuantificationType::Exists) && !parsed_hyperltl_f->q_list.empty()) {
                for(auto trace_var: q.trace_vars) {
                    built_aut_ = existential_projection({trace_var});
                }
            }
            else if(q.type == static_cast<unsigned int>(QuantificationType::Exists) && parsed_hyperltl_f->q_list.empty()) {
                for(int i = q.trace_vars.size() - 1; i >= 0; i--) {
                    built_aut_ = existential_projection({q.trace_vars[i]});
                }
                spot::generic_emptiness_check(built_aut_)
                ? std::cout << "UNSAT" << std::endl
                : std::cout << "SAT" << std::endl;

                return;
            }
            else if(q.type == static_cast<unsigned int>(QuantificationType::Forall) && parsed_hyperltl_f->q_list.empty()) {
                // n-fold self composition and inclusion,
                // although when there is only seq. of existential, we can
                // perform emptiness check on the fly
                auto aut_A = n_fold_self_composition(q.trace_vars);
                kofola::inclusionTest inclusion;
                inclusion.test(aut_A, built_aut_)
                ? std::cout << "SAT" << std::endl
                : std::cout << "UNSAT" << std::endl;

                return;
            }
            else {
                // TODO when innermost, the formula can be negated only ??
                built_aut_ = kofola::complement_tela(built_aut_);
                built_aut_ = existential_projection(q.trace_vars);
                built_aut_ = kofola::complement_tela(built_aut_);
            }
        }
    }

    void hyperltl_mc::mark_redundant_aps_formula() {
        std::vector<int> to_remove;
        for(const auto& keys: parsed_hyperltl_f_->aps_map) {
            bool remove = true;
            for(const auto& sys_ap: system_->ap()) {
                if(sys_ap.ap_name() == keys.second.atomic_prop) {
                    remove = false;
                    break;
                }
            }
            if(remove) {
                auto ap = built_aut_->register_ap(keys.first);
                to_remove.emplace_back(ap);
            }
        }
        aps_not_in_system_ = bdd_makeset(to_remove.data(), to_remove.size());
    }

    std::vector<unsigned> hyperltl_mc::system_successors(spot::kripke_graph_state* s) {
        std::vector<unsigned> res;

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
        aut_to_system_ = res;
        return res;
    }

    bdd hyperltl_mc::get_relevant_aut_aps(const std::vector<std::string>& exist_trac_vars, spot::twa_graph_ptr &projected) {
        auto dict = built_aut_->get_dict();
        auto system_dict = system_->get_dict();

        // obtain all aps from the automaton
        std::vector<int> aps_to_remove;
        std::vector<int> aps_to_keep;
        for(const auto& ap: built_aut_->ap()) {
            bool keep = true;

            for(const auto& trace_var: exist_trac_vars) {
                if(parsed_hyperltl_f_->aps_map[ap.ap_name()].trace_var == trace_var) {
                    keep = false;
                    break;
                }
            }

            if(keep)
            {aps_to_keep.emplace_back(dict->var_map[ap]); projected->register_ap(ap);}
            else
            {aps_to_remove.emplace_back(dict->var_map[ap]);}
        }

        if(aps_to_keep.empty())
            return bdd_false();

        aps_to_remove_ = bdd_makeset(aps_to_remove.data(), aps_to_remove.size());
        return bdd_makeset(aps_to_keep.data(), aps_to_keep.size());
    }

    spot::twa_graph_ptr hyperltl_mc::existential_projection(const std::vector<std::string>& exist_trac_vars) {
        auto init_formula = built_aut_->get_init_state_number();
        auto init_system = system_->get_init_state_number();
        mc_macrostate initial_macrostate = std::make_pair(init_system, init_formula);

        auto projected = spot::make_twa_graph(built_aut_->get_dict());
        unsigned spot_state = projected->new_state();

        std::queue<std::pair<mc_macrostate, unsigned>> macrostates_spotstate;
        macrostates_spotstate.emplace(std::make_pair(initial_macrostate, spot_state));

        std::pair<mc_macrostate, unsigned> s;
        std::map<mc_macrostate, unsigned> used_states;
        used_states[initial_macrostate] = spot_state;

        while(!macrostates_spotstate.empty()) {
            s = macrostates_spotstate.front();
            macrostates_spotstate.pop();

            unsigned system_src = s.first.first;
            unsigned aut_src = s.first.second;
            unsigned new_aut_src = s.second;

            auto kripke_state = system_->state_from_number(system_src);
            auto system_cond = system_->state_condition(kripke_state);
            auto system_succs = system_successors(kripke_state);

            auto not_interesting_ap = get_relevant_aut_aps(exist_trac_vars, projected); // aps to remove
            for (const auto &t : built_aut_->out(aut_src)) {
                auto relevant_aut_aps = bdd_exist(bdd_exist(t.cond, aps_not_in_system_), not_interesting_ap);
                auto to_compare = bdd_replace(relevant_aut_aps,get_bdd_pair_aut_to_system());
                if (bdd_implies(system_cond, to_compare)) {
                    for(auto s_succ: system_succs) {
                        mc_macrostate mstate = std::make_pair(s_succ, t.dst);
                        if(used_states.count(mstate) == 0) {
                            spot_state = projected->new_state();
                            macrostates_spotstate.push(std::make_pair(mstate, spot_state));
                            used_states[mstate] = spot_state;
                        }
                        else {
                            spot_state = used_states[mstate];
                        }
                        projected->new_edge(new_aut_src, spot_state, bdd_exist(bdd_exist(t.cond, aps_not_in_system_), aps_to_remove_), t.acc);
                    }
                }
            }
        }
        projected->remove_unused_ap();
        return projected;
    }

    bdd hyperltl_mc::get_bdd_pair_system_to_aut(const std::string& trace_var, spot::twa_graph_ptr &composed, bdd cond) {
        auto system_dict = system_->get_dict();

        std::vector<int> new_aps;
        std::vector<int> old_aps;

        for(const auto &system_ap: system_->ap()) {
            for(const auto &composed_ap: composed->ap()) {
                if(parsed_hyperltl_f_->aps_map[composed_ap.ap_name()].atomic_prop ==  system_ap.ap_name()) {
                    new_aps.emplace_back(composed->get_dict()->var_map[composed_ap]);
                    old_aps.emplace_back(system_dict->var_map[system_ap]);
                    break;
                }
            }
        }

        auto pair = bdd_newpair();
        bdd_setpairs(pair, old_aps.data(), new_aps.data(), old_aps.size());
        auto res = bdd_replace(cond, pair);
        return res;
    }

    std::vector<std::vector<unsigned>> hyperltl_mc::prod(const std::vector<std::vector<unsigned>>& sets) {
        std::vector<std::vector<unsigned>> res = {{}};
        for(const auto & set : sets) {
            std::vector<std::vector<unsigned>> new_res;
            for(const auto& product: res) {
                for(auto val: set) {
                    auto tmp = product;
                    tmp.emplace_back(val);
                    new_res.emplace_back(tmp);
                }
            }
            res = new_res;
        }

        return res;
    }

    spot::twa_graph_ptr hyperltl_mc::n_fold_self_composition(std::vector<std::string> trac_vars) {
        auto dict = spot::make_bdd_dict();
        auto self_composition = spot::make_twa_graph(dict);
        for(const auto &system_ap: system_->ap()) {
            for(const auto& trace_var: trac_vars) {
                self_composition->register_ap("{" + system_ap.ap_name() + "}_{" + trace_var + "}");
            }
        }
        auto cond = spot::acc_cond::acc_code::inf({0});
        self_composition->set_acceptance(cond);

        auto n = trac_vars.size();
        auto init_system = system_->get_init_state_number();

        std::vector<unsigned> init_self_comp_state(n, init_system);
        unsigned spot_state = self_composition->new_state();
        std::queue<std::pair<std::vector<unsigned>, unsigned>> macrostates_spotstate;
        macrostates_spotstate.emplace(std::make_pair(init_self_comp_state, spot_state));

        std::map<std::vector<unsigned>, unsigned> used_states;
        used_states[init_self_comp_state] = spot_state;
        std::pair<std::vector<unsigned>, unsigned> s;
        while(!macrostates_spotstate.empty()) {
            s = macrostates_spotstate.front();
            macrostates_spotstate.pop();

            auto src_state_vec = s.first;
            unsigned new_aut_src = s.second;

            std::vector<std::vector<unsigned>> to_cross_prod;
            std::vector<std::vector<unsigned>> succs;
            bdd trans = bddtrue;

            size_t i = 0;
            for(auto src_state: src_state_vec) {
                auto kripke_state = system_->state_from_number(src_state);
                // translate system AP to automaton AP
                auto transl_aps = get_bdd_pair_system_to_aut(trac_vars[i], self_composition, kripke_state->cond()); // conjunction of aps of all states;
                i++;
                trans = bdd_and(trans, transl_aps);
                to_cross_prod.emplace_back(system_successors(kripke_state));
            }
            auto all_succs = prod(to_cross_prod);
            for(const auto& succ: all_succs) {
                if(used_states.count(succ) == 0) {
                    spot_state = self_composition->new_state();
                    macrostates_spotstate.push(std::make_pair(succ, spot_state));
                    used_states[succ] = spot_state;
                }
                else {
                    spot_state = used_states[succ];
                }

                std::vector<unsigned> acceptance(1, 0);
                spot::acc_cond::mark_t spot_cols(acceptance.begin(), acceptance.end());
                self_composition->new_edge(new_aut_src, spot_state, trans, spot_cols);
            }

// TODO used state na zaklade spot_state => moze sa vyuzit unordered_map
        }
        return self_composition;
    }
} // kofola
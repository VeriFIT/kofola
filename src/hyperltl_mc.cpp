#include "hyperltl_mc.hpp"
#include "inclusion_check.hpp"
#include "emptiness_check.hpp"
#include <spot/twa/twagraph.hh>
#include <utility>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/genem.hh>

#include <map>

namespace kofola {
    hyperltl_mc::hyperltl_mc(const parsed_hyperltl_form_ptr& parsed_hyperltl_f, std::vector<spot::kripke_graph_ptr> kripke_structs):
    parsed_hyperltl_f_(parsed_hyperltl_f),
    built_aut_(parsed_hyperltl_f->aut),
    kripke_structs_(std::move(kripke_structs))
    {
        //mark_redundant_aps_formula();
        Quantification q;

        one_system_only_ = (kripke_structs_.size() == 1);
        if(one_system_only_) {
            // make n copies of system
            while(kripke_structs_.size() != parsed_hyperltl_f_->qantifiers) {
                kripke_structs_.emplace_back(kripke_structs_.back());
            }
        }

        bool sat = false;

        while(!(parsed_hyperltl_f->q_list.empty())) {
            q = parsed_hyperltl_f->q_list.back();
            parsed_hyperltl_f->q_list.pop_back();

            if(q.type == static_cast<unsigned int>(QuantificationType::Exists) && !parsed_hyperltl_f->q_list.empty()) {
                use_last_n_kripke_structs(q.trace_vars.size());
                built_aut_ = existential_projection(q.trace_vars);
            }
            else if(q.type == static_cast<unsigned int>(QuantificationType::Exists) && parsed_hyperltl_f->q_list.empty()) {
                use_last_n_kripke_structs(q.trace_vars.size());

                emptiness_check emptiness_checker(this);

                if(emptiness_checker.empty())
                    sat = false;
                else
                    sat = true;

                break;
            }
            else if(q.type == static_cast<unsigned int>(QuantificationType::Forall) && parsed_hyperltl_f->q_list.empty()) {
                // n-fold self composition and inclusion,
                // although when there is only seq. of existential, we can
                // perform emptiness check on the fly
                use_last_n_kripke_structs(q.trace_vars.size());

                auto aut_A = n_fold_self_composition(q.trace_vars);
                kofola::inclusion_check inclusion_checker(aut_A, built_aut_);
                // spot::print_hoa(std::cout, aut_A); spot::print_hoa(std::cout, built_aut_);

                if(inclusion_checker.inclusion())
                    sat = true;
                else
                    sat = false;

                break;
            }
            else {
                // TODO when innermost, the formula can be negated only ??
                //built_aut_ = ;
                kofola::OPTIONS.output_type = "buchi";
                built_aut_ = kofola::complement_tela(built_aut_);
                built_aut_ = existential_projection(q.trace_vars);
                built_aut_ = kofola::complement_tela(built_aut_);
            }
        }

        if((!sat && !parsed_hyperltl_f_->negate) || (sat && parsed_hyperltl_f_->negate))
            std::cout << "UNSAT\n";
        else if((sat && !parsed_hyperltl_f_->negate) || (!sat && parsed_hyperltl_f_->negate))
            std::cout << "SAT\n";
    }

    void hyperltl_mc::use_last_n_kripke_structs(unsigned n) {
        curr_kripke_structs_.clear();
        for(unsigned i = 0; i < n; i++) {
            curr_kripke_structs_.insert(curr_kripke_structs_.begin(), kripke_structs_.back());
            kripke_structs_.pop_back();
        }
    }

    std::vector<bdd> hyperltl_mc::get_sys_state_conds(std::vector<unsigned> sys_state) {
        std::vector<bdd> res;

        for(unsigned i = 0; i < sys_state.size(); i++) {
            auto ks_num = sys_state[i];
            auto ks_s = curr_kripke_structs_[i]->state_from_number(ks_num);
            auto cond = curr_kripke_structs_[i]->state_condition(ks_s);

            res.emplace_back(cond);
        }

        return res;
    }

    std::vector<std::vector<unsigned>> hyperltl_mc::system_successors(std::vector<unsigned> src) {
        std::vector<std::vector<unsigned>> res;

        // obtain kripke states reprs
        for(unsigned i = 0; i < src.size(); i++) {
            auto ks_num = src[i];
            auto ks_s = curr_kripke_structs_[i]->state_from_number(ks_num);

            std::vector<unsigned> partial_res;
            for (auto t: curr_kripke_structs_[i]->succ(ks_s))
            {
                auto tmp = t->dst();
                partial_res.emplace_back(curr_kripke_structs_[i]->state_number(tmp));
            }
            res.emplace_back(partial_res);
        }

        return res;
    }

    bdd hyperltl_mc::remove_bdd_vars(bdd cond, std::vector<int> vars) {
        bdd res = bdd_makeset(vars.data(), vars.size());
        res = bdd_exist(cond, res);

        return res;
    }

    bool hyperltl_mc::is_accepting(spot::acc_cond::mark_t cond) {
        return built_aut_->get_acceptance().accepting(cond);
    }

    std::vector<std::unique_ptr<hyperltl_mc_mstate>> hyperltl_mc::get_succs_internal(std::vector<unsigned> src, std::vector<std::string> exist_trac_vars) {
        std::vector<std::unique_ptr<hyperltl_mc_mstate>> res;

        std::vector<unsigned> system_src_states(src.begin(), src.end() - 1);
        std::vector<std::vector<unsigned>> sets_of_sys_succs = system_successors(system_src_states);
        auto state_conds = get_sys_state_conds(system_src_states);

        auto aut_src = src.back();
        for (const auto &t : built_aut_->out(aut_src)) {
            bdd restricted = t.cond;
            for(int i = curr_kripke_structs_.size() - 1; i >= 0; i--) {
                auto to_restrict = get_bdd_pair_system_to_aut(exist_trac_vars[i], built_aut_, state_conds[i], i);
                restricted = bdd_restrict(restricted, to_restrict);
            }

            if (restricted != bddfalse) {
                auto to_prod = sets_of_sys_succs;
                std::vector<unsigned> aut_dst = {t.dst};
                to_prod.emplace_back(aut_dst);
                auto product = prod(to_prod);
                for(const auto& mstate: product) {
                    auto ptr = std::make_unique<hyperltl_mc_mstate>();
                    ptr->state_ = mstate; ptr->acc_ = t.acc; ptr->trans_cond_ = restricted;
                    res.emplace_back(std::move(ptr));
                }
            }
        }

        return res;
    }

    std::vector<std::unique_ptr<abstr_succ::abstract_successor::mstate>> hyperltl_mc::get_initial_states() {
        std::vector<std::unique_ptr<abstract_successor::mstate>> res;

        std::vector<unsigned> init_aut = {built_aut_->get_init_state_number()};
        auto sys_init = get_systems_init();

        std::vector<std::vector<unsigned>> to_prod = sys_init;
        to_prod.emplace_back(init_aut);
        auto initial_macrostates = prod(to_prod);
        for(auto init: initial_macrostates) {
            auto ptr = std::make_unique<hyperltl_mc_mstate>();
            ptr->state_ = init;
            res.emplace_back(std::move(ptr));
        }

        return res;
    }


    std::vector<std::unique_ptr<abstr_succ::abstract_successor::mstate>> hyperltl_mc::get_succs(const std::unique_ptr<abstract_successor::mstate> &src) {
        auto casted_src = dynamic_cast<hyperltl_mc_mstate*>(src.get());
        auto succs = get_succs_internal(casted_src->state_, this->parsed_hyperltl_f_->q_list.front().trace_vars);
        std::vector<std::unique_ptr<abstract_successor::mstate>> res;

        for (auto& succ : succs) {
            res.emplace_back(std::move(succ));
        }

        return res;
    }

    std::vector<std::vector<unsigned>> hyperltl_mc::get_systems_init() {
        std::vector<std::vector<unsigned>> res;
        for(auto & curr_kripke_struct : curr_kripke_structs_) { // TODO podmienka
            auto init = curr_kripke_struct->get_init_state_number();
            auto kripke_state = curr_kripke_struct->state_from_number(init);
            std::vector<unsigned> initial_sys_states = {init};
            if(kripke_state->cond() == bddtrue) { // this should be done properly, just hack
                initial_sys_states = system_successors({init}).front();
            }
            res.emplace_back(initial_sys_states);
        }

        return res;
    }

    spot::twa_graph_ptr hyperltl_mc::existential_projection(const std::vector<std::string>& exist_trac_vars) {
        std::vector<unsigned> init_aut = {built_aut_->get_init_state_number()};
        auto sys_init = get_systems_init();

        auto projected = spot::make_twa_graph(built_aut_->get_dict());
        for(const auto& orig_ap: built_aut_->ap()) {
            projected->register_ap(orig_ap);
        }
        projected->set_acceptance(built_aut_->get_acceptance());

        std::queue<std::pair<mc_macrostate, unsigned>> macrostates_spotstate;
        std::vector<std::vector<unsigned>> to_prod = sys_init;
        to_prod.emplace_back(init_aut);
        auto initial_macrostates = prod(to_prod);

        unsigned spot_state = projected->new_state();
        projected->set_init_state(spot_state);

        std::map<mc_macrostate, unsigned> used_states;
        for(const auto& init_mstate: initial_macrostates) {
            macrostates_spotstate.emplace(std::make_pair(init_mstate, spot_state));
            used_states[init_mstate] = spot_state;
        }

        std::pair<mc_macrostate, unsigned> s;
        while(!macrostates_spotstate.empty()) {
            s = macrostates_spotstate.front();
            macrostates_spotstate.pop();

            auto src_mstate = s.first;
            std::vector<unsigned> system_srcs(src_mstate.begin(), src_mstate.end() - 1);
            unsigned new_aut_src = s.second;

            auto succs = get_succs_internal(src_mstate, exist_trac_vars);
            for(const auto& mstate: succs) {
                if(used_states.count(mstate->state_) == 0) {
                    spot_state = projected->new_state();
                    macrostates_spotstate.push(std::make_pair(mstate->state_, spot_state));
                    used_states[mstate->state_] = spot_state;
                }
                else {
                    spot_state = used_states[mstate->state_];
                }
                projected->new_edge(new_aut_src, spot_state, mstate->trans_cond_, mstate->acc_);
            }
        }
        projected->remove_unused_ap();
        return projected;
    }

    bdd hyperltl_mc::get_bdd_pair_system_to_aut(const std::string& trace_var, spot::twa_graph_ptr &composed, bdd cond, unsigned ith_sys) {
        auto system_dict = curr_kripke_structs_[ith_sys]->get_dict();

        std::vector<int> new_aps;
        std::vector<int> old_aps;

        std::vector<int> new_domain;

        bool has_eq;
        for(const auto &system_ap: curr_kripke_structs_[ith_sys]->ap()) {
            has_eq = false;
            for(const auto &composed_ap: composed->ap()) {
                if(parsed_hyperltl_f_->aps_map[composed_ap.ap_name()].atomic_prop ==  system_ap.ap_name() &&
                        parsed_hyperltl_f_->aps_map[composed_ap.ap_name()].trace_var == trace_var) {
                    new_aps.emplace_back(composed->get_dict()->var_map[composed_ap]);
                    old_aps.emplace_back(system_dict->var_map[system_ap]);

                    has_eq = true;
                    break;
                }
            }
            if(!has_eq)
                new_domain.emplace_back(system_dict->var_map[system_ap]);
        }

        auto pair = bdd_newpair();
        bdd_setpairs(pair, old_aps.data(), new_aps.data(), old_aps.size());
        auto res = bdd_replace(remove_bdd_vars(cond, new_domain), pair);

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
        auto self_composition = spot::make_twa_graph(built_aut_->get_dict());
        for(unsigned i = 0; i < curr_kripke_structs_.size(); i++) {
            for (const auto &system_ap: curr_kripke_structs_[i]->ap()) {
                self_composition->register_ap("{" + system_ap.ap_name() + "}_{" + trac_vars[i] + "}");
            }
        }

        auto cond = spot::acc_cond::acc_code::inf({0});
        self_composition->set_acceptance(cond);

        auto init_sys_states = get_systems_init();
        auto all_comb_init = prod(init_sys_states);

        unsigned spot_state = self_composition->new_state();
        std::queue<std::pair<std::vector<unsigned>, unsigned>> macrostates_spotstate;
        std::map<std::vector<unsigned>, unsigned> used_states;
        for(const auto& init_self_comp_state: all_comb_init) {
            macrostates_spotstate.emplace(std::make_pair(init_self_comp_state, spot_state));
            used_states[init_self_comp_state] = spot_state;
        }


        std::pair<std::vector<unsigned>, unsigned> s;
        while(!macrostates_spotstate.empty()) {
            s = macrostates_spotstate.front();
            macrostates_spotstate.pop();

            auto src_state_vec = s.first;
            unsigned new_aut_src = s.second;

            std::vector<std::vector<unsigned>> to_cross_prod;
            std::vector<std::vector<unsigned>> succs;
            bdd trans = bddtrue;

            for(unsigned i = 0; i < src_state_vec.size(); i++) {
                auto ks_num = src_state_vec[i];
                auto ks_s = curr_kripke_structs_[i]->state_from_number(ks_num);

                auto transl_aps = get_bdd_pair_system_to_aut(trac_vars[i], self_composition, ks_s->cond(), i); // conjunction of aps of all states;
                trans = bdd_and(trans, transl_aps);

                std::vector<unsigned> partial_res;
                for (auto t: curr_kripke_structs_[i]->succ(ks_s))
                {
                    auto tmp = t->dst();
                    partial_res.emplace_back(curr_kripke_structs_[i]->state_number(tmp));
                }
                to_cross_prod.emplace_back(partial_res);
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
        }
        return self_composition;
    }
} // kofola
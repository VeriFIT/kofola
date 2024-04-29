#include "emptiness_check.hpp"
#include <tuple>

// kofola
#include "complement_tela.hpp"
#include "util.hpp"
#include "decomposer.hpp"
#include "complement_sync.hpp"
#include "inclusion_check.hpp"
#include "hyperltl_mc.hpp"

// Spot
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complete.hh>

// standard library
#include <queue>
#include <utility>

namespace kofola {
    emptiness_check::emptiness_check(abstract_successor *as):
    abstr_succ_(as)
    {
    }

    bool emptiness_check::empty() {
        auto init_states = abstr_succ_->get_initial_states();
        for (const auto& state: init_states) {
            dfs_num_.insert({state, UNDEFINED});
            on_stack_.insert({state, false});
        }

        for(const auto& init: init_states) {
            if (dfs_num_.at(init) == UNDEFINED) {
                if(kofola::OPTIONS.params.count("use_tough_opt") != 0 && kofola::OPTIONS.params["use_tough_opt"] == "yes")
                    couvrer_edited(init, spot::acc_cond::mark_t());
                else
                    tarjan_is_empty(init, spot::acc_cond::mark_t());
                if(decided_) {
                    std::cout << cnt_ << "\n";
                    return empty_;
                }
            }
        }

        std::cout << cnt_ << "\n";
        return empty_;
    }

    bool emptiness_check::check_simul_less(const std::shared_ptr<abstract_successor::mstate> &dst_mstate) {
        auto cond = dst_mstate->get_acc();
        for (auto it = dfs_acc_stack_.rbegin(); it != dfs_acc_stack_.rend() && cond == 0; ++it) {
            const auto &s = (*it).first;
            if (abstr_succ_->subsum_less_early(dst_mstate, s)) {
                for (auto it2 = dfs_acc_stack_.rbegin(); it2 != dfs_acc_stack_.rend() && (*it2).first != s; ++it2) {
                    auto s_between = (*it2).first;
                    if(state_jumps_to_cutoffs_.count(s_between) == 0) state_jumps_to_cutoffs_.insert({s_between, {dst_mstate}}); // init
                    else state_jumps_to_cutoffs_[s_between].insert(dst_mstate);
                }
                return true;
            }

            cond |= (*it).second;
        }

        return false;
    }

    void emptiness_check::couvrer_edited(const std::shared_ptr<abstract_successor::mstate> &src_mstate, spot::acc_cond::mark_t path_cond) {
        /// STRONGCONNECT
//        if(dfs_num_[src_mstate] == UNDEFINED)
//            cnt_++;
        cnt_++;
        if(abstr_succ_->is_accepting(path_cond) && simulation_prunning(src_mstate))
            return;

        SCCs_.push(src_mstate);
        dfs_acc_stack_.emplace_back(src_mstate, src_mstate->get_acc());
        dfs_num_[src_mstate] = index_;
        index_++;
        tarjan_stack_.push_back(src_mstate);
        on_stack_[src_mstate] = true;

        auto succs = abstr_succ_->get_succs(src_mstate);
        for (auto &dst_mstate: succs) {
            // checking emptiness of curr dst_mstate
            if(empty_lang(dst_mstate))
                continue;

            // init structures
            if(dfs_num_.count(dst_mstate) == 0)
            {
                dfs_num_.insert({dst_mstate, UNDEFINED});
                on_stack_.insert({dst_mstate, false});
            }

            if (dfs_num_[dst_mstate] == UNDEFINED && !check_simul_less(dst_mstate))
            {
                couvrer_edited( dst_mstate, (path_cond | dst_mstate->get_acc()) );
                if(decided_)
                    return;
            }
            else if(dfs_num_[dst_mstate] != UNDEFINED) {
                if(on_stack_[dst_mstate] && merge_acc_marks(dst_mstate))
                    return;

                if(state_jumps_to_cutoffs_.count(dst_mstate) == 0)
                    continue;
                for(auto &jumping_dst_mstate: state_jumps_to_cutoffs_[dst_mstate]) {
                    if(dfs_num_[jumping_dst_mstate] == UNDEFINED && !check_simul_less(jumping_dst_mstate))
                        couvrer_edited( jumping_dst_mstate, (path_cond | jumping_dst_mstate->get_acc()) );
                    else if(on_stack_[jumping_dst_mstate] && merge_acc_marks(jumping_dst_mstate))
                        return;
                    if(decided_)
                        return;
                }
            }
        }

        if (SCCs_.top() == (src_mstate)) {
            remove_SCC(src_mstate);
        }
    }

    void emptiness_check::tarjan_is_empty(const std::shared_ptr<abstract_successor::mstate> &src_mstate, spot::acc_cond::mark_t path_cond) {
        /// STRONGCONNECT
        // abstr_succ_->print_mstate(src_mstate);
        cnt_++;
        if(abstr_succ_->is_accepting(path_cond) && simulation_prunning(src_mstate))
            return;

        SCCs_.push(src_mstate);
        dfs_acc_stack_.emplace_back(src_mstate, src_mstate->get_acc());
        dfs_num_[src_mstate] = index_;
        index_++;
        tarjan_stack_.push_back(src_mstate);
        on_stack_[src_mstate] = true;

        auto succs = abstr_succ_->get_succs(src_mstate);
        for (auto &dst_mstate: succs) {
            // checking emptiness of curr dst_mstate
            if(empty_lang(dst_mstate))
                continue;

            // init structures
            if(dfs_num_.count(dst_mstate) == 0)
            {
                dfs_num_.insert({dst_mstate, UNDEFINED});
                on_stack_.insert({dst_mstate, false});
            }

            if (dfs_num_[dst_mstate] == UNDEFINED)
            {
                tarjan_is_empty( dst_mstate, (path_cond | dst_mstate->get_acc()) );
                if(decided_)
                    return;
            } else if(on_stack_[dst_mstate] && merge_acc_marks(dst_mstate)) {
                return;
            }
        }

        if (SCCs_.top() == (src_mstate)) {
            remove_SCC(src_mstate);
        }
    }

    void emptiness_check::remove_SCC(const std::shared_ptr<abstract_successor::mstate> & src_mstate) {
        SCCs_.pop();
        std::shared_ptr<abstract_successor::mstate> tmp;

        do {
            tmp = tarjan_stack_.back(); tarjan_stack_.pop_back();
            dfs_acc_stack_.pop_back();
            on_stack_[tmp] = false;
            empty_lang_states_.emplace_back(tmp); // when here, each state has empty language, otherwise we would have ended
        } while (src_mstate != tmp);
    }

    bool emptiness_check::empty_lang(const std::shared_ptr<abstract_successor::mstate> & dst_mstate) {
        if(kofola::OPTIONS.params.count("early_sim") != 0 && kofola::OPTIONS.params["early_sim"] == "yes") {
            for (const auto &empty_state: empty_lang_states_) {
                if (abstr_succ_->subsum_less_early(dst_mstate, empty_state)) {
                    return true;
                }
            }
        }

        if(kofola::OPTIONS.params.count("early_plus_sim") != 0 && kofola::OPTIONS.params["early_plus_sim"] == "yes")
        {
            for (const auto &empty_state: empty_lang_states_) {
                if (abstr_succ_->subsum_less_early_plus(dst_mstate, empty_state)) {
                    return true;
                }
            }
        }

        return false;
    }

    bool emptiness_check::merge_acc_marks(const std::shared_ptr<abstract_successor::mstate> &dst_mstate) {
        spot::acc_cond::mark_t cond = dst_mstate->get_acc();

        // merge acc. marks
        std::shared_ptr<abstract_successor::mstate> tmp;
        do {
            tmp = SCCs_.top(); SCCs_.pop();
            bool root_encountered = dst_mstate->get_encountered();
            if(root_encountered || dfs_num_[tmp] > dfs_num_[dst_mstate])
                cond = (cond |= tmp->get_acc());
            if(abstr_succ_->is_accepting(cond)){
                decided_ = true;
                empty_ = false;
                return true;
            }
        } while(dfs_num_[tmp] > dfs_num_[dst_mstate]);
        dst_mstate->set_encountered(true); // mark visited root
        tmp->set_acc(cond);
        SCCs_.push(tmp);

        return false;
    }

    bool emptiness_check::simulation_prunning(const std::shared_ptr<abstract_successor::mstate> & src_mstate) {
        if(kofola::OPTIONS.params.count("early_sim") != 0 && kofola::OPTIONS.params["early_sim"] == "yes") {
            auto cond = src_mstate->get_acc();
            for (auto it = dfs_acc_stack_.rbegin(); it != dfs_acc_stack_.rend(); ++it) {
                const auto &s = (*it).first;
                if (abstr_succ_->is_accepting(cond) && abstr_succ_->subsum_less_early(s, src_mstate)) {
                    decided_ = true;
                    empty_ = false;
                    return true;
                }
                cond |= (*it).second;
            }
        }

        if(kofola::OPTIONS.params.count("early_plus_sim") != 0 && kofola::OPTIONS.params["early_plus_sim"] == "yes") {
            auto cond1 = src_mstate->get_acc();
            auto cond2 = spot::acc_cond::mark_t();
            for (auto it = dfs_acc_stack_.rbegin(); it != dfs_acc_stack_.rend(); ++it) {
                const auto &s = (*it).first;
                if (abstr_succ_->is_accepting(cond1) && abstr_succ_->is_accepting(cond2) && abstr_succ_->subsum_less_early_plus(s, src_mstate)) {
                    decided_ = true;
                    empty_ = false;
                    return true;
                }

                if(cond1.operator&((*it).second))
                    cond2 |= (*it).second; // already is in cond1, therefore composing cond2
                cond1 |= (*it).second;
            }
        }

        return false;
    }
}// namespace KOFOLA
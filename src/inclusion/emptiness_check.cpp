/**
 * @file emptiness_check.cpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief Implementation of on the fly emptiness check procedure
 * @version 0.1
 * @date 2024-05-03
 * 
 * 
 */

#include "emptiness_check.hpp"
#include <tuple>

// kofola
#include "../complement/complement_tela.hpp"
#include "../util/util.hpp"
#include "../complement/decomposer.hpp"
#include "../complement/complement_sync.hpp"
#include "inclusion_check.hpp"

// Spot
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complete.hh>

// standard library
#include <queue>
#include <utility>

namespace kofola {
    emptiness_check::emptiness_check(inclusion_check *incl_checker):
        incl_checker_(incl_checker),
        dfs_num_(),
        on_stack_(),
        tarjan_stack_(),
        SCCs_(),
        dfs_acc_stack_(),
        empty_lang_states_(),
        state_jumps_to_cutoffs_()
    {
        if(kofola::OPTIONS.params.count("early_sim") != 0 && kofola::OPTIONS.params["early_sim"] == "yes") {
            early_prune_ = true;
        }
        if(kofola::OPTIONS.params.count("early_plus_sim") != 0 && kofola::OPTIONS.params["early_plus_sim"] == "yes") {
            early_prune_ = true;
        }
    }

    bool emptiness_check::empty() {
        entry_states_ = incl_checker_->get_initial_states();
        for (const auto& state: entry_states_) {
            dfs_num_.insert({state, UNDEFINED});
            on_stack_.insert({state, false});
        }

        // Allow entry_states_ to be enriched by gs()/gs_edited
        size_t i = 0;
        auto intersect_aut_acc_code = incl_checker_->get_acc_cond();
        spot::acc_cond code(intersect_aut_acc_code);
        auto fin_mark = intersect_aut_acc_code.fin_unit();

        auto inf_marks_container = intersect_aut_acc_code.inf_unit().sets();
        for(const auto& inf_mark : inf_marks_container) {
            infs_pos_.insert({inf_mark, {}});
        }
        fin_pos_ = {};

        #ifdef ENABLE_COUNTER
            cnt_ = 0;
        #endif

        while (i < entry_states_.size()) {
            const auto& entry = entry_states_[i];
            if (dfs_num_[entry] == UNDEFINED) {
                bool empty;
                #ifdef ENABLE_COUNTER
                    cnt_++;
                #endif

                if(kofola::OPTIONS.params.count("gfee") != 0 && kofola::OPTIONS.params["gfee"] == "yes")
                    empty = gs_edited(entry);
                else {
                    if(code.is_generalized_buchi()) {
                        empty = gs(entry, fin_mark);
                    } else {
                        empty = gen_rabin(entry, fin_mark);
                        // empty = gs(entry, fin_mark);
                    }                     
                } 
                if(!empty){
                    #ifdef ENABLE_COUNTER
                        std::cout << cnt_ << "\n";
                    #endif
                    return false;
                }
            }
            ++i;
        }

        #ifdef ENABLE_COUNTER
            std::cout << cnt_ << "\n";
        #endif
        return true;
    }

    bool emptiness_check::check_simul_less(const std::shared_ptr<inclusion_mstate> &dst_mstate) {
        auto cond = dst_mstate->get_acc();
        spot::acc_cond::mark_t zero{0};
        // traversing 'stack' in reversed order, while according to theorem in thesis, the cond has to be 0
        for (auto it = dfs_acc_stack_.rbegin(); it != dfs_acc_stack_.rend() && cond == zero; ++it) {
            const auto &s = (*it).first;
            if (incl_checker_->subsum_less_early(dst_mstate, s)) {
                for (auto it2 = dfs_acc_stack_.rbegin(); it2 != dfs_acc_stack_.rend() && (*it2).first != s; ++it2) {
                    auto s_between = (*it2).first;
                    // marking destinations of jumps from the states between s and dst_mstate
                    if(state_jumps_to_cutoffs_.count(s_between) == 0) state_jumps_to_cutoffs_.insert({s_between, {dst_mstate}}); // init
                    else state_jumps_to_cutoffs_[s_between].insert(dst_mstate);
                }
                return true;
            }

            cond |= (*it).second;
        }

        return false;
    }

    bool emptiness_check::gs_edited(std::shared_ptr<inclusion_mstate> src_mstate) {
        // stacks to replace recursion
        std::stack<std::shared_ptr<inclusion_mstate>> src_mstates;
        src_mstates.push(nullptr); // default
        std::stack<std::vector<std::shared_ptr<inclusion_mstate>>> successors;
        successors.push({}); // default
        std::stack<spot::acc_cond::mark_t> path_conds;
        path_conds.push({}); // default
        auto succs = incl_checker_->get_succs(src_mstate);
        auto path_cond = spot::acc_cond::mark_t();

        update_structures(src_mstate);

        while(src_mstate != nullptr) {
            // early(+1) simul can decide nonemptiness
            if(early_prune_ && incl_checker_->is_accepting(path_cond) && simulation_prunning(src_mstate))
                return false;

            bool recursion_like = false;
            while(!succs.empty()) {
                auto dst_mstate = succs.back();
                succs.pop_back();

                if(empty_lang(dst_mstate))
                    continue;

                // init structure
                if(dfs_num_.count(dst_mstate) == 0)
                {
                    dfs_num_.insert({dst_mstate, UNDEFINED});
                    on_stack_.insert({dst_mstate, false});
                }
                // dfs_num_ initialised for dst_mstate

                if (dfs_num_[dst_mstate] == UNDEFINED && !check_simul_less(dst_mstate))
                {
                    dfs_acc_stack_.emplace_back(src_mstate, src_mstate->get_acc());
                    // recursion nesting
                    #ifdef ENABLE_COUNTER
                        cnt_++;
                    #endif
                    path_conds.push(path_cond);
                    path_cond |= dst_mstate->get_acc();
                    src_mstates.push(src_mstate);
                    src_mstate = dst_mstate;
                    successors.push(succs);
                    succs = incl_checker_->get_succs(src_mstate);

                    update_structures(src_mstate); 
                    
                    recursion_like = true; // to be able to 'jump'
                    break; 
                } else if(dfs_num_[dst_mstate] != UNDEFINED) {
                    if(on_stack_[dst_mstate] && merge_acc_marks(dst_mstate))
                        return false;

                    // new approach
                    if(state_jumps_to_cutoffs_.count(dst_mstate) == 0) // if there are some jumps
                        continue;
                    for (auto it = state_jumps_to_cutoffs_[dst_mstate].begin(); it != state_jumps_to_cutoffs_[dst_mstate].end();) {
                        auto &jumping_dst_mstate = *it;
                        // same scenarios as for the exploration in original GS alg.
                        if(dfs_num_[jumping_dst_mstate] == UNDEFINED && !check_simul_less(jumping_dst_mstate))
                        {
                            #ifdef ENABLE_COUNTER
                                cnt_++;
                            #endif
                            // recursion nesting
                            path_conds.push(path_cond);
                            path_cond |= jumping_dst_mstate->get_acc();
                            src_mstates.push(src_mstate);
                            src_mstate = jumping_dst_mstate;
                            successors.push(succs);
                            succs = incl_checker_->get_succs(src_mstate);

                            update_structures(src_mstate); 
                            
                            it = state_jumps_to_cutoffs_[dst_mstate].erase(it);
                            recursion_like = true; // to be able to 'jump'
                            break;
                        }
                        else if(on_stack_[jumping_dst_mstate] && merge_acc_marks(jumping_dst_mstate))
                            return false;
                    }

                    if(recursion_like)
                        break;
                    // end of new appraoch
                }
            }

            if(recursion_like)
                continue;

            if (SCCs_.top() == (src_mstate)) {
                remove_SCC(src_mstate);
            }

            // backtracking from recursion
            if(!dfs_acc_stack_.empty())
                dfs_acc_stack_.pop_back();
            src_mstate = src_mstates.top();
            src_mstates.pop();
            succs = successors.top();
            successors.pop();
            path_cond = path_conds.top();
            path_conds.pop();
        }

        return true;
    }

    void emptiness_check::update_structures(const std::shared_ptr<inclusion_mstate>& src_mstate) {
        SCCs_.push(src_mstate);
        dfs_num_[src_mstate] = index_;
        lowlink_[src_mstate] = index_;
        index_++;
        tarjan_stack_.push_back(src_mstate);
        on_stack_[src_mstate] = true;
    }

    bool emptiness_check::gen_rabin(std::shared_ptr<inclusion_mstate> src_mstate, spot::acc_cond::mark_t fin_mark) {
        // stacks to replace recursion
        std::stack<std::shared_ptr<inclusion_mstate>> src_mstates;
        src_mstates.push(nullptr); // default
        std::stack<std::vector<std::shared_ptr<inclusion_mstate>>> successors;
        successors.push({}); // default
        std::stack<spot::acc_cond::mark_t> path_conds;
        path_conds.push({}); // default
        auto succs = incl_checker_->get_succs(src_mstate);
        auto path_cond = spot::acc_cond::mark_t();

        lowlink_[src_mstate] = index_;

        update_structures(src_mstate);

        while(src_mstate != nullptr) {
            // early(+1) simul can decide nonemptiness
            if(early_prune_ && incl_checker_->is_accepting(path_cond) && simulation_prunning(src_mstate))
                return false;
            
            bool recursion_like = false;
            while(!succs.empty()) {
                auto dst_mstate = succs.back();
                // std::cout << "--------------------------------------------------------------------\n";
                // incl_checker_->print_mstate(src_mstate);
                // std::cout << " --to--> \n";
                // incl_checker_->print_mstate(dst_mstate);
                // std::cout << "--------------------------------------------------------------------\n";
                succs.pop_back();

                if(empty_lang(dst_mstate)) {
                    continue;
                }

                // init structure
                if(dfs_num_.count(dst_mstate) == 0)
                {
                    dfs_num_.insert({dst_mstate, UNDEFINED});
                    on_stack_.insert({dst_mstate, false});
                    prefix_.insert({dst_mstate, spot::acc_cond::mark_t()});
                }
                // dfs_num_ initialised for dst_mstate

                if (dfs_num_[dst_mstate] == UNDEFINED)
                {
                    dfs_acc_stack_.emplace_back(src_mstate, src_mstate->get_acc());
                    #ifdef ENABLE_COUNTER
                        cnt_++;
                    #endif
                    path_conds.push(path_cond);
                    path_cond |= dst_mstate->get_acc();
                    src_mstates.push(src_mstate);
                    src_mstate = dst_mstate;
                    successors.push(succs);
                    succs = incl_checker_->get_succs(src_mstate);

                    update_structures(src_mstate); 
                    // mark last positions of 
                    auto marks_container = src_mstate->get_acc().sets();
                    for(const auto& mark: marks_container) {
                        if(infs_pos_.count(mark) != 0) {
                            // store the position (dfs number) for this infinitary mark
                            infs_pos_[mark].push_back(dfs_num_[src_mstate]);
                        } else {
                            fin_pos_.push_back(dfs_num_[src_mstate]); // might be off by one error (which state src or dst??)
                        }
                    }
                    // end of marking
                    
                    recursion_like = true; // to be able to 'jump'
                    break; 
                // } else if(on_stack_[dst_mstate] && merge_acc_marks(dst_mstate)) {
                } else if(on_stack_[dst_mstate] && !(fin_mark & dst_mstate->get_acc())) { // here we can ignore the edge if fin is closing the cycle
                    lowlink_[src_mstate] = std::min(lowlink_[src_mstate], dfs_num_[dst_mstate]);
                    signed most_recent_fin_pos = fin_pos_.empty() ? -1 : fin_pos_.back();
                    bool fin_not_in_cycle = !(dfs_num_[src_mstate] > most_recent_fin_pos && most_recent_fin_pos > dfs_num_[dst_mstate]) || (most_recent_fin_pos == -1);
                    if(fin_not_in_cycle) {
                        auto all_marks = prefix_[dst_mstate] | dst_mstate->get_acc();
                        for (const auto& [inf, all_pos] : infs_pos_) {
                            signed most_recent_inf_pos = all_pos.empty() ? -1 : all_pos.back();
                            bool is_inf_in_cycle = (dfs_num_[src_mstate] > most_recent_inf_pos && most_recent_inf_pos > dfs_num_[dst_mstate]) && (most_recent_inf_pos != -1);
                            if(is_inf_in_cycle) {
                                all_marks.set(inf);
                                // todo early exit if all_marks is accepting
                            }
                        }
                        if(incl_checker_->is_accepting(all_marks)) {
                            return false;
                        }
                    }
                }
            }

            if(recursion_like)
                continue;

            if (SCCs_.top() == (src_mstate)) {
                remove_SCC(src_mstate);
            }
            // backtracking from recursion
            if(!dfs_acc_stack_.empty())
                dfs_acc_stack_.pop_back();

            if(src_mstates.top() == nullptr)
                break;

            auto backtrack_to = src_mstates.top(); 
            // if(dfs_num_[backtrack_to] < dfs_num_[src_mstate] && !(src_mstate->get_acc() & fin_mark)) {
            if(!(src_mstate->get_acc() & fin_mark) && lowlink_[src_mstate] <= lowlink_[backtrack_to]) {
                prefix_[backtrack_to] = (src_mstate->get_acc() | prefix_[src_mstate]);
                lowlink_[backtrack_to] = std::min(lowlink_[src_mstate], lowlink_[backtrack_to]);
            }

            // removes marking
            auto marks_container = src_mstate->get_acc().sets();
            for(const auto& mark: marks_container) {
                if(infs_pos_.count(mark) != 0) {
                    infs_pos_[mark].pop_back();
                } else {
                    fin_pos_.pop_back(); // might be off by one error (which state src or dst??)
                }
            }
            // end of removing

            src_mstate = backtrack_to;
            src_mstates.pop();
            succs = successors.top();
            successors.pop();
            path_cond = path_conds.top();
            path_conds.pop();
        }

        return true;
    }

    bool emptiness_check::gs(std::shared_ptr<inclusion_mstate> src_mstate, spot::acc_cond::mark_t fin_mark) {
        // stacks to replace recursion
        std::stack<std::shared_ptr<inclusion_mstate>> src_mstates;
        src_mstates.push(nullptr); // default
        std::stack<std::vector<std::shared_ptr<inclusion_mstate>>> successors;
        successors.push({}); // default
        std::stack<spot::acc_cond::mark_t> path_conds;
        path_conds.push({}); // default
        auto succs = incl_checker_->get_succs(src_mstate);
        auto path_cond = spot::acc_cond::mark_t();

        update_structures(src_mstate);

        while(src_mstate != nullptr) {
            // early(+1) simul can decide nonemptiness
            if(early_prune_ && incl_checker_->is_accepting(path_cond) && simulation_prunning(src_mstate))
                return false;
            
            bool recursion_like = false;
            while(!succs.empty()) {
                auto dst_mstate = succs.back();
                succs.pop_back();

                if(empty_lang(dst_mstate)) {
                    continue;
                }

                // init structure
                if(dfs_num_.count(dst_mstate) == 0)
                {
                    dfs_num_.insert({dst_mstate, UNDEFINED});
                    on_stack_.insert({dst_mstate, false});
                }
                // dfs_num_ initialised for dst_mstate
                if(fin_mark & dst_mstate->get_acc()) {
                    if(dfs_num_[dst_mstate] == UNDEFINED) {
                        entry_states_.push_back(dst_mstate);
                    }
                    continue;
                }

                if (dfs_num_[dst_mstate] == UNDEFINED)
                {
                    dfs_acc_stack_.emplace_back(src_mstate, src_mstate->get_acc());
                    #ifdef ENABLE_COUNTER
                        cnt_++;
                    #endif
                    path_conds.push(path_cond);
                    path_cond |= dst_mstate->get_acc();
                    src_mstates.push(src_mstate);
                    src_mstate = dst_mstate;
                    successors.push(succs);
                    succs = incl_checker_->get_succs(src_mstate);

                    update_structures(src_mstate); 
                    
                    recursion_like = true; // to be able to 'jump'
                    break; 
                } else if(on_stack_[dst_mstate] && merge_acc_marks(dst_mstate)) {
                    return false;
                }
            }

            if(recursion_like)
                continue;

            if (SCCs_.top() == (src_mstate)) {
                remove_SCC(src_mstate);
            }
            // backtracking from recursion
            if(!dfs_acc_stack_.empty())
                dfs_acc_stack_.pop_back();
            src_mstate = src_mstates.top();
            src_mstates.pop();
            succs = successors.top();
            successors.pop();
            path_cond = path_conds.top();
            path_conds.pop();
        }

        return true;
    }

    void emptiness_check::remove_SCC(const std::shared_ptr<inclusion_mstate> & src_mstate) {
        SCCs_.pop();
        std::shared_ptr<inclusion_mstate> tmp;
        do {
            tmp = tarjan_stack_.back(); tarjan_stack_.pop_back();
            on_stack_[tmp] = false;
            empty_lang_states_[tmp->get_intersect_state().first].emplace_back(tmp); // when here, each state has empty language, otherwise we would have ended
        } while (src_mstate != tmp);
    }

    bool emptiness_check::empty_lang(const std::shared_ptr<inclusion_mstate> & dst_mstate) {
        if(kofola::OPTIONS.params.count("early_sim") != 0 && kofola::OPTIONS.params["early_sim"] == "yes") {
            const auto& col = empty_lang_states_[dst_mstate->get_intersect_state().first];
            if(col.size() > MAX_SUBSUM_BUCKET) {
                return false;
            }
            for (const auto &empty_state: col) {
                if (incl_checker_->subsum_less_early(dst_mstate, empty_state)) {
                    return true;
                }
            }
        }

        if(kofola::OPTIONS.params.count("early_plus_sim") != 0 && kofola::OPTIONS.params["early_plus_sim"] == "yes") {
            const auto& col = empty_lang_states_[dst_mstate->get_intersect_state().first];
            if(col.size() > MAX_SUBSUM_BUCKET) {
                return false;
            }
            for (const auto &empty_state: col) {
                if (incl_checker_->subsum_less_early_plus(dst_mstate, empty_state)) {
                    return true;
                }
            }
        }

        return false;
    }

    bool emptiness_check::merge_acc_marks(const std::shared_ptr<inclusion_mstate> &dst_mstate) {
        spot::acc_cond::mark_t cond = dst_mstate->get_acc();
        std::shared_ptr<inclusion_mstate> tmp;
        do {
            tmp = SCCs_.top(); SCCs_.pop();
            if(dfs_num_[tmp] > dfs_num_[dst_mstate])
                cond |= tmp->get_acc();
            else {
                cond |= tmp->accumulator_;
            }
            if(incl_checker_->is_accepting(cond)){
                decided_ = true;
                empty_ = false;
                return true;
            }
        } while(dfs_num_[tmp] > dfs_num_[dst_mstate]);
        tmp->accumulator_ |= cond;
        SCCs_.push(tmp);

        return false;
    }

    bool emptiness_check::simulation_prunning(const std::shared_ptr<inclusion_mstate> & src_mstate) {
        if(kofola::OPTIONS.params.count("early_sim") != 0 && kofola::OPTIONS.params["early_sim"] == "yes") {
            auto cond = src_mstate->get_acc();
            for (auto it = dfs_acc_stack_.rbegin(); it != dfs_acc_stack_.rend(); ++it) {
                const auto &s = (*it).first;
                // there is a path from s to src_mstate while witnessing acc. cond (and s is simul. < than src_mstate)
                if (incl_checker_->is_accepting(cond) && incl_checker_->subsum_less_early(s, src_mstate)) {
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
                // there is a path from s to src_mstate while witnessing 2 acc. conds (and s is simul. < than src_mstate)
                if (incl_checker_->is_accepting(cond1) && incl_checker_->is_accepting(cond2) && incl_checker_->subsum_less_early_plus(s, src_mstate)) {
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
// Copyright (C) 2017-2019 Laboratoire de Recherche et Développement
// de l'Epita.
// Copyright (C) 2022  The COLA Authors
//
// COLA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// COLA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "elevatorization.hpp"
#include <queue>
#include "complement_sync.hpp"

spot::acc_cond::mark_t get_all_fins_in_dnf(const kofola::CondDNF& dnf) {
    spot::acc_cond::mark_t all_fins{};
    
    for(auto disj : dnf) {
        for(auto& fin : disj.fins) {
            all_fins |= fin;
        }
    }

    return all_fins;
}

kofola::Elevatorization::Elevatorization(const spot::twa_graph_ptr& aut)
    : aut_(aut),
      old_aut_num_states_(aut->num_states()),
      si_(aut_, spot::scc_info_options::ALL),
      support_(old_aut_num_states_),
      compat_(old_aut_num_states_)
{
    // if we work with TELA, we need to properly determine SCC acceptance
    // Spot's is_acceptance might say unknown for Fin conditions
    si_.determine_unknown_acceptance();
    
    // Generate bdd supports and compatible options for each state.
    // Also check if all its transitions are accepting.
    for (unsigned i = 0; i < old_aut_num_states_; ++i) {
        bdd res_support = bddtrue;
        bdd res_compat = bddfalse;
        for (const auto &out: aut->out(i)) {
            res_support &= bdd_support(out.cond);
            res_compat |= out.cond;
        }
        support_[i] = res_support;
        compat_[i] = res_compat;
    }
    // obtain the types of each SCC
    scc_types_ = helpers::get_scc_types(si_);
    // find out the DACs and NACs

    partitions_ = helpers::tnba_complement::create_partitions(si_, kofola::OPTIONS);
    info_ = std::make_unique<kofola::cmpl_info>(
            this->aut_,             // automaton
            std::get<0>(partitions_),         // number of partitions
            std::get<1>(partitions_),       // partition types
            std::get<2>(partitions_),            // state to partition map
            ReachableVector{},// vector of reachable states
            helpers::tnba_complement::create_part_to_scc_map(std::get<3>(partitions_)),        // map of partitions to sets of SCCS they contain
            SCCToSCCSetMap{},   // maps SCCs to the sets of their predecessors
            std::get<4>(partitions_),            // partitions to acceptance condition map
            si_,              // SCC information
            Simulation{},         // direct simulation
            std::vector<bool>{},    // vector for acceptance of states
            false);
}

void kofola::Elevatorization::remove_all_cols_within_scc(size_t scc_idx) 
{
    spot::acc_cond::mark_t empty_mark({});

    auto scc_states = info_->scc_info_.states_of(scc_idx);
    for(auto s : scc_states) {
        for (auto &t : aut_->out(s)) {
            if (s < old_aut_num_states_ && t.dst < old_aut_num_states_ && info_->scc_info_.scc_of(s) == info_->scc_info_.scc_of(t.dst)) {
                t.acc = empty_mark; // remove all colors
            }
        }
    }
}

void kofola::Elevatorization::add_cols_within_scc(size_t part_index, spot::acc_cond::mark_t fins) 
{
    auto scc_states = info_->scc_info_.states_of(part_index);
    for(auto s : scc_states) {
        for (auto &t : aut_->out(s)) {
            if (s < old_aut_num_states_ && t.dst < old_aut_num_states_ && info_->scc_info_.scc_of(s) == info_->scc_info_.scc_of(t.dst)) {
                t.acc = t.acc | fins; // add fins
            }
        }
    }
}

std::set<unsigned> kofola::Elevatorization::get_succ_excluding_colors(
    const std::set<unsigned>& states,
    const bdd& bdd,
    const spot::acc_cond::mark_t& col
) 
{
  std::set<unsigned> succ_states;

  DEBUG_PRINT_LN("states: " + std::to_string(states) + ", letter: " + std::to_string(bdd) + ", col: " + std::to_string(col));
  for (unsigned s : states) {
      for (const auto &t : aut_->out(s)) {
          if (info_->scc_info_.scc_of(s) == info_->scc_info_.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
              if (t.acc & col)
                continue;
              if (t.dst < old_aut_num_states_) { // exclude colors
                    DEBUG_PRINT_LN(std::to_string(s) + " --> " + std::to_string(t.dst) + ", acc: " + std::to_string(t.acc));
                    succ_states.insert(t.dst);
                }
            }
        }
    }

    return succ_states;
}

std::set<unsigned> kofola::Elevatorization::get_succ_including_colors(
    const std::set<unsigned>& states,
    const bdd& bdd,
    const spot::acc_cond::mark_t& col,
    const spot::acc_cond::mark_t& without_cols
) 
{
  std::set<unsigned> succ_states;
  
  DEBUG_PRINT_LN("states: " + std::to_string(states) + ", letter: " + std::to_string(bdd) + ", col: " + std::to_string(col));
    for (unsigned s : states) {
        for (const auto &t : aut_->out(s)) {
            if (info_->scc_info_.scc_of(s) == info_->scc_info_.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
                if (t.acc & without_cols)
                    continue;

                if ((t.acc & col) == col && t.dst < old_aut_num_states_) { // include colors
                    DEBUG_PRINT_LN(std::to_string(s) + " --> " + std::to_string(t.dst) + ", acc: " + std::to_string(t.acc));
                    succ_states.insert(t.dst);
                }
            }
        }
    }

    return succ_states;
}


void kofola::Elevatorization::create_deter_part(size_t scc_idx, kofola::AccClause disjunct) 
{
    disjunct.simplify();

    
    auto scc_states = info_->scc_info_.states_of(scc_idx);
    std::queue<RBL> to_process;
    std::map<RBL, unsigned> det_state_to_spot_state;
    spot::acc_cond::mark_t fins = get_all_fins_in_dnf({disjunct});
    std::vector<spot::acc_cond::mark_t> infs = disjunct.infs;
    unsigned l_mod = infs.size();
    unsigned l_init = 0;
    
    if(l_mod == 0)
        l_mod = 1;

    for(auto s : scc_states) {
        for (auto &t : aut_->out(s)) {
            if (info_->scc_info_.scc_of(s) != info_->scc_info_.scc_of(t.dst) || t.dst >= old_aut_num_states_)
                continue; // only transitions within the SCC

            auto det_state = RBL{{t.dst}, {}, l_init};
            
            
            if(det_state_to_spot_state.find(det_state) == det_state_to_spot_state.end()) {
                auto new_state_num = aut_->new_state();
                det_state_to_spot_state[det_state] = new_state_num;
            }
            
            to_process.push(det_state);
            aut_->new_edge(s, det_state_to_spot_state[det_state], t.cond, {});
        }
    }

    while(!to_process.empty()) {
        RBL current = to_process.front();
        to_process.pop();

        DEBUG_PRINT_LN("Processing det state: " + std::to_string(current) + " as spot state " + std::to_string(det_state_to_spot_state[current]));
        
        // symbols processing
        bdd msupport = bddtrue;
        bdd n_s_compat = bddfalse;
        const std::set<unsigned> &reach_set = current.R;

        // TODO move to disjuncts
        for (unsigned s: reach_set) {
            if(old_aut_num_states_ <= s) {
                continue; // skip newly created states
            }
            msupport &= support_[s];
            n_s_compat |= compat_[s];
        }

        bdd all = n_s_compat;
        
        // iterate over all symbols
        while (all != bddfalse) {
            bdd letter = bdd_satoneset(all, msupport, bddfalse);
            all -= letter;

            DEBUG_PRINT_LN("symbol: " + std::to_string(letter));
            
            auto all_succs_R = get_succ_excluding_colors(current.R, letter, fins);
            auto all_succs_R_visited_inf = all_succs_R;

            if(infs.size() > 0) {
                all_succs_R_visited_inf = get_succ_including_colors(current.R, letter, infs[current.l], fins);
            }

            if(all_succs_R.empty()) {
                continue; // no successors on this letter
            }

            auto all_succs_B = get_succ_excluding_colors(current.B, letter, fins);

            RBL next;
            next.R = all_succs_R;
            next.B = kofola::get_set_union(all_succs_B, all_succs_R_visited_inf);
            next.l = current.l;

            bool emit_acc = false;
            // move to next inf level 
            if (next.R == next.B) {
                next.B.clear();
                next.l = (next.l + 1) % l_mod;
                emit_acc = true;
            }
            
            if(det_state_to_spot_state.find(next) == det_state_to_spot_state.end()) {
                auto new_state_num = aut_->new_state();
                det_state_to_spot_state[next] = new_state_num;
                
                to_process.push(next);
            }
            
            
            if(emit_acc) {
                aut_->new_edge(det_state_to_spot_state[current], det_state_to_spot_state[next], letter, {new_inf_col_, new_fin_col_});
            } else {
                aut_->new_edge(det_state_to_spot_state[current], det_state_to_spot_state[next], letter, {new_fin_col_});
            }
            DEBUG_PRINT_LN("Created transition from " + std::to_string(current) + " to " + std::to_string(next) + " on " + std::to_string(letter) + " as a state " + std::to_string(det_state_to_spot_state[next]) + (emit_acc ? (" with acc " + std::to_string(new_inf_col_)) : "") + " infs: " + std::to_string(infs) + ", fins: " + std::to_string(fins));
        }
    }
}

void kofola::Elevatorization::limit_deter(size_t part_index, size_t scc_idx) 
{
    auto acc = info_->part_to_acc_map_.at(part_index).get_acceptance();
    auto dnf = kofola::cmpl_info::preserve_acc_code_dnf(acc); // TODO check for multiple occurences of one color within one DNF clause 
    
    for(auto disj : dnf) {
        create_deter_part(scc_idx, disj);
    }

    // make sure the nondet. component is non-accepting
    remove_all_cols_within_scc(scc_idx);
    add_cols_within_scc(scc_idx, spot::acc_cond::mark_t{new_fin_col_});
    // should be non-acc. now

    aut_->merge_edges();
}

const spot::twa_graph_ptr& kofola::Elevatorization::elevatorize(bool only_non_buchi) 
{
    auto aut_acc = aut_->get_acceptance();
    auto dnf_aut_acc = kofola::cmpl_info::preserve_acc_code_dnf(aut_acc);

    // ELEVATORIZE nondet. accepting components
    bool elevatorize_needed = false;
    
    auto old_colors_cnt = aut_->acc().num_sets();
    new_inf_col_ = old_colors_cnt;
    new_fin_col_ = old_colors_cnt + 1;
    
    for(unsigned i = 0; i < info_->num_partitions_; i++) {
        if (info_->part_to_type_map_.at(i) != PartitionType::NONDETERMINISTIC)
            continue;
        
        bool is_buchi = info_->part_to_acc_map_.at(i).is_buchi();
        if(is_buchi && only_non_buchi)
            continue;
        
        auto scc_idx = *(info_->part_to_scc_map_.at(i).begin()); // nondet. partitions are singletons
        
        limit_deter(i, scc_idx);
        elevatorize_needed = true;
    }

    std::cerr << "Elevatorization needed: " << (elevatorize_needed ? "yes" : "no") << "\n";

    if(elevatorize_needed) {
        // new acc marks for deter. components (TODO: might use the only one for each originally nonodet. component)
        auto old_acc = aut_->get_acceptance();
        old_acc &= spot::acc_cond::acc_code::fin({new_fin_col_});
        aut_->set_acceptance(old_acc);
        
        old_acc |= spot::acc_cond::acc_code::inf({new_inf_col_});
        aut_->set_acceptance(old_acc);

        aut_->prop_reset(); // elevatorization might violate for instance completeness
    }

    return aut_;
}
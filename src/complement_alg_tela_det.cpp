
#include "complement_alg_tela_det.hpp"
#include "complement_sync.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

std::set<unsigned> union_vecs(std::vector<std::set<unsigned>> sets) {
    std::set<unsigned> result;

    for(unsigned i = 0; i < sets.size(); i++) {
        result = kofola::get_set_union(result, sets[i]);
    }

    return result;
}

/// computes the Cartesian product of a vector of sets (no repetitions
        /// assumed in the inputs)
template<class A>
std::vector<std::vector<A>> compute_cartesian_prod(
        const std::vector<std::vector<A>> vec_of_sets){ // {{{
    const size_t length = vec_of_sets.size();
    std::vector<std::vector<A>> result;

    // this vector will iterate over all possible tuples of indices
    std::vector<size_t> indices(length, 0);

    while (true) {
        std::vector<A> vec;
        for (size_t i = 0; i < length; ++i) {
            assert(indices[i] < vec_of_sets[i].size());
            vec.push_back(vec_of_sets[i][indices[i]]);
        }

        assert(vec.size() == length);
        result.push_back(std::move(vec));

        // generate the next vector of indices, if possible
        bool generated = false;
        for (size_t j = 0; j < length; ++j) {
            ++(indices[j]);
            if (indices[j] < vec_of_sets[j].size()) { // indices is set
                generated = true;
                break;
            } else { // we need to move into the next index
                indices[j] = 0;
            }
        }

        if (!generated) { break; }
    }

    return result;
} // compute_cartesian_prod() }}}

struct PairCompareByContentConj {
    bool operator()(const std::pair<std::shared_ptr<conj_mstate>, unsigned>& a,
                    const std::pair<std::shared_ptr<conj_mstate>, unsigned>& b) const {
        if (*(a.first) < *(b.first)) return true;
        if (*(b.first) < *(a.first)) return false;
        return a.second < b.second;  // tie-breaker
    }
};

struct PairCompareByContentDisj {
    bool operator()(const std::pair<std::shared_ptr<disj_mstate>, unsigned>& a,
                    const std::pair<std::shared_ptr<disj_mstate>, unsigned>& b) const {
        if (*(a.first) < *(b.first)) return true;
        if (*(b.first) < *(a.first)) return false;
        return a.second < b.second;  // tie-breaker
    }
};

conj_mstate::conj_mstate(spot::acc_cond cond) {
    auto all_conjs = cond.top_conjuncts();
    std::vector<spot::acc_cond> disjs;

    for(unsigned i = 0; i < all_conjs.size(); i++) {
        if(all_conjs[i].top_disjuncts().size() == 1) { // atom here
            auto inf = all_conjs[i].inf_unit();
            auto fin = all_conjs[i].fin_unit();

            if(inf.count() != 0) {
                safes_.emplace_back(); // add empty set
                inf_colors_.emplace_back(inf.min_set() - 1); // min_set returns incremented value 
                infs_ = true;
            } else if(fin.count() != 0) {
                m_check_.emplace_back(); // add empty set
                fin_colors_.emplace_back(fin.min_set() - 1); // min_set returns incremented value 
                fins_ = true;
            }
        } else {
            disjs.emplace_back(all_conjs[i]);
        }
    }

    // inner disjuncts
    for(auto disj: disjs) {
        std::shared_ptr<disj_mstate> tmp(new disj_mstate(disj));
        disjuncts_.emplace_back(tmp);
    }
}

std::shared_ptr<conj_mstate> conj_mstate::clone() const {
    // Deep copy disjuncts
    std::vector<std::shared_ptr<disj_mstate>> copied_disjuncts;
    for (const auto& d : disjuncts_) {
        copied_disjuncts.push_back(d->clone());
    }

    return std::make_shared<conj_mstate>(
            check_,
            safes_,
            m_check_,
            breakpoint_,
            copied_disjuncts,
            inf_colors_,
            fin_colors_,
            infs_,
            fins_,
            rr_pointer_,
            active_
    );
}

std::set<unsigned> conj_mstate::get_all_states() {
    std::set<unsigned> result;
    auto all_safes = union_vecs(safes_);
    auto all_Ms = union_vecs(m_check_);

    result = kofola::get_set_union(check_, all_safes);
    result = kofola::get_set_union(result, all_Ms);

    for(auto disj: disjuncts_) {
        result = kofola::get_set_union(result, disj->get_all_states());
    }

    return result;
}

std::string conj_mstate::to_str() {
    std::string result = std::to_string("[ (") + ((this->active_)? "Act" : "Pas") + std::to_string(") ");

    // Add check
    result += "C=" + std::to_string(check_);

    // Add safes
    for (size_t i = 0; i < safes_.size(); ++i) {
        result += ",S" + std::to_string(i) + "=" + std::to_string(safes_[i]);
    }

    // Add m_check
    for (size_t i = 0; i < m_check_.size(); ++i) {
        result += ",M" + std::to_string(i) + "=" + std::to_string(m_check_[i]);
    }

    // Add breakpoint
    result += ",B=" + std::to_string(breakpoint_);

    for(size_t i = 0; i < disjuncts_.size(); i++) {
        result += ",Disj" + std::to_string(i) + "=[" + disjuncts_[i]->to_str() + "]";
    }

    result += " | " + std::to_string(rr_pointer_) + "]";
    return result; 
}

void conj_mstate::passivate() {
    breakpoint_ = {};
    active_ = false;
    for(unsigned i = 0; i < disjuncts_.size(); i++) {
        disjuncts_[i]->passivate();
    }
}

void conj_mstate::activate() {
    rr_pointer_ = 0;
    active_ = true;
    if(infs_ || fins_) {
        breakpoint_ = kofola::get_set_union(check_, union_vecs(m_check_));

    } else {
        if(disjuncts_.size() != 0) {
            disjuncts_[0]->activate();
        }
    }
}

void conj_mstate::move_rr_ptr() {
    unsigned offset = 0;
    if(infs_ || fins_)
        offset = 1;

    if(rr_pointer_ == 0 && (fins_ || infs_)) {
        breakpoint_ = {};
    } else {
        if(disjuncts_.size() != 0)
            disjuncts_[rr_pointer_ - offset]->passivate();
    }

    rr_pointer_ = (rr_pointer_ + 1) % (offset + disjuncts_.size());

    if(rr_pointer_ == 0 && (fins_ || infs_)) {
        breakpoint_ = kofola::get_set_union(check_, union_vecs(m_check_));
    } else {
        if(disjuncts_.size() != 0)
            disjuncts_[rr_pointer_ - offset]->activate();
    }
}

unsigned conj_mstate::get_rr_ptr() {
    return rr_pointer_;
}

bool contains_outgoing_transitions_in_scc_given_color(
  const spot::const_twa_graph_ptr&    aut,
  const spot::scc_info&               scc_info,
  const std::set<unsigned>&           states,
  const bdd&                          symbol,
  const spot::acc_cond::mark_t        color)
{ // {{{
  for (unsigned s : states) {
    for (const auto &t : aut->out(s)) {
      if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(symbol, t.cond)) {
        if (t.acc == color) { return true; }
      }
    }
  }

  return false;
}

template <class T>
  std::set<unsigned> get_all_successors_in_scc_without_color(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const T&                          current_states,
    const bdd&                        symbol,
    const spot::acc_cond::mark_t      colors)
  { // {{{
    std::set<unsigned> successors;

    for (unsigned s : current_states) {
      for (const auto &t : aut->out(s)) {
        if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(symbol, t.cond) && !(t.acc & colors)) {
          successors.insert(t.dst); }
      }
    }

    return successors;
  } // get_all_successors_in_scc() }}}

std::vector<std::vector<std::set<unsigned>>> nondeter_scatter(std::vector<std::set<unsigned>> safes, std::vector<unsigned> runs) {
    unsigned runs_cnt = runs.size();
    std::vector<unsigned> run_map(runs.size(), 0);
    unsigned maxVal = safes.size() - 1;

    std::vector<std::vector<std::set<unsigned>>> result;

    while (true) {
        std::vector<std::set<unsigned>> new_safes = safes;
        for(unsigned i = 0; i < runs_cnt; i++) {
            new_safes[run_map[i]].insert(runs[i]);
        }
        result.emplace_back(new_safes);

        // Increment the vector like a number in base (maxVal+1)
        int idx = runs_cnt - 1;
        while (idx >= 0) {
            if (run_map[idx] < maxVal) {
                run_map[idx]++;
                break;
            } else {
                run_map[idx] = 0;
                idx--;
            }
        }

        // If we tried to increment beyond the first element, we're done
        if (idx < 0)
            break;
    }

    return result;
}

bool contains(
        const std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>>& vec,
        const std::pair<std::shared_ptr<conj_mstate>, unsigned>& value
) {
    for (const auto& item : vec) {
        if (*item.first == *value.first && item.second == value.second) {
            return true;
        }
    }
    return false;
}

bool contains(
        const std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned>>& vec,
        const std::pair<std::shared_ptr<disj_mstate>, unsigned>& value
) {
    for (const auto& item : vec) {
        if (*item.first == *value.first && item.second == value.second) {
            return true;
        }
    }
    return false;
}

std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> conj_mstate::succs(
    const std::vector<unsigned>&  new_runs,
    const bdd&                  symbol,
    const kofola::cmpl_info& info) {
    // NCSB
    std::vector<std::set<unsigned>> S_nexts;
    std::set<unsigned> S_nexts_all;

    for(unsigned i = 0; i < safes_.size(); i++) {
        auto S = safes_[i];
        if(contains_outgoing_transitions_in_scc_given_color(info.aut_, info.scc_info_, safes_[i], symbol, spot::acc_cond::mark_t{inf_colors_[i]})) {
            return {}; // violated safe runs
        }
        auto S_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, safes_[i], symbol);
        S_nexts.emplace_back(S_next);
        S_nexts_all = kofola::get_set_union(S_nexts_all, S_next);
    }

    auto C_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, check_, symbol);
    C_next = kofola::get_set_difference(C_next, S_nexts_all);

    auto B_from_NCSB = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, kofola::get_set_intersection(breakpoint_, check_), symbol);

    std::vector<unsigned> B_from_NCSB_vec(B_from_NCSB.begin(), B_from_NCSB.end());
    auto guessing_check = kofola::get_set_difference(C_next, B_from_NCSB);
    auto guessing_break = guessing_check;
    auto all_guesing_safes = nondeter_scatter(S_nexts, B_from_NCSB_vec);
    if(rr_pointer_ > 1)
        all_guesing_safes = {};

    // MH
    std::vector<std::set<unsigned>> m_check_nexts;

    for(unsigned i = 0; i < m_check_.size(); i++) {
        auto m_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, m_check_[i], symbol);
        m_check_nexts.emplace_back(m_next);
    }
    auto B_from_MH = get_all_successors_in_scc_without_color(info.aut_, info.scc_info_, kofola::get_set_difference(breakpoint_, check_), symbol, spot::acc_cond::mark_t(fin_colors_.begin(), fin_colors_.end()));
    // final 
    unsigned runs_cnt = new_runs.size();
    std::vector<unsigned> run_map(new_runs.size(), 0);
    unsigned maxVal = 1 + disjuncts_.size();

//    std::set<std::pair<std::shared_ptr<conj_mstate>, unsigned>, PairCompareByContentConj> result;
    std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> result;

    std::vector<unsigned> result_acc;
    while (true) {
        std::set<unsigned> enrich_C;
        std::vector<unsigned> enrich_M;
        std::vector<std::set<unsigned>> new_runs_disj(disjuncts_.size());

        unsigned scattered = 0;

        for(unsigned i = 0; i < runs_cnt; i++) {
            if(run_map[i] == 0) { // this is not optimal at all - generating for infs if there are none
                if(infs_) {
                    enrich_C.insert(new_runs[i]);
                    scattered++;
                }
            } else if(run_map[i] == 1) {
                if(fins_) {
                    enrich_M.emplace_back(new_runs[i]);
                    scattered++;
                }
            } else {
                if(disjuncts_.size() > 0) {
                    new_runs_disj[run_map[i] - 2].insert(new_runs[i]);
                    scattered++;
                }
            }
        }

        // Increment the vector like a number in base (maxVal+1)
        int idx = runs_cnt - 1;
        while (idx >= 0) {
            if (run_map[idx] < maxVal) {
                run_map[idx]++;
                break;
            } else {
                run_map[idx] = 0;
                idx--;
            }
        }

        if(scattered != runs_cnt)
            continue;

        auto all_guesing_Ms = nondeter_scatter(m_check_nexts, enrich_M);

        auto B = kofola::get_set_union(B_from_MH,B_from_NCSB);
        auto C = kofola::get_set_union(C_next, enrich_C);
        bool move_ncsbm_ptr = false;
        if(B.empty() && (infs_ || fins_) && rr_pointer_ == 0 && active_) {
            move_ncsbm_ptr = true;
        }
        
        std::vector<std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned>>> new_disjs;
        for(unsigned i = 0; i < disjuncts_.size(); i++) {
            auto new_disj = disjuncts_[i]->succs(std::vector<unsigned>(new_runs_disj[i].begin(), new_runs_disj[i].end()),symbol,info);
            if(new_disj.empty()) {
                return {};
            }
            new_disjs.emplace_back(new_disj);
        }

        auto prod_of_disj_succs = compute_cartesian_prod(new_disjs);
        if(prod_of_disj_succs.size() == 0) // disjuncts_.size() == 0
            prod_of_disj_succs.emplace_back();

        if(!fins_)
            all_guesing_Ms = {{}}; // one empty vector to allow one cycle

        for(unsigned i = 0; i < all_guesing_Ms.size(); i++) {
            auto Ms_flatt = union_vecs(all_guesing_Ms[i]);
            if(move_ncsbm_ptr == true) { // this is redundant (to resample B)
                B = kofola::get_set_union(C,Ms_flatt);
            }

            for(unsigned j = 0; j < prod_of_disj_succs.size(); j++) {
                auto tmp_disj = prod_of_disj_succs[j];

                // extract disjs and acc.
                std::vector<std::shared_ptr<disj_mstate>> succ_disj;
                std::vector<unsigned> accs;
                for(auto next_pair: tmp_disj) {
                    succ_disj.emplace_back(next_pair.first->clone());
                    accs.emplace_back(next_pair.second);
                }
                
                std::shared_ptr<conj_mstate> tmp(new conj_mstate(C, S_nexts, all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_, active_));
                bool move_disj_ptr = false;
                if(rr_pointer_ >= (infs_ || fins_) && accs[rr_pointer_ - (infs_ || fins_)] == 1) {
                    move_disj_ptr = true;
                }
                
                if((move_ncsbm_ptr || move_disj_ptr) && active_) tmp->move_rr_ptr();

                unsigned curr_acc = 0;
                if(tmp->get_rr_ptr() == 0 && (move_ncsbm_ptr || move_disj_ptr) && active_)
                    curr_acc = 1;
                if(!contains(result, {tmp, curr_acc}))
                    result.emplace_back(tmp, curr_acc);

            }
            // }
                
            // guessed
            for(unsigned k = 0; k < all_guesing_safes.size() && active_ && rr_pointer_ == 0 && infs_ && !B_from_NCSB.empty(); k++) {
                auto B = kofola::get_set_union(B_from_MH, guessing_break);
                auto C = kofola::get_set_union(guessing_check, enrich_C);

                B = kofola::get_set_union(C,Ms_flatt); // also redundant (since moving rr_ptr)

                for(unsigned j = 0; j < prod_of_disj_succs.size(); j++) {
                    auto tmp_disj = prod_of_disj_succs[j];
                
                    // extract disjs and acc.
                    std::vector<std::shared_ptr<disj_mstate>> succ_disj;
                    for(auto next_pair: tmp_disj) {
                        succ_disj.emplace_back(next_pair.first->clone());
                    }

                    std::shared_ptr<conj_mstate> tmp2(new conj_mstate(C, all_guesing_safes[k], all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_, active_));
                    
                    tmp2->move_rr_ptr();
                    unsigned curr_acc = 0;
                    if(tmp2->get_rr_ptr() == 0)
                        curr_acc = 1;

                    if(!contains(result, {tmp2, curr_acc}))
                        result.emplace_back(tmp2, curr_acc);
                }
                // }
            }
        }
        // If we tried to increment beyond the first element, we're done
        if (idx < 0)
            break;
    }

//    std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned >> v(result.begin(),result.end());
    return result;
}

bool conj_mstate::operator==(const conj_mstate& other) const {
    return (
        this->check_ == other.check_ &&
        this->safes_ == other.safes_ &&
        this->m_check_ == other.m_check_ &&
        this->breakpoint_ == other.breakpoint_ &&
        this->disjuncts_ == other.disjuncts_ &&
        this->inf_colors_ == other.inf_colors_ &&
        this->fin_colors_ == other.fin_colors_ &&
        this->rr_pointer_ == other.rr_pointer_ &&
        this->active_ == other.active_
    );
}

bool conj_mstate::operator<(const conj_mstate& other) const {

    if (this->active_ != other.active_) return this->active_ < other.active_;
    if (this->inf_colors_ != other.inf_colors_) return this->inf_colors_ < other.inf_colors_;
    if (this->fin_colors_ != other.fin_colors_) return this->fin_colors_ < other.fin_colors_;
    if (this->rr_pointer_ != other.rr_pointer_) return this->rr_pointer_ < other.rr_pointer_;
    if (this->check_ != other.check_) return this->check_ < other.check_;
    if (this->safes_ != other.safes_) return this->safes_ < other.safes_;
    if (this->m_check_ != other.m_check_) return this->m_check_ < other.m_check_;
    if (this->breakpoint_ != other.breakpoint_) return this->breakpoint_ < other.breakpoint_;

    const size_t length = this->disjuncts_.size();
    for (size_t i = 0; i < length; ++i) {
        if ( !(*(this->disjuncts_[i]) == *(other.disjuncts_[i])) ) {
            return *(this->disjuncts_[i]) < *(other.disjuncts_[i]);
        }
    }

    return false; // all equal
}

disj_mstate::disj_mstate(spot::acc_cond cond) {
    auto all_disjs = cond.top_disjuncts();
    std::vector<spot::acc_cond> conjs;

    for(unsigned i = 0; i < all_disjs.size(); i++) {
        if(all_disjs[i].top_conjuncts().size() == 1) { // atom here
            auto inf = all_disjs[i].inf_unit();
            auto fin = all_disjs[i].fin_unit();

            if(inf.count() != 0) {
                // safes_ = {}; // add empty set
                inf_colors_.emplace_back(inf.min_set() - 1); // min_set returns incremented value
                infs_ = true;
            } else if(fin.count() != 0) {
                fin_colors_.emplace_back(fin.min_set() - 1); // min_set returns incremented value
                fins_ = true;
            }
        } else {
            conjs.emplace_back(all_disjs[i]);
        }
    }

    // inner conjuncts
    for(auto conj: conjs) {
        std::shared_ptr<conj_mstate> tmp(new conj_mstate(conj));
        conjuncts_.emplace_back(tmp);
    }
}

std::shared_ptr<disj_mstate> disj_mstate::clone() const {
    // Deep copy conjuncts
    std::vector<std::shared_ptr<conj_mstate>> copied_conjuncts;
    for (const auto& c : conjuncts_) {
        copied_conjuncts.push_back(c->clone());
    }

    return std::make_shared<disj_mstate>(
            check_,
            safes_,
            breakpoint_,
            copied_conjuncts,
            inf_colors_,
            fin_colors_,
            infs_,
            fins_,
            rr_pointer_,
            active_
    );
}

std::set<unsigned> disj_mstate::get_all_states() {
    std::set<unsigned> result;

    result = kofola::get_set_union(check_, safes_);

    for(auto conj: conjuncts_) {
        result = kofola::get_set_union(result, conj->get_all_states());
    }

    return result;
}

std::string disj_mstate::to_str() {
    std::string result = std::to_string("[ (") + ((this->active_)? "Act" : "Pas") + std::to_string(") ");

    // Add check
    result += "C=" + std::to_string(check_);

    // Add safes
    result += ",S=" + std::to_string(safes_);

    // Add breakpoint
    result += ",B=" + std::to_string(breakpoint_);

    for(size_t i = 0; i < conjuncts_.size(); i++) {
        result += ",Conj" + std::to_string(i) + "=[" + conjuncts_[i]->to_str() + "]";
    }

    result += " | " + std::to_string(rr_pointer_) + "]";
    return result; 
}

void disj_mstate::passivate() {
    breakpoint_ = {};
    active_ = false;
    rr_pointer_ = 0;
    for(unsigned i = 0; i < conjuncts_.size(); i++) {
        conjuncts_[i]->passivate();
    }
}

void disj_mstate::activate() {
    rr_pointer_ = 0;
    active_ = true;
    if(infs_) {
        breakpoint_ = check_;
    } 
    else if(fins_) {
        breakpoint_ = kofola::get_set_union(check_, safes_);
    }
    else {
        if(conjuncts_.size() != 0) {
            conjuncts_[0]->activate();
        }
    }
}

void disj_mstate::move_rr_ptr() {
    unsigned infs_and_fins_cnt = inf_colors_.size() + fin_colors_.size();

    if(rr_pointer_ < inf_colors_.size()) {
        breakpoint_ = {};
        rr_pointer_ = inf_colors_.size() - 1;
    } else if(rr_pointer_ < infs_and_fins_cnt) {
        rr_pointer_ = infs_and_fins_cnt - 1;
    } else {
        if(conjuncts_.size() != 0)
            conjuncts_[rr_pointer_ - infs_and_fins_cnt]->passivate();
    }
    rr_pointer_ = (rr_pointer_ + 1) % (infs_and_fins_cnt + conjuncts_.size());
    if(rr_pointer_ == 0 && infs_) {
        breakpoint_ = check_;
    } else if(rr_pointer_ == 0 && fins_) {
        breakpoint_ = kofola::get_set_union(check_, safes_);
    } else {
        if(conjuncts_.size() != 0)
            conjuncts_[rr_pointer_ - infs_and_fins_cnt]->activate();
    }
}

unsigned disj_mstate::get_rr_ptr() {
    return rr_pointer_;
}
std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned>> disj_mstate::succs(const std::vector<unsigned>&  new_runs,
    const bdd&                 symbol,
    const kofola::cmpl_info& info) {
    // NCSB
    // check violation of safe
    for(unsigned i = 0; i < inf_colors_.size(); i++) {
        if(contains_outgoing_transitions_in_scc_given_color(info.aut_, info.scc_info_, safes_, symbol, spot::acc_cond::mark_t{inf_colors_[i]})) {
            return {}; // violated safe runs
        }
    }
    auto S_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, safes_, symbol);
    auto C_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, check_, symbol);
    C_next = kofola::get_set_difference(C_next, S_next);
    C_next = kofola::get_set_union(C_next, std::set<unsigned>(new_runs.begin(), new_runs.end()));
    auto B_next = kofola::get_all_successors_in_scc(info.aut_, info.scc_info_, breakpoint_, symbol);

    auto guessing_C = kofola::get_set_difference(C_next, B_next);
    auto guessing_B = guessing_C;
    auto guessing_S = kofola::get_set_union(S_next, B_next);

    // MH
    std::set<unsigned> B_from_MH = {};
    std::set<unsigned> C_from_MH = {};
    if(fins_) {
        auto B_from_MH = get_all_successors_in_scc_without_color(info.aut_, info.scc_info_, breakpoint_, symbol, spot::acc_cond::mark_t{fin_colors_[rr_pointer_]});
        auto C_from_MH = kofola::get_set_union(C_next, S_next);
    }

    // final 
    std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned>> result;
    std::vector<unsigned> result_acc;

    std::vector<std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>>> new_conjs;
    for(unsigned i = 0; i < conjuncts_.size(); i++) {
        auto new_conj = conjuncts_[i]->succs(new_runs,symbol,info);
        if(new_conj.empty()) {
            return {};
        }
        new_conjs.emplace_back(new_conj);
    }

    auto prod_of_conj_succs = compute_cartesian_prod(new_conjs);
    if(prod_of_conj_succs.size() == 0)
        prod_of_conj_succs.emplace_back();

    for(unsigned i = 0; i < prod_of_conj_succs.size(); i++) {
         auto tmp_conj = prod_of_conj_succs[i];
                
        // extract disjs and acc.
        std::vector<std::shared_ptr<conj_mstate>> succ_conj;
        std::vector<std::shared_ptr<conj_mstate>> succ_conj2;
        std::vector<unsigned> accs;
        for(auto next_pair: tmp_conj) {
            succ_conj.emplace_back(next_pair.first->clone());
            succ_conj2.emplace_back(next_pair.first->clone());
            accs.emplace_back(next_pair.second);
        }
        
        std::set<unsigned> C;
        std::set<unsigned> B;

        unsigned infs_and_fins_cnt = inf_colors_.size() + fin_colors_.size();
        bool move_ptr = false;

        if(rr_pointer_ < inf_colors_.size()) {
            C = C_next;
            B = B_next;
            move_ptr = (B.size() == 0);
        }
        else if (rr_pointer_ < infs_and_fins_cnt) {
            C = C_from_MH;
            B = B_from_MH;
            move_ptr = (B.size() == 0);
        } else {
            C = C_from_MH; 
            B = {};
        }
        
        std::shared_ptr<disj_mstate> tmp(new disj_mstate(C, S_next, B, succ_conj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_, active_));
        
        if(rr_pointer_ >= infs_and_fins_cnt && accs[rr_pointer_ - infs_and_fins_cnt] == 1 && active_) {
            move_ptr = true;
        }

        if(move_ptr && active_) tmp->move_rr_ptr();

        unsigned curr_acc = 0;
        if(tmp->get_rr_ptr() == 0 && move_ptr && active_)
            curr_acc = 1;

        if(!contains(result, {tmp, curr_acc}))
            result.emplace_back(tmp, curr_acc);

        // guessing
        if(rr_pointer_ >= inf_colors_.size() || !active_ || !infs_)
            continue;
        
        std::shared_ptr<disj_mstate> tmp2(new disj_mstate(guessing_C, guessing_S, guessing_B, succ_conj2, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_, active_));
        
        tmp2->move_rr_ptr();

        curr_acc = 0;
        if(tmp2->get_rr_ptr() == 0)
            curr_acc = 1;

        if(!contains(result, {tmp2, curr_acc}))
            result.emplace_back(tmp2, curr_acc);
    }

//    std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned >> v(result.begin(),result.end());
    return result;
}

bool disj_mstate::operator==(const disj_mstate& other) const {
    return (
        this->check_ == other.check_ &&
        this->safes_ == other.safes_ &&
        this->breakpoint_ == other.breakpoint_ &&
        this->conjuncts_ == other.conjuncts_ &&
        this->inf_colors_ == other.inf_colors_ &&
        this->fin_colors_ == other.fin_colors_ &&
        this->rr_pointer_ == other.rr_pointer_ &&
        this->active_ == other.active_
    );
}

bool disj_mstate::operator<(const disj_mstate& other) const {
    if (this->active_ != other.active_) return this->active_ < other.active_;
    if (this->inf_colors_ != other.inf_colors_) return this->inf_colors_ < other.inf_colors_;
    if (this->fin_colors_ != other.fin_colors_) return this->fin_colors_ < other.fin_colors_;
    if (this->rr_pointer_ != other.rr_pointer_) return this->rr_pointer_ < other.rr_pointer_;
    if (this->check_ != other.check_) return this->check_ < other.check_;
    if (this->safes_ != other.safes_) return this->safes_ < other.safes_;
    if (this->breakpoint_ != other.breakpoint_) return this->breakpoint_ < other.breakpoint_;

    const size_t length = this->conjuncts_.size();
    for (size_t i = 0; i < length; ++i) {
        if ( !(*(this->conjuncts_[i]) == *(other.conjuncts_[i])) ) {
            return *(this->conjuncts_[i]) < *(other.conjuncts_[i]);
        }
    }

    return false; // all equal
}

namespace { // {{{
    /// partial macrostate for the given component
    class mstate_tela_det : public abstract_complement_alg::mstate
    { // {{{
    private: // DATA MEMBERS

    std::shared_ptr<conj_mstate> conj_;
    bool active_;
    std::set<unsigned> breakpoint_;

    public: // METHODS

    /// constructor
    mstate_tela_det(spot::acc_cond cond) : conj_(std::make_shared<conj_mstate>(cond))
    { 
    }

    mstate_tela_det(std::shared_ptr<conj_mstate> conj) : conj_(conj)
    { 
    }

    virtual std::string to_string() const override;
    virtual bool is_active() const override { return this->conj_->get_activity(); }
    virtual bool eq(const mstate& rhs) const override;
    virtual bool lt(const mstate& rhs) const override;
    virtual ~mstate_tela_det() override { }

    virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
    virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override {  }

    friend class kofola::complement_tela_det;
    }; // mstate_tela_det }}}


    std::string mstate_tela_det::to_string() const
    {
        std::string res = std::string("[TELA_DET(") + ((this->conj_->get_activity())? "A" : "T") + "): ";
        res += conj_->to_str() + "]";

        return res;
    }

    bool mstate_tela_det::eq(const mstate& rhs) const
    {
    const mstate_tela_det* rhs_mh = dynamic_cast<const mstate_tela_det*>(&rhs);
    assert(rhs_mh);
    return ((*this->conj_) == (*rhs_mh->conj_));
    }

    bool mstate_tela_det::lt(const mstate& rhs) const
    {
    const mstate_tela_det* rhs_mh = dynamic_cast<const mstate_tela_det*>(&rhs);
    assert(rhs_mh);

    return (*this->conj_) < (*rhs_mh->conj_);
    }

    } // anonymous namespace }}}

    complement_tela_det::complement_tela_det(const cmpl_info& info, unsigned part_index)
    : abstract_complement_alg(info, part_index)
    {  }

    mstate_set complement_tela_det::get_init()
    { // {{{
        std::set<unsigned> init_state;

        unsigned orig_init = this->info_.aut_->get_init_state_number();

        if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
            init_state.insert(orig_init);
        }

        mstate_set result;

        std::shared_ptr<mstate_tela_det> ms(new mstate_tela_det(info_.aut_->acc()));
        if(init_state.size() == 0) {
            return {ms};
        }

        auto all_inits = ms->conj_->succs(std::vector<unsigned>(init_state.begin(),init_state.end()), bddtrue, this->info_);
        for(auto init: all_inits) {
            std::shared_ptr<mstate_tela_det> ms_succ(new mstate_tela_det(init.first));
            ms_succ->conj_->passivate();
//            ms_succ->conj_->activate();
            std::cout << init.first->to_str() << "\n";
            result.push_back(ms_succ);
        }

        return result;
    } // get_init() }}}

    mstate_col_set complement_tela_det::get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol)
    { // {{{
        const mstate_tela_det* src_tela_det = dynamic_cast<const mstate_tela_det*>(src);
        assert(src_tela_det->conj_->get_activity());

        auto tmp = kofola::get_set_difference(glob_reached,src_tela_det->conj_->get_all_states());
        std::vector<unsigned> new_runs;
        for(auto state: tmp) {
            if(this->info_.st_to_part_map_.at(state) == static_cast<int>(this->part_index_)) {
                new_runs.emplace_back(state);
            }
        }

        auto succ_states = src_tela_det->conj_->succs(new_runs,symbol,this->info_);
        mstate_col_set res;
        for(auto suc: succ_states) {
            suc.first->passivate();
            std::shared_ptr<mstate> ms(new mstate_tela_det(suc.first));
            if(suc.second == 0)
                res.emplace_back(ms, std::set<unsigned>{});
            else
                res.emplace_back(std::make_pair(ms, std::set<unsigned>{suc.second}));
        }

        return res;
    }
    
    
    mstate_set complement_tela_det::lift_track_to_active(const mstate* src) {
        const mstate_tela_det* src_tela_det = dynamic_cast<const mstate_tela_det*>(src);

        auto tmp = src_tela_det->conj_->clone();
        tmp->activate();

        std::shared_ptr<mstate> ms(new mstate_tela_det(tmp));

        return {ms};
    }
    
    // lift_track_to_active() }}}

    mstate_col_set complement_tela_det::get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol,
    bool resample)
    {
        const mstate_tela_det* src_tela_det = dynamic_cast<const mstate_tela_det*>(src);

        auto in_src = kofola::get_all_successors_in_scc(this->info_.aut_, this->info_.scc_info_, src_tela_det->conj_->get_all_states(), symbol);
        auto tmp = kofola::get_set_difference(glob_reached,in_src);
        std::vector<unsigned> new_runs;
        for(auto state: tmp) {
            if(this->info_.st_to_part_map_.at(state) == static_cast<int>(this->part_index_)) {
                new_runs.emplace_back(state);
            }
        }

        auto succ_states = src_tela_det->conj_->succs(new_runs,symbol,this->info_);
        mstate_col_set res;
        for(auto suc: succ_states) {
            std::shared_ptr<mstate> ms(new mstate_tela_det(suc.first));
            if(suc.second == 0)
                res.emplace_back(ms, std::set<unsigned>{});
            else
                res.emplace_back(std::make_pair(ms, std::set<unsigned>{0}));
        }

        return res;
    }

complement_tela_det::~complement_tela_det()
{ }
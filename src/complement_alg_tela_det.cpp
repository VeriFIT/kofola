
#include "complement_alg_tela_det.hpp"
#include "complement_sync.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

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
            } else if(fin.count() != 0) {
                m_check_.emplace_back(); // add empty set
                fin_colors_.emplace_back(fin.min_set() - 1); // min_set returns incremented value 
            }
        } else {
            disjs.emplace_back(all_conjs[i]);
        }
    }

    // inner disjuncts
    for(auto disj: disjs) {
        disjuncts_.emplace_back(disj_mstate(disj));
    }
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
    if(rr_pointer_ == 0 && (fins_ || infs_)) {
        breakpoint_ = {};
    } else {
        if(disjuncts_.size() != 0)
            disjuncts_[rr_pointer_]->passivate();
    }
    
    unsigned offset = 0;
    if(infs_ || fins_)
        offset = 1;
    rr_pointer_ = (rr_pointer_ + 1) % (offset + disjuncts_.size());

    if(rr_pointer_ == 0 && (fins_ || infs_)) {
        breakpoint_ = kofola::get_set_union(check_, union_vecs(m_check_));
    } else {
        if(disjuncts_.size() != 0)
            disjuncts_[rr_pointer_]->activate();
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

std::set<unsigned> union_vecs(std::vector<std::set<unsigned>> sets) {
    std::set<unsigned> result;

    for(unsigned i = 0; i < sets.size(); i++) {
        result = kofola::get_set_union(result, sets[i]);
    }

    return result;
}

std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> conj_mstate::succs(
    const std::vector<unsigned>&  new_runs,
    const bdd&                  symbol,
    kofola::cmpl_info& info) {

    // NCSB
    std::vector<std::set<unsigned>> S_nexts;
    std::set<unsigned> S_nexts_all;

    for(unsigned i = 0; i < safes_.size(); i++) {
        auto S = safes_[i];
        if(contains_outgoing_transitions_in_scc_given_color(info.aut_, info.scc_info_, safes_[i], symbol, spot::acc_cond::mark_t(inf_colors_[i]))) {
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
    unsigned maxVal = 2 + disjuncts_.size();

    std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> result;
    std::vector<unsigned> result_acc;

    while (true) {
        std::set<unsigned> enrich_C;
        std::vector<unsigned> enrich_M;
        std::vector<std::set<unsigned>> new_runs_disj;
        for(unsigned i = 0; i < runs_cnt; i++) {
            if(run_map[i] == 0) { // this is not optimal at all - generating for infs if there are none
                if(infs_)
                    enrich_C.insert(new_runs[i]);
            } else if(run_map[i] == 1) {
                if(fins_)
                    enrich_M.emplace_back(new_runs[i]);
            } else {
                new_runs_disj[run_map[i]].insert(new_runs[i]);
            }
        }

        auto all_guesing_Ms = nondeter_scatter(m_check_nexts, enrich_M);

        auto B = kofola::get_set_union(B_from_MH,B_from_NCSB);
        auto C = kofola::get_set_union(C_next, enrich_C);
        bool move_ncsbm_ptr = false;
        if(B.empty() && (infs_ || fins_) && rr_pointer_ == 0) {
            move_ncsbm_ptr = true;
        }

        std::vector<std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>>> new_disjs;
        for(unsigned i = 0; i < disjuncts_.size(); i++) {
            new_disjs.emplace_back(disjuncts_[i]->succs(new_runs_disj[i],symbol,info));
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
                    succ_disj.emplace_back(next_pair.first);
                    accs.emplace_back(next_pair.second);
                }
                
                
                std::shared_ptr<conj_mstate> tmp(new conj_mstate(C, S_nexts, all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_));
                
                bool move_disj_ptr = false;
                if(rr_pointer_ > 1 && accs[rr_pointer_] == 1) {
                    move_disj_ptr = true;
                }

                if(move_ncsbm_ptr || move_disj_ptr) tmp->move_rr_ptr();

                unsigned curr_acc = 0;
                if(tmp->get_rr_ptr() == 0 && (move_ncsbm_ptr || move_disj_ptr))
                    curr_acc = 1;
                
                result.emplace_back(tmp,curr_acc);
            }
            // }
            
            // guessed
            for(unsigned k = 0; k < all_guesing_safes.size(); k++) {
                auto B = kofola::get_set_union(B_from_MH, guessing_break);
                auto C = kofola::get_set_union(guessing_check, enrich_C);
                auto next_rr_ptr = rr_pointer_ + 1;
                
                B = kofola::get_set_union(C,Ms_flatt); // also redundant (since moving rr_ptr)

                for(unsigned j = 0; j < prod_of_disj_succs.size(); j++) {
                    auto tmp_disj = prod_of_disj_succs[j];
                
                    // extract disjs and acc.
                    std::vector<std::shared_ptr<disj_mstate>> succ_disj;
                    for(auto next_pair: tmp_disj) {
                        succ_disj.emplace_back(next_pair.first);
                    }

                    std::shared_ptr<conj_mstate> tmp(new conj_mstate(C, S_nexts, all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_));
                    
                    tmp->move_rr_ptr();
                    unsigned curr_acc = 0;
                    if(tmp->get_rr_ptr() == 0)
                        curr_acc = 1;
                    
                    result.emplace_back(tmp,curr_acc);
                }
                // }
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

        // If we tried to increment beyond the first element, we're done
        if (idx < 0)
            break;
    }

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
    return (
        this->check_ < other.check_ &&
        this->safes_ < other.safes_ &&
        this->m_check_ < other.m_check_ &&
        this->breakpoint_ < other.breakpoint_ &&
        this->disjuncts_ < other.disjuncts_ &&
        this->inf_colors_ < other.inf_colors_ &&
        this->fin_colors_ < other.fin_colors_ &&
        this->rr_pointer_ < other.rr_pointer_ &&
        this->active_ < other.active_
    );
}

disj_mstate::disj_mstate(spot::acc_cond cond) {
    auto all_disjs = cond.top_disjuncts();
    std::vector<spot::acc_cond> conjs;

    for(unsigned i = 0; i < all_disjs.size(); i++) {
        if(all_disjs[i].top_conjuncts().size() == 1) { // atom here
            auto inf = all_disjs[i].inf_unit();
            auto fin = all_disjs[i].fin_unit();

            if(inf.count() != 0) {
                safes_.emplace_back(); // add empty set
                inf_colors_.emplace_back(inf.min_set() - 1); // min_set returns incremented value 
            } else if(fin.count() != 0) {
                fin_colors_.emplace_back(fin.min_set() - 1); // min_set returns incremented value 
            }
        } else {
            conjs.emplace_back(all_disjs[i]);
        }
    }

    // inner disjuncts
    for(auto conj: conjs) {
        conjuncts_.emplace_back(conj_mstate(conj));
    }
}

void disj_mstate::passivate() {
    breakpoint_ = {};
    active_ = false;
    for(unsigned i = 0; i < conjuncts_.size(); i++) {
        conjuncts_[i]->passivate();
    }
}

void disj_mstate::activate() {
    rr_pointer_ = 0;
    active_ = true;
    if(infs_ || fins_) {
        breakpoint_ = check_;

    } else {
        if(conjuncts_.size() != 0) {
            conjuncts_[0]->activate();
        }
    }
}

std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> disj_mstate::succs(const std::set<unsigned>&  new_runs,
    const bdd&                 symbol,
    kofola::cmpl_info& info) {
    
    // NCSB
    std::vector<std::set<unsigned>> S_nexts;
    std::set<unsigned> S_nexts_all;

    for(unsigned i = 0; i < safes_.size(); i++) {
        auto S = safes_[i];
        if(contains_outgoing_transitions_in_scc_given_color(info.aut_, info.scc_info_, safes_[i], symbol, spot::acc_cond::mark_t(inf_colors_[i]))) {
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
    unsigned maxVal = 2 + disjuncts_.size();

    std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> result;
    std::vector<unsigned> result_acc;

    while (true) {
        std::set<unsigned> enrich_C;
        std::vector<unsigned> enrich_M;
        std::vector<std::set<unsigned>> new_runs_disj;
        for(unsigned i = 0; i < runs_cnt; i++) {
            if(run_map[i] == 0) { // this is not optimal at all - generating for infs if there are none
                if(infs_)
                    enrich_C.insert(new_runs[i]);
            } else if(run_map[i] == 1) {
                if(fins_)
                    enrich_M.emplace_back(new_runs[i]);
            } else {
                new_runs_disj[run_map[i]].insert(new_runs[i]);
            }
        }

        auto all_guesing_Ms = nondeter_scatter(m_check_nexts, enrich_M);

        auto B = kofola::get_set_union(B_from_MH,B_from_NCSB);
        auto C = kofola::get_set_union(C_next, enrich_C);
        bool move_ncsbm_ptr = false;
        if(B.empty() && (infs_ || fins_) && rr_pointer_ == 0) {
            move_ncsbm_ptr = true;
        }

        std::vector<std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>>> new_disjs;
        for(unsigned i = 0; i < disjuncts_.size(); i++) {
            new_disjs.emplace_back(disjuncts_[i]->succs(new_runs_disj[i],symbol,info));
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
                    succ_disj.emplace_back(next_pair.first);
                    accs.emplace_back(next_pair.second);
                }
                
                
                std::shared_ptr<conj_mstate> tmp(new conj_mstate(C, S_nexts, all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_));
                
                bool move_disj_ptr = false;
                if(rr_pointer_ > 1 && accs[rr_pointer_] == 1) {
                    move_disj_ptr = true;
                }

                if(move_ncsbm_ptr || move_disj_ptr) tmp->move_rr_ptr();

                unsigned curr_acc = 0;
                if(tmp->get_rr_ptr() == 0 && (move_ncsbm_ptr || move_disj_ptr))
                    curr_acc = 1;
                
                result.emplace_back(tmp,curr_acc);
            }
            // }
            
            // guessed
            for(unsigned k = 0; k < all_guesing_safes.size(); k++) {
                auto B = kofola::get_set_union(B_from_MH, guessing_break);
                auto C = kofola::get_set_union(guessing_check, enrich_C);
                auto next_rr_ptr = rr_pointer_ + 1;
                
                B = kofola::get_set_union(C,Ms_flatt); // also redundant (since moving rr_ptr)

                for(unsigned j = 0; j < prod_of_disj_succs.size(); j++) {
                    auto tmp_disj = prod_of_disj_succs[j];
                
                    // extract disjs and acc.
                    std::vector<std::shared_ptr<disj_mstate>> succ_disj;
                    for(auto next_pair: tmp_disj) {
                        succ_disj.emplace_back(next_pair.first);
                    }

                    std::shared_ptr<conj_mstate> tmp(new conj_mstate(C, S_nexts, all_guesing_Ms[i], B, succ_disj, inf_colors_, fin_colors_, infs_, fins_, rr_pointer_));
                    
                    tmp->move_rr_ptr();
                    unsigned curr_acc = 0;
                    if(tmp->get_rr_ptr() == 0)
                        curr_acc = 1;
                    
                    result.emplace_back(tmp,curr_acc);
                }
                // }
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

        // If we tried to increment beyond the first element, we're done
        if (idx < 0)
            break;
    }

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
    return (
        this->check_ < other.check_ &&
        this->safes_ < other.safes_ &&
        this->breakpoint_ < other.breakpoint_ &&
        this->conjuncts_ < other.conjuncts_ &&
        this->inf_colors_ < other.inf_colors_ &&
        this->fin_colors_ < other.fin_colors_ &&
        this->rr_pointer_ < other.rr_pointer_ &&
        this->active_ < other.active_
    );
}

// namespace { // {{{
//     /// partial macrostate for the given component
//     class mstate_tela_det : public abstract_complement_alg::mstate
//     { // {{{
//     private: // DATA MEMBERS

//     conj_mstate conj_;
//     bool active_;
//     std::set<unsigned> breakpoint_;

//     public: // METHODS

//     /// constructor
//     mstate_tela_det(spot::acc_cond cond) : conj_(cond)
//     { }

//     virtual std::string to_string() const override;
//     virtual bool is_active() const override { return this->active_; }
//     virtual bool eq(const mstate& rhs) const override;
//     virtual bool lt(const mstate& rhs) const override;
//     virtual ~mstate_tela_det() override { }

//     virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
//     virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override {  }

//     friend class kofola::complement_tela_det;
//     }; // mstate_tela_det }}}


//     std::string mstate_tela_det::to_string() const
//     {
//     std::string res = std::string("[TELA_DET(") + ((this->active_)? "A" : "T") + "): ";
//     //   res += "C=" + std::to_string(this->states_);
//     //   if (this->active_) {
//     //     res += ", B=" + std::to_string(this->breakpoint_);
//     //   }
//     //   res += "]";
//     return res;
//     }

//     bool mstate_tela_det::eq(const mstate& rhs) const
//     {
//     const mstate_tela_det* rhs_mh = dynamic_cast<const mstate_tela_det*>(&rhs);
//     assert(rhs_mh);
//     return (this->conj_ == rhs_mh->conj_);
//     }

//     bool mstate_tela_det::lt(const mstate& rhs) const
//     {
//     const mstate_tela_det* rhs_mh = dynamic_cast<const mstate_tela_det*>(&rhs);
//     assert(rhs_mh);

//     return this->conj_ < rhs_mh->conj_;
//     }

//     } // anonymous namespace }}}

//     complement_tela_det::complement_tela_det(const cmpl_info& info, unsigned part_index)
//     : abstract_complement_alg(info, part_index), ncsb_alg_(info,part_index)
//     {  }

//     mstate_set complement_tela_det::get_init()
//     { // {{{
//     std::set<unsigned> init_state;

//     unsigned orig_init = this->info_.aut_->get_init_state_number();

//     if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
//         init_state.insert(orig_init);
//     }

//     mstate_set result;
//     std::shared_ptr<mstate> ms(new mstate_mh(init_state, {}, false));
//     result.push_back(ms);

//     return result;
//     } // get_init() }}}

//     mstate_col_set complement_tela_det::get_succ_track(
//     const std::set<unsigned>&  glob_reached,
//     const mstate*              src,
//     const bdd&                 symbol)
//     { // {{{
//     DEBUG_PRINT_LN("Miyano-Hayashi successor");
//     DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
//     DEBUG_PRINT_LN("src = " + std::to_string(*src));
//     DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

//     assert(src_mh);
//     assert(!src_mh->active_);

//     std::set<unsigned> states;
//     for (unsigned st : glob_reached) {
//         if (this->info_.st_to_part_map_.at(st) == static_cast<int>(this->part_index_)) {
//         states.insert(st);
//         }
//     }

//     std::shared_ptr<mstate> ms(new mstate_mh(states, {}, false));
//     return {{ms, {}}};
//     } // get_succ_track() }}}

//     mstate_set complement_tela_det::lift_track_to_active(const mstate* src)
//     { // {{{
//     const mstate_mh* src_mh = dynamic_cast<const mstate_mh*>(src);
//     assert(src_mh);
//     assert(!src_mh->active_);

//     std::shared_ptr<mstate> ms(new mstate_mh(src_mh->states_, src_mh->states_, true));
//     return {ms};
//     } // lift_track_to_active() }}}

//     mstate_col_set complement_tela_det::get_succ_active(
//     const std::set<unsigned>&  glob_reached,
//     const mstate*              src,
//     const bdd&                 symbol,
//     bool resample)
//     {
//     const mstate_mh* src_mh = dynamic_cast<const mstate_mh*>(src);
//     assert(src_mh);
//     assert(src_mh->active_);

//     DEBUG_PRINT_LN("tracking successor of: " + std::to_string(*src_mh));
//     mstate_mh tmp(src_mh->states_, {}, false);
//     mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

//     if (track_succ.size() == 0) { return {};}
//     assert(track_succ.size() == 1);

//     const mstate_mh* track_ms = dynamic_cast<const mstate_mh*>(track_succ[0].first.get());
//     assert(track_ms);

//     DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

//     std::set<unsigned> succ_break = kofola::get_all_successors_in_scc(
//         this->info_.aut_, this->info_.scc_info_, src_mh->breakpoint_, symbol);

//     // intersect with what is really reachable (for simulation pruning)
//     succ_break = kofola::get_set_intersection(succ_break, glob_reached);

//     mstate_col_set result;
//     if (succ_break.empty() && resample) { // hit breakpoint
//         if (this->use_round_robin()) {
//         std::shared_ptr<mstate> ms(new mstate_mh(track_ms->states_, {}, false));
//         result.push_back({ms, {0}});
//         } else { // no round robin
//         std::shared_ptr<mstate> ms(new mstate_mh(track_ms->states_, track_ms->states_, true));
//         result.push_back({ms, {0}});
//         }
//     }
//     else { // no breakpoint
//         std::shared_ptr<mstate> ms(new mstate_mh(track_ms->states_, succ_break, true));
//         result.push_back({ms, {}});
//     }

//     return result;
// }

// complement_tela_det::~complement_tela_det()
// { }

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

        /// removes duplicit values (warning: can change order!)
        template <class T>
        static void remove_duplicit(T& t){ // {{{
            std::sort(t.begin(), t.end());
            t.erase(std::unique(t.begin(), t.end()), t.end());
        } // remove_duplicit }}}
#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

class disj_mstate;

class conj_mstate
{
private:
    std::set<unsigned> all_states_;
    std::set<unsigned> check_;
    std::vector<std::set<unsigned>> safes_;
    std::vector<std::set<unsigned>> m_check_;
    std::set<unsigned> breakpoint_;
    bool active_;
    std::vector<std::shared_ptr<disj_mstate>> disjuncts_;
    std::vector<unsigned> inf_colors_;
    std::vector<unsigned> fin_colors_;
    unsigned rr_pointer_ = 0;

    bool infs_ = false;
    bool fins_ = false;

public:
    conj_mstate(spot::acc_cond cond);
    conj_mstate(std::set<unsigned> check,
                std::vector<std::set<unsigned>> safes,
                std::vector<std::set<unsigned>> m_check,
                std::set<unsigned> breakpoint,
                std::vector<std::shared_ptr<disj_mstate>> disjuncts,
                std::vector<unsigned> inf_colors,
                std::vector<unsigned> fin_colors,
                bool infs,
                bool fins,
                unsigned rr_ptr,
                bool active
              ) {
                  check_ = std::move(check);
                  safes_ = std::move(safes);
                  m_check_ = std::move(m_check);
                  breakpoint_ = std::move(breakpoint);
                  active_ = active;
                  disjuncts_ = disjuncts;
                  inf_colors_ = std::move(inf_colors);
                  fin_colors_ = std::move(fin_colors);
                  infs_ = infs;
                  fins_ = fins;    
                  rr_pointer_ = rr_ptr;              
              }

    std::set<unsigned> get_all_states();

    std::shared_ptr<conj_mstate> clone() const;

    std::string to_str();

    bool get_activity() {
        return active_;
    }

    void passivate();
    void activate();
    void move_rr_ptr();
    unsigned get_rr_ptr();

    std::vector<std::pair<std::shared_ptr<conj_mstate>, unsigned>> succs(
    const std::vector<unsigned>&  new_runs,
    const bdd&                 symbol,
    const kofola::cmpl_info& info);

    bool operator==(const conj_mstate& other) const;

    bool operator<(const conj_mstate& other) const;
}; // conj_mstate

class disj_mstate
{
private:
    std::set<unsigned> check_;
    std::set<unsigned> safes_;
    std::set<unsigned> breakpoint_;
    bool active_;
    std::vector<std::shared_ptr<conj_mstate>> conjuncts_;
    std::vector<unsigned> inf_colors_;
    std::vector<unsigned> fin_colors_;
    unsigned bins;
    unsigned rr_pointer_;

    bool infs_ = false;
    bool fins_ = false;

    std::vector<spot::acc_cond::mark_t> fins_marks_set_;

public:
    disj_mstate(spot::acc_cond cond);
    disj_mstate(std::set<unsigned> check,
                std::set<unsigned> safes,
                std::set<unsigned> breakpoint,
                std::vector<std::shared_ptr<conj_mstate>> conjuncts,
                std::vector<unsigned> inf_colors,
                std::vector<unsigned> fin_colors,
                bool infs,
                bool fins,
                unsigned rr_ptr,
                bool active
              ) {
                  check_ = std::move(check);
                  safes_ = std::move(safes);
                  breakpoint_ = std::move(breakpoint);
                  active_ = active;
                  conjuncts_ = conjuncts;
                  inf_colors_ = std::move(inf_colors);
                  fin_colors_ = std::move(fin_colors);
                  infs_ = infs;
                  fins_ = fins;    
                  rr_pointer_ = rr_ptr;

                  for(unsigned i = 0; i < fin_colors_.size(); i++) {
                      fins_marks_set_.emplace_back(spot::acc_cond::mark_t{fin_colors_[i]});
                  }
              }

    std::set<unsigned> get_all_states();

    std::shared_ptr<disj_mstate> clone() const;


    std::string to_str();

    void passivate();
    void activate();

    void move_rr_ptr();
    unsigned get_rr_ptr();

    std::vector<std::pair<std::shared_ptr<disj_mstate>, unsigned>> succs(
    const std::vector<unsigned>&  new_runs,
    const bdd&                 symbol,
    const kofola::cmpl_info& info);

    bool operator==(const disj_mstate& other) const;

    bool operator<(const disj_mstate& other) const;

}; // disj_mstate

/// implementation of the Miyano & Hayashi complementation algorithm for
/// inherently weak SCCs
class complement_tela_det : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_tela_det(const cmpl_info& info, unsigned part_index);

  virtual mstate_set get_init() override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) override;

  virtual mstate_set lift_track_to_active(const mstate* src) override;

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol,
    bool resample = true) override;

  virtual bool use_round_robin() const override { return false; }

  virtual bool use_shared_breakpoint() const override { return this->info_.shared_breakpoint_; }

  virtual spot::acc_cond get_acc_cond() override
  { return spot::acc_cond(1, spot::acc_cond::inf({0})); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~complement_tela_det() override;
}; // complement_tela_det }}}
} // namespace kofola }}}


#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

class determinisation_acc_cond 
{ // {{{
private:
    spot::acc_cond::acc_code acc_code_;
    std::vector<unsigned> additional_fins_;
    unsigned disj_size_;
    
public:
    determinisation_acc_cond(const spot::acc_cond::acc_code& acc_cond, unsigned disjuncts);

    spot::acc_cond get_acc_cond() const;

    unsigned get_min_colour() const;

    spot::acc_cond::mark_t get_all_discontinuation_colours() const;

    unsigned get_fin_mark(unsigned disjunct_index) const {
        assert(disjunct_index < additional_fins_.size());
        return additional_fins_[disjunct_index];
    }

    unsigned map_colour(unsigned colour, unsigned disjunct_index) const {
        return colour + disj_size_ * disjunct_index;
     }
};

class complement_init_almost_det : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_init_almost_det(const cmpl_info& info, unsigned part_index);

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

  // no breakpoint at all
  virtual bool use_shared_breakpoint() const override { return false; }

  virtual spot::acc_cond get_acc_cond() override
  { return acc_cond_.get_acc_cond().get_acceptance().complement(); }

  virtual unsigned get_min_colour() const override { return acc_cond_.get_min_colour(); }

  static unsigned count_part_states(const cmpl_info& info, unsigned part_index);

  virtual ~complement_init_almost_det() override;

private:
    unsigned runs_bound_;
    determinisation_acc_cond acc_cond_;
}; // complement_init_almost_det }}}
} // namespace kofola }}}


#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

class complement_init_det : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_init_det(const cmpl_info& info, unsigned part_index);

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
  { return spot::acc_cond(2, spot::acc_cond::fin({1})); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~complement_init_det() override;
}; // complement_init_det }}}
} // namespace kofola }}}


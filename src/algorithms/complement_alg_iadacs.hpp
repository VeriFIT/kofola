#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

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
  { return dualized_acc_cond_; }

  virtual unsigned get_min_colour() const override { return min_colour_; }

  virtual ~complement_init_almost_det() override;

private:
    spot::acc_cond dualized_acc_cond_;
    unsigned min_colour_;
}; // complement_init_almost_det }}}
} // namespace kofola }}}


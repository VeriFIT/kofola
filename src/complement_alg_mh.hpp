// implementation of the Miyano & Hayashi complementation algorithm for
// inherently weak SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of the Miyano & Hayashi complementation algorithm for
/// inherently weak SCCs
class complement_mh : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_mh(const cmpl_info& info, unsigned part_index);

  virtual mstate_set get_init() const override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual mstate_set lift_track_to_active(const mstate* src) const override;

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual bool use_round_robin() const override { return false; }

  virtual spot::acc_cond get_acc_cond() const override
  { return spot::acc_cond(1, spot::acc_cond::inf({0})); }

  virtual ~complement_mh() override;
}; // complement_mh }}}
} // namespace kofola }}}


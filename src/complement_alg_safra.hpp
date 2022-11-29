// implementation of determinization-based complementation algorithm for
// nondeterministic accepting SCCs based on Safra's trees and producing parity
// condition in the similar way as in
// Yong Li, Andrea Turrini, Weizhi Feng,
// Moshe Y. Vardi, Lijun Zhang: Divide-and-Conquer Determinization of Büchi
// Automata Based on SCC Decomposition. CAV (2) 2022: 152-173

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of determinization-based complementation algorithm for
/// nondeterministic accepting SCCs based on Safra's trees and producing parity
/// condition
class complement_safra : public abstract_complement_alg
{ // {{{
private:// DATA MEMBERS

  // to keep track of minimum/maximum colours (mutable to be usable in const
  // methods)
  mutable int min_colour_ = INT_MAX;
  mutable int max_colour_ = -1;

public: // METHODS

  /// constructor
  complement_safra(const cmpl_info& info, unsigned part_index);

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

  virtual bool use_shared_breakpoint() const override { return false; }

  /// note: should be called only after the construction is finished (otherwise
  /// some colours might be missing)
  virtual spot::acc_cond get_acc_cond() override;

  virtual unsigned get_min_colour() const override;

  virtual ~complement_safra() override;
}; // complement_safra }}}
} // namespace kofola }}}

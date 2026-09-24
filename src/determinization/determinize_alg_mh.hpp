// implementation of the Miyano & Hayashi determinization algorithm for
// inherently weak SCCs

#pragma once

#include "abstract_determinize_alg.hpp"

namespace kofola { // {{{

/// Determinization of an inherently weak (accepting) partition via the
/// Miyano & Hayashi breakpoint construction.
///
/// A word is accepted through the partition iff some run of the input automaton
/// eventually stays inside the partition forever (all SCCs of the partition are
/// accepting and inherently weak).  This is detected by the breakpoint: the
/// breakpoint is emptied only finitely often iff some run settles in the
/// partition, hence the acceptance condition is Fin(0) - the dual of the Inf(0)
/// used by @p complement_mh.
class determinize_mh : public abstract_determinize_alg
{ // {{{
public: // METHODS

  /// constructor
  determinize_mh(const cmpl_info& info, unsigned part_index);

  virtual mstate_p get_init() override;

  virtual mstate_col get_succ(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) override;

  virtual spot::acc_cond get_acc_cond() override
  { return spot::acc_cond(1, spot::acc_cond::fin({0})); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~determinize_mh() override;
}; // determinize_mh }}}
} // namespace kofola }}}

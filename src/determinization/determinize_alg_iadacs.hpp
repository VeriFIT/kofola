// partial determinization algorithm for initial almost deterministic
// accepting components (IADACs)

#pragma once

#include "abstract_determinize_alg.hpp"
#include "../algorithms/complement_alg_iadacs.hpp"

namespace kofola { // {{{

/// Determinization of a partition made of initial almost deterministic
/// accepting components (IADACs).
///
/// This is the determinization that the IADACs complementation
/// (@p complement_init_almost_det) is based on, taken without the final
/// complementation of the acceptance condition: the macrostates, the
/// successors, the colours and the (uncomplemented) acceptance condition are
/// exactly those of complement_init_almost_det::get_succ_active().
///
/// # Macrostate
///
/// A macrostate is a vector of 'runs_bound' slots, each holding a state of the
/// partition or UINT_MAX; the live runs always occupy a prefix of the slots.
/// The position of a run in the vector is its label.  A step moves every run
/// along its (unique) successor within its SCC, drops the runs that leave the
/// SCC or merge into a run with a smaller label, appends the runs entering the
/// partition in the order of the state numbers, and finally compresses the
/// vector back to a prefix.
///
/// # Acceptance
///
/// The labels are the disjuncts of @p determinisation_acc_cond: the label i
/// owns a copy of the colours of the input automaton plus a discontinuation
/// colour, and its disjunct is the acceptance of the input automaton over that
/// copy conjoined with Fin(discontinuation colour).  In a step, the label i
/// sees the colours of the transition taken by its run as long as every label
/// up to i keeps carrying the same run; from the first label that does not, all
/// the labels are discontinued.
///
/// Unlike @p determinize_dac, which discovers the number of labels on the fly,
/// the number of slots (and thus the colour range) is fixed upfront by
/// helpers::tnba_complement::max_runs_in_partition(), exactly as in the
/// complementation.
class determinize_iadacs : public abstract_determinize_alg
{ // {{{
private: // DATA MEMBERS

  /// the number of slots of a macrostate
  unsigned runs_bound_;

  /// the acceptance condition, shared with the IADACs complementation
  determinisation_acc_cond acc_cond_;

public: // METHODS

  /// constructor
  determinize_iadacs(const cmpl_info& info, unsigned part_index, unsigned runs_bound);

  virtual mstate_p get_init() override;

  virtual mstate_col get_succ(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) override;

  virtual spot::acc_cond get_acc_cond() override
  { return acc_cond_.get_acc_cond(); }

  virtual unsigned get_min_colour() const override { return acc_cond_.get_min_colour(); }

  virtual ~determinize_iadacs() override;
}; // determinize_iadacs }}}
} // namespace kofola }}}

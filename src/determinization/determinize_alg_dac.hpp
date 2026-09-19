// partial determinization algorithm for deterministic accepting components (DACs)

#pragma once

#include <unordered_map>
#include <vector>

#include "abstract_determinize_alg.hpp"

namespace kofola { // {{{

/// Determinization of a partition made of deterministic accepting components
/// (DACs), i.e. accepting SCCs that are deterministic *inside* (jumps leaving
/// an SCC may still be nondeterministic).
///
/// This is an Emerson-Lei generalization of the DAC part of the
/// divide-and-conquer Buchi determinization of Li, Turrini, Feng, Vardi and
/// Zhang (CAV'22).
///
/// # Run labelling
///
/// Inside the partition the runs are deterministic, so they can only *merge*,
/// never branch; new runs enter from the outside.  Every run alive in the
/// partition carries a label (a natural number) maintained by two rules:
///
///   (*)  the greatest dies first - when two runs merge, only the smaller label
///        survives and the larger one is discontinued;
///   (**) the elder takes precedence - runs already inside the partition get
///        smaller labels than the runs entering in the current step (ties among
///        the entering ones are broken by the state number).
///
/// Consequently the label of an infinite run staying in the partition is
/// non-increasing, hence it stabilizes, and the stable label identifies the run.
/// So "a run stays in the partition forever" is equivalent to "its label changes
/// (or the run dies) only finitely often", which is turned into an acceptance
/// condition by giving every label its own private block of colours.
///
/// # Macrostate
///
/// By the two invariants of the labelling - the support is exactly the set of
/// reached states of the partition, and the labels form an initial segment of
/// the naturals - a macrostate is just an ordered sequence of distinct states of
/// the partition, the position in the sequence being the label.
class determinize_dac : public abstract_determinize_alg
{ // {{{
private: // DATA MEMBERS

  /// states of the partition, ordered by the (fixed) state number order;
  /// determines the tie-breaking order of simultaneously entering runs
  std::vector<unsigned> part_states_;

  /// number of labels the algorithm may use (= number of states of the
  /// partition, the maximum number of simultaneously live runs)
  unsigned num_labels_;

  /// colours occurring inside the SCCs of the partition
  spot::acc_cond::mark_t used_colours_;

  /// maps a colour of used_colours_ to its index in {0,...,num_colours_-1}
  std::unordered_map<unsigned, unsigned> colour_rank_;

  /// number of colours of the local acceptance condition
  unsigned num_colours_;

  /// the local acceptance condition, over the colours {0,...,num_colours_-1}
  spot::acc_cond::acc_code local_acc_;

  /// The largest number of simultaneously live runs seen during the
  /// construction, i.e. the number of labels that are actually used.  It is
  /// discovered on the fly and only read by get_acc_cond(), which the top-level
  /// algorithm calls once the state space is complete; get_succ() emits the
  /// death colours for all the num_labels_ labels and the ones belonging to a
  /// label that turned out to be unused are dropped afterwards (see the
  /// contract of abstract_determinize_alg::get_acc_cond()).  This matters: the
  /// condition has one disjunct of num_colours_+1 colours per label, so using
  /// |D| labels instead of the reachable maximum quickly exhausts Spot's colour
  /// budget.
  unsigned max_live_;

public: // METHODS

  /// constructor
  determinize_dac(const cmpl_info& info, unsigned part_index);

  virtual mstate_p get_init() override;

  virtual mstate_col get_succ(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) override;

  virtual spot::acc_cond get_acc_cond() override;

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~determinize_dac() override;

private: // METHODS

  /// the colour of the original automaton 'orig_colour' as seen by the label
  /// 'label'
  unsigned shift(unsigned label, unsigned colour) const
  { return label * (this->num_colours_ + 1) + colour; }

  /// the colour signalling that the run with the label 'label' was discontinued
  unsigned die(unsigned label) const
  { return this->shift(label, this->num_colours_); }

  /// the unique successor of 'state' over 'symbol' within its SCC; returns
  /// false if there is none (the run leaves the SCC and is discontinued)
  bool succ_in_scc(
    unsigned                 state,
    const bdd&               symbol,
    unsigned&                succ,
    spot::acc_cond::mark_t&  cols) const;
}; // determinize_dac }}}
} // namespace kofola }}}

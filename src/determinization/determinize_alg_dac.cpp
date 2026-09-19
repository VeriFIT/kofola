// partial determinization algorithm for deterministic accepting components (DACs)

#include "determinize_alg_dac.hpp"

#include <algorithm>

using namespace kofola;
using mstate_p = abstract_determinize_alg::mstate_p;
using mstate_col = abstract_determinize_alg::mstate_col;


namespace { // {{{

/// partial macrostate: an ordered sequence of distinct states of the partition,
/// the position in the sequence being the label of the run
class mstate_dac : public abstract_determinize_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::vector<unsigned> labels_;

public: // METHODS

  /// constructor
  explicit mstate_dac(const std::vector<unsigned>& labels) : labels_(labels) { }

  virtual std::string to_string() const override;
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_dac() override { }

  friend class kofola::determinize_dac;
}; // mstate_dac }}}


std::string mstate_dac::to_string() const
{
  return std::string("[DAC: ") + std::to_string(this->labels_) + "]";
}

bool mstate_dac::eq(const mstate& rhs) const
{
  const mstate_dac* rhs_dac = dynamic_cast<const mstate_dac*>(&rhs);
  assert(rhs_dac);
  return this->labels_ == rhs_dac->labels_;
}

bool mstate_dac::lt(const mstate& rhs) const
{
  const mstate_dac* rhs_dac = dynamic_cast<const mstate_dac*>(&rhs);
  assert(rhs_dac);
  return this->labels_ < rhs_dac->labels_;   // lexicographic
}

} // anonymous namespace }}}


determinize_dac::determinize_dac(const cmpl_info& info, unsigned part_index)
  : abstract_determinize_alg(info, part_index),
    part_states_(),
    num_labels_(0),
    used_colours_({}),
    colour_rank_(),
    num_colours_(0),
    local_acc_(spot::acc_cond::acc_code::f()),
    max_live_(0)
{ // {{{
  // the states of the partition in the fixed order used to break ties among
  // simultaneously entering runs (any fixed order works, but it has to stay the
  // same during the whole construction, otherwise macrostates lose canonicity)
  for (unsigned st = 0; st < this->info_.aut_->num_states(); ++st) {
    auto it = this->info_.st_to_part_map_.find(st);
    if (this->info_.st_to_part_map_.end() != it &&
        it->second == static_cast<int>(this->part_index_)) {
      this->part_states_.push_back(st);
    }
  }
  this->num_labels_ = this->part_states_.size();

  // Colours occurring on the transitions the algorithm follows.  Note that a
  // partition may consist of several DACs; transitions between two of them are
  // *not* followed (see succ_in_scc()), so only the colours inside the SCCs
  // count here.
  for (unsigned st : this->part_states_) {
    unsigned scc = this->info_.scc_info_.scc_of(st);
    for (const auto& t : this->info_.aut_->out(st)) {
      if (this->info_.scc_info_.scc_of(t.dst) == scc) {
        this->used_colours_ |= t.acc;
      }
    }
  }

  // The local acceptance condition: the condition of the input automaton with
  // every colour that does not occur inside the partition declared missing
  // (Inf(c) becomes false, Fin(c) becomes true).  strip() additionally
  // renumbers the remaining colours down to {0,...,num_colours_-1}, packing
  // them while preserving their order - which is exactly how colour_rank_
  // below maps the original colours.
  spot::acc_cond::mark_t unused =
    this->info_.aut_->acc().all_sets() - this->used_colours_;
  this->local_acc_ = this->info_.aut_->acc().get_acceptance().strip(unused, true);

  for (unsigned col : this->used_colours_.sets()) {
    this->colour_rank_[col] = this->num_colours_;
    ++this->num_colours_;
  }

  DEBUG_PRINT_LN("DAC determinization of partition " + std::to_string(part_index) +
                 ": states = " + std::to_string(this->part_states_) +
                 ", colours = " + std::to_string(this->used_colours_) +
                 ", local acc = " + std::to_string(this->local_acc_));
} // constructor }}}

determinize_dac::~determinize_dac()
{ }

bool determinize_dac::succ_in_scc(
  unsigned                 state,
  const bdd&               symbol,
  unsigned&                succ,
  spot::acc_cond::mark_t&  cols) const
{ // {{{
  // Only transitions that stay within the SCC of 'state' are followed.  For a
  // partition holding several DACs this means that a run moving from one of
  // them to another is treated as discontinued and immediately re-enters with a
  // fresh (larger) label.  That is sound because a run staying in the partition
  // forever eventually stays in a single SCC of it, so it moves only finitely
  // often; and it is what keeps the followed transition relation deterministic.
  unsigned scc = this->info_.scc_info_.scc_of(state);
  bool found = false;

  for (const auto& t : this->info_.aut_->out(state)) {
    if (this->info_.scc_info_.scc_of(t.dst) != scc) { continue; }
    if (!bdd_implies(symbol, t.cond)) { continue; }

    assert(!found);   // the SCCs of the partition are deterministic inside
    found = true;
    succ = t.dst;
    cols = t.acc;
  }

  return found;
} // succ_in_scc() }}}

mstate_p determinize_dac::get_init()
{ // {{{
  std::vector<unsigned> labels;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  auto it = this->info_.st_to_part_map_.find(orig_init);
  if (this->info_.st_to_part_map_.end() != it &&
      it->second == static_cast<int>(this->part_index_)) {
    labels.push_back(orig_init);
  }

  this->max_live_ = std::max<unsigned>(this->max_live_, labels.size());

  // the partition is typically entered only later, so the initial macrostate is
  // usually the empty labelling - which is a perfectly normal macrostate here
  return std::make_shared<mstate_dac>(labels);
} // get_init() }}}

mstate_col determinize_dac::get_succ(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{ // {{{
  const mstate_dac* src_dac = dynamic_cast<const mstate_dac*>(src);
  assert(src_dac);

  DEBUG_PRINT_LN("DAC determinization successor");
  DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("src = " + std::to_string(*src_dac));
  DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

  const std::vector<unsigned>& labels = src_dac->labels_;
  std::vector<unsigned> succ_labels;
  std::set<unsigned> cols;

  // 'lost' is the smallest label lost in this step (num_labels_ if none was);
  // while no label has been lost yet, a surviving run keeps its own label, so
  // the output index equals the input index
  unsigned lost = this->num_labels_;

  // Surviving runs, walked in the order of their labels.  Appending them in
  // this order realizes rule (*) - the first (hence smallest) label reaching a
  // state wins the merge - and the compression of the labels back into an
  // initial segment, both in a single pass.
  for (unsigned i = 0; i < labels.size(); ++i) {
    unsigned succ = 0;
    spot::acc_cond::mark_t succ_cols {};

    if (!this->succ_in_scc(labels[i], symbol, succ, succ_cols) ||
        succ_labels.end() != std::find(succ_labels.begin(), succ_labels.end(), succ)) {
      // the run either left the SCC or merged into a smaller label
      if (this->num_labels_ == lost) { lost = i; }
      continue;
    }

    if (this->num_labels_ == lost) { // the run survives at its own label 'i'
      for (unsigned col : (succ_cols & this->used_colours_).sets()) {
        cols.insert(this->shift(i, this->colour_rank_.at(col)));
      }
    }

    succ_labels.push_back(succ);
  }

  // Entering runs, appended after all the survivors in the fixed state order -
  // this is rule (**).  'glob_reached' is ordered by the state number, so the
  // tie-breaking order comes for free.  Testing against succ_labels is enough:
  // after the loop above it holds exactly the successors reached from within
  // the partition.
  for (unsigned st : glob_reached) {
    auto it = this->info_.st_to_part_map_.find(st);
    if (this->info_.st_to_part_map_.end() == it ||
        it->second != static_cast<int>(this->part_index_)) { continue; }
    if (succ_labels.end() != std::find(succ_labels.begin(), succ_labels.end(), st)) { continue; }

    succ_labels.push_back(st);
  }

  assert(succ_labels.size() <= this->num_labels_);
  this->max_live_ = std::max<unsigned>(this->max_live_, succ_labels.size());

  // Death colours.  A label is reported as discontinued whenever it does not
  // carry the same run before and after the step: that is the case for every
  // label from 'lost' on (those are shifted down by the compression, if they
  // survive at all) and for every label that is not live in the successor.
  //
  // [deviation] The specification only emits die(i) for lost <= i, which is
  // unsound for a condition satisfied by the empty colour set (e.g. Fin(0)):
  // once all the runs are gone, no label is ever lost again, so no die(i) is
  // emitted and the disjunct of a label that is not live would hold vacuously.
  // Note that succ_labels.size() >= lost whenever lost < num_labels_, so this
  // only strengthens the case where nothing was lost.
  unsigned dead_from = std::min<unsigned>(lost, succ_labels.size());
  for (unsigned i = dead_from; i < this->num_labels_; ++i) {
    cols.insert(this->die(i));
  }

  return {std::make_shared<mstate_dac>(succ_labels), cols};
} // get_succ() }}}

spot::acc_cond determinize_dac::get_acc_cond()
{ // {{{
  // Every label owns a private block of num_colours_+1 colours: a copy of the
  // local condition's colours plus its own 'die' colour.  The word is accepted
  // within the partition iff some label eventually stabilizes on a run that
  // never dies and that run satisfies the local condition; as the blocks are
  // pairwise disjoint, the disjuncts are checked completely independently.
  // Only the labels that have actually been live during the construction get a
  // disjunct; a label that was never used carries no run, so its disjunct could
  // only be satisfied vacuously.  The death colours emitted for those labels are
  // dropped by the top-level algorithm.
  spot::acc_cond::acc_code res = spot::acc_cond::acc_code::f();

  for (unsigned i = 0; i < this->max_live_; ++i) {
    spot::acc_cond::acc_code disjunct = this->local_acc_;
    disjunct <<= this->shift(i, 0);                    // rename j to shift(i, j)
    disjunct &= spot::acc_cond::acc_code::fin({this->die(i)});
    res |= disjunct;
  }

  return spot::acc_cond(this->max_live_ * (this->num_colours_ + 1), res);
} // get_acc_cond() }}}

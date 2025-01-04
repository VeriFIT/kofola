// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_ncsb.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // anonymous namespace {{{

/// partial macrostate for the given component
class mstate_ncsb : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::set<unsigned> check_;       // states for runs that need to be checked
  std::set<unsigned> safe_;        // safe states (cannot see accepting transition)
  std::set<unsigned> breakpoint_;
  bool active_;                    // true = active ; false = track

public: // METHODS

  /// constructor
  mstate_ncsb(
    const std::set<unsigned>&  check,
    const std::set<unsigned>&  safe,
    const std::set<unsigned>&  breakpoint,
    bool                       active
  ) : check_(check),
    safe_(safe),
    breakpoint_(breakpoint),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_ncsb() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { this->breakpoint_ = get_set_intersection(breakpoint, this->check_); }

  virtual bool subsum_less_early(const mstate& rhs) override {
    auto rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);

    auto S_subs = std::includes(this->safe_.begin(), this->safe_.end(), rhs_ncsb->safe_.begin(), rhs_ncsb->safe_.end());
    if(!S_subs) {
      return false;
    }

    std::set<unsigned> S_and_B;
    set_union(safe_.begin(), safe_.end(), breakpoint_.begin() , breakpoint_.end(), std::inserter(S_and_B, S_and_B.begin()));

    auto B_subs = std::includes(S_and_B.begin(), S_and_B.end(), rhs_ncsb->breakpoint_.begin(), rhs_ncsb->breakpoint_.end());
    
    return B_subs;
  };

  friend class kofola::complement_ncsb;
}; // mstate_ncsb }}}


/// returns true of there is at least one outgoing accepting transition from
/// a set of states over the given symbol in the SCC the source state is in
bool contains_accepting_outgoing_transitions_in_scc(
  const spot::const_twa_graph_ptr&    aut,
  const spot::scc_info&               scc_info,
  const std::set<unsigned>&           states,
  const bdd&                          symbol)
{ // {{{
  for (unsigned s : states) {
    for (const auto &t : aut->out(s)) {
      if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(symbol, t.cond)) {
        if (t.acc) { return true; }
      }
    }
  }

  return false;
} // contains_accepting_outgoing_transitions() }}}

} // anonymous namespace }}}



std::string mstate_ncsb::to_string() const
{
  std::string res = std::string("[NCSB(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->check_);
  res += ", S=" + std::to_string(this->safe_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

bool mstate_ncsb::eq(const mstate& rhs) const
{
  const mstate_ncsb* rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);
  assert(rhs_ncsb);
  return (this->active_ == rhs_ncsb->active_) &&
    (this->check_ == rhs_ncsb->check_) &&
    (this->safe_ == rhs_ncsb->safe_) &&
    (this->breakpoint_ == rhs_ncsb->breakpoint_);
}


bool mstate_ncsb::lt(const mstate& rhs) const
{ // {{{
  const mstate_ncsb* rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);
  assert(rhs_ncsb);

  if (this->active_ != rhs_ncsb->active_) { return this->active_ < rhs_ncsb->active_; }
  if (this->check_ != rhs_ncsb->check_) { return this->check_ < rhs_ncsb->check_; }
  if (this->safe_ != rhs_ncsb->safe_) { return this->safe_ < rhs_ncsb->safe_; }
  if (this->breakpoint_ != rhs_ncsb->breakpoint_) { return this->breakpoint_ < rhs_ncsb->breakpoint_; }

  return false;   // if all are equal
} // lt() }}}


complement_ncsb::complement_ncsb(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
{ }

mstate_set complement_ncsb::get_init()
{ // {{{
  DEBUG_PRINT_LN("init NCSB for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new mstate_ncsb(init_state, {}, {}, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

mstate_col_set complement_ncsb::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(!src_ncsb->active_);

  // check that safe states do not see accepting transition in the same SCC
  if (contains_accepting_outgoing_transitions_in_scc(
      this->info_.aut_,
      this->info_.scc_info_,
      src_ncsb->safe_, symbol)) {
    return {};
  }

  std::set<unsigned> succ_safe = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_ncsb->safe_, symbol);

  std::set<unsigned> succ_states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == static_cast<int>(this->part_index_)) {
      if (succ_safe.find(st) == succ_safe.end()) { // if not in safe
        succ_states.insert(st);
      }
    }
  }

  // intersect with what is really reachable (for simulation pruning)
  // TODO: make intersection with glob_reached()

  std::shared_ptr<mstate> ms(new mstate_ncsb(succ_states, succ_safe, {}, false));
  mstate_col_set result = {{ms, {}}}; return result;
} // get_succ_track() }}}

mstate_set complement_ncsb::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(!src_ncsb->active_);

  std::shared_ptr<mstate> ms(new mstate_ncsb(src_ncsb->check_, src_ncsb->safe_, src_ncsb->check_, true));
  return {ms};
} // lift_track_to_active() }}}

mstate_col_set complement_ncsb::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol,
  bool resample)
{
  DEBUG_PRINT_LN("computing successor for glob_reached = " + std::to_string(glob_reached) +
    ", " + std::to_string(*src) + " over " + std::to_string(symbol));
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(src_ncsb->active_);

  DEBUG_PRINT_LN("tracking successor of: " + std::to_string(*src_ncsb));
  mstate_ncsb tmp(src_ncsb->check_, src_ncsb->safe_, {}, false);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_ncsb* track_ms = dynamic_cast<const mstate_ncsb*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  std::set<unsigned> tmp_break = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_ncsb->breakpoint_, symbol);
  tmp_break = get_set_intersection(tmp_break, glob_reached);

  DEBUG_PRINT_LN("tmp_break = " + std::to_string(tmp_break));

  std::set<unsigned> succ_break = get_set_difference(tmp_break, track_ms->safe_);
  if (succ_break.empty() && resample) { // if we hit breakpoint
    mstate_col_set result;
    if (this->use_round_robin()) {
      std::shared_ptr<mstate> ms(new mstate_ncsb(track_ms->check_, track_ms->safe_, {}, false));
      result.push_back({ms, {0}});
    } else { // no round robing
      std::shared_ptr<mstate> ms(new mstate_ncsb(track_ms->check_, track_ms->safe_, track_ms->check_, true));
      result.push_back({ms, {0}});
    }
    return result;
  } else { // not breakpoint
    mstate_col_set result;
    std::shared_ptr<mstate> ms(new mstate_ncsb(track_ms->check_, track_ms->safe_, succ_break, true));
    DEBUG_PRINT_LN("standard successor: " + ms->to_string());
    result.push_back({ms, {}});

    // let us generate decreasing successor if the following three conditions hold:
    //   1) src_ncsb->breakpoint_ contains no accepting state
    //   2) succ_break contains no accepting state
    //   3) delta(src_ncsb->breakpoint_, symbol) contains no accepting condition

    // 1) check src_ncsb->breakpoint_ is not accepting
    if (succ_break.empty() || kofola::set_contains_accepting_state(src_ncsb->breakpoint_,
      this->info_.state_accepting_)) {
      return result;
    }

    // 2) check succ_break contains no accepting state
    if (kofola::set_contains_accepting_state(succ_break,
      this->info_.state_accepting_)) {
      return result;
    }

    // 3) delta(src_ncsb->breakpoint_, symbol) contains no accepting condition
    if (contains_accepting_outgoing_transitions_in_scc(
        this->info_.aut_,
        this->info_.scc_info_,
        src_ncsb->breakpoint_, symbol)) {
      return result;
    }

    // add the decreasing successor
    std::set<unsigned> decr_safe = get_set_union(track_ms->safe_, succ_break);
    std::set<unsigned> decr_check = get_set_difference(track_ms->check_, decr_safe);
    std::shared_ptr<mstate> decr_ms;
    if(this->use_shared_breakpoint()) {
      decr_ms = std::shared_ptr<mstate>(new mstate_ncsb(decr_check, decr_safe, std::set<unsigned>(), true));
      DEBUG_PRINT_LN("SB decreasing successor: " + decr_ms->to_string());
    }
    else {
      decr_ms = std::shared_ptr<mstate>(new mstate_ncsb(decr_check, decr_safe, decr_check, true));
      DEBUG_PRINT_LN("decreasing successor: " + decr_ms->to_string());
    }
    result.push_back({decr_ms, {0}});

    return result;
  }
}

complement_ncsb::~complement_ncsb()
{ }

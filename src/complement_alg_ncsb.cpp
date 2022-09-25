// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_ncsb.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

complement_ncsb::mstate_ncsb::mstate_ncsb(
  const std::set<unsigned>&  check,
  const std::set<unsigned>&  safe,
  const std::set<unsigned>&  breakpoint,
  bool                       active
  ) :
  check_(check),
  safe_(safe),
  breakpoint_(breakpoint),
  active_(active)
{ }

std::string complement_ncsb::mstate_ncsb::to_string() const
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

bool complement_ncsb::mstate_ncsb::eq(const mstate& rhs) const
{
  const mstate_ncsb* rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);
  assert(rhs_ncsb);
  return (this->active_ == rhs_ncsb->active_) &&
    (this->check_ == rhs_ncsb->check_) &&
    (this->safe_ == rhs_ncsb->safe_) &&
    (this->breakpoint_ == rhs_ncsb->breakpoint_);
}

bool complement_ncsb::mstate_ncsb::lt(const mstate& rhs) const
{
  assert(false);
}

complement_ncsb::mstate_ncsb::~mstate_ncsb()
{ }

complement_ncsb::complement_ncsb(const cmpl_info& info, unsigned scc_index)
  : abstract_complement_alg(info, scc_index)
{ }

mstate_set complement_ncsb::get_init() const
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.scc_info_.scc_of(orig_init) == this->scc_index_) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new mstate_ncsb(init_state, {}, {}, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

mstate_col_set complement_ncsb::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(!src_ncsb->active_);

  // check that safe states do not see accepting transition
  std::set<unsigned> succ_safe;
  for (unsigned s : src_ncsb->safe_) {
    for (const auto &t : this->info_.aut_->out(s)) {
      if (t.acc) { // TODO: this says *any* colour is set
        return {};
      }
      if (bdd_implies(symbol, t.cond)) { succ_safe.insert(t.dst); }
    }
  }

  std::set<unsigned> succ_states;
  for (unsigned st : glob_reached) {
    if (this->info_.scc_info_.scc_of(st) == this->scc_index_) {
      if (succ_safe.find(st) == succ_safe.end()) { // if not in safe
        succ_states.insert(st);
      }
    }
  }

  std::shared_ptr<mstate> ms(new mstate_ncsb(succ_states, succ_safe, {}, false));   // FIXME: activity
  mstate_col_set result = {{ms, {}}}; return result;
} // get_succ_track() }}}

mstate_set complement_ncsb::lift_track_to_active(const mstate* src) const
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
  const bdd&                 symbol) const
{
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(src_ncsb->active_);

  mstate_ncsb tmp(src_ncsb->check_, src_ncsb->safe_, {}, false);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);
  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_ncsb* track_ms = dynamic_cast<const mstate_ncsb*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  std::set<unsigned> tmp_break = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, this->scc_index_, src_ncsb->breakpoint_, symbol);

  std::set<unsigned> succ_break = get_set_difference(tmp_break, track_ms->safe_);
  if (succ_break.empty()) { // if we hit breakpoint
    std::shared_ptr<mstate> ms(new mstate_ncsb(track_ms->check_, track_ms->safe_, succ_break, false));
    return {{ms, {0}}};
  } else {
    // FIXME: add decreasing transitions
    assert(false);
  }
}

complement_ncsb::~complement_ncsb()
{ }

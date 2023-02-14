// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_ncsb_delay.hpp"

#include <stack>

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

complement_ncsb_delay::complement_ncsb_delay(const cmpl_info& info, unsigned part_index)
  : complement_ncsb(info, part_index)
{ }

mstate_set complement_ncsb_delay::get_init()
{ // {{{
  DEBUG_PRINT_LN("init NCSB for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new mstate_ncsb(init_state, {}, {}, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

mstate_col_set complement_ncsb_delay::get_succ_track(
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
      this->info_.st_to_part_map_,
      this->info_.scc_info_,
      src_ncsb->safe_, symbol)) {
    return {};
  }

  std::set<unsigned> succ_safe = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_ncsb->safe_, symbol);

  std::set<unsigned> succ_states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == this->part_index_) {
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

mstate_set complement_ncsb_delay::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_ncsb* src_ncsb = dynamic_cast<const mstate_ncsb*>(src);
  assert(src_ncsb);
  assert(!src_ncsb->active_);

  std::shared_ptr<mstate> ms(new mstate_ncsb(src_ncsb->check_, src_ncsb->safe_, src_ncsb->check_, true));
  return {ms};
} // lift_track_to_active() }}}

mstate_col_set complement_ncsb_delay::get_succ_active(
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
      auto new_state = new mstate_ncsb(track_ms->check_, track_ms->safe_, {}, false);
      std::shared_ptr<mstate> ms(new_state);
      result.push_back({ms, {0}});
      active_mstates_.insert(new_state);
    } else { // no round robing
      auto new_state = new mstate_ncsb(track_ms->check_, track_ms->safe_, track_ms->check_, true);
      std::shared_ptr<mstate> ms(new_state);
      result.push_back({ms, {0}});
      active_mstates_.insert(new_state);
      assert(new_state->active_);
      successors_[src_ncsb].insert(new_state);
      succ_ncsb_[*src_ncsb].insert(*new_state);
    }
    return result;
  } else { // not breakpoint
    mstate_col_set result;
    auto new_state = new mstate_ncsb(track_ms->check_, track_ms->safe_, succ_break, true);
    std::shared_ptr<mstate> ms(new_state);
    DEBUG_PRINT_LN("standard successor: " + ms->to_string());
    result.push_back({ms, {}});

    // let us generate decreasing successor if the following three conditions hold:
    //   1) src_ncsb->breakpoint_ contains no accepting state
    //   2) succ_break contains no accepting state
    //   3) delta(src_ncsb->breakpoint_, symbol) contains no accepting condition

    // 1) check src_ncsb->breakpoint_ is not accepting
    if (succ_break.empty() || kofola::set_contains_accepting_state(src_ncsb->breakpoint_,
      this->info_.state_accepting_)) {
      active_mstates_.insert(new_state);
      successors_[src_ncsb].insert(new_state);
      succ_ncsb_[*src_ncsb].insert(*new_state);
      return result;
    }

    // 2) check succ_break contains no accepting state
    if (kofola::set_contains_accepting_state(succ_break,
      this->info_.state_accepting_)) {
      active_mstates_.insert(new_state);
      succ_ncsb_[*src_ncsb].insert(*new_state);
      return result;
    }

    // 3) delta(src_ncsb->breakpoint_, symbol) contains no accepting condition
    if (contains_accepting_outgoing_transitions_in_scc(
        this->info_.aut_,
        this->info_.st_to_part_map_,
        this->info_.scc_info_,
        src_ncsb->breakpoint_, symbol)) {
      active_mstates_.insert(new_state);
      succ_ncsb_[*src_ncsb].insert(*new_state);
      return result;
    }

    succ_ncsb_[*src_ncsb].insert(*new_state);

    // add the decreasing successor
    if (closes_a_cycle(src_ncsb, new_state)) {
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
      auto decr = new mstate_ncsb(decr_check, decr_safe, decr_check, true);
      succ_ncsb_[*src_ncsb].insert(*decr);
    }

    return result;
  }
}

bool complement_ncsb_delay::closes_a_cycle(const mstate_ncsb *src, const mstate_ncsb *dst)
{
    assert(src->active_);
    assert(dst->active_);
    //if (successors_.find(dst) == successors_.end())
    //    return false;

    // is there a path from dst to src?
    std::stack<mstate_ncsb> stack;
    auto comparator = [](const mstate_ncsb *const &a,
                         const mstate_ncsb *const &b)
    { return a->lt(*b); };
    // std::set<const mstate_ncsb *, decltype(comparator)> visited(comparator);
    std::set<mstate_ncsb> visited;

    stack.push(*dst);

    while (not stack.empty())
    {
        mstate_ncsb current = stack.top();
        stack.pop();
        // assert(current->active_);

        visited.insert(current);

        for (mstate_ncsb succ : succ_ncsb_[current]/*successors_[current]*/)
        {
            if (succ.eq(*src) /**succ == *src*/)
            {
                return true; 
            }

            if (visited.find(succ) == visited.end())
            {
                //std::cerr << current->to_string() << std::endl;
                //std::cerr << succ->to_string() << std::endl;
                stack.push(succ);
            }
        }
    }

    return false;
}

complement_ncsb_delay::~complement_ncsb_delay()
{ }

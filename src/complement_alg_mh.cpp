// implementation of the Miyano & Hayashi complementation algorithm for
// inherently weak SCCs

#include "complement_alg_mh.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

complement_mh::mstate_mh::mstate_mh(
  const std::set<unsigned>&  states,
  const std::set<unsigned>&  breakpoint
  ) :
  states_(states),
  breakpoint_(breakpoint)
{ }

std::string complement_mh::mstate_mh::to_string() const
{
  std::string res = "[MH: ";
  res += "C=" + std::to_string(this->states_);
  res += ", B=" + std::to_string(this->breakpoint_);
  res += "]";
  return res;
}

bool complement_mh::mstate_mh::eq(const mstate& rhs) const
{
  const mstate_mh* rhs_mh = dynamic_cast<const mstate_mh*>(&rhs);
  assert(rhs_mh);
  return (this->states_ == rhs_mh->states_) &&
    (this->breakpoint_ == rhs_mh->breakpoint_);
}

bool complement_mh::mstate_mh::lt(const mstate& rhs) const
{
  const mstate_mh* rhs_mh = dynamic_cast<const mstate_mh*>(&rhs);
  assert(rhs_mh);

  if (this->states_ != rhs_mh->states_) { return this->states_ < rhs_mh->states_; }
  if (this->breakpoint_ != rhs_mh->breakpoint_) { return this->breakpoint_ < rhs_mh->breakpoint_; }

  return false;   // if all are equal
}

complement_mh::mstate_mh::~mstate_mh()
{ }

complement_mh::complement_mh(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
{ }

mstate_set complement_mh::get_init() const
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new mstate_mh(init_state, {}));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

mstate_col_set complement_mh::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  const mstate_mh* src_mh = dynamic_cast<const mstate_mh*>(src);
  assert(src_mh);

  std::set<unsigned> states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == this->part_index_) {
      states.insert(st);
    }
  }

  std::set<unsigned> succ_break = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_mh->breakpoint_, symbol);

  // intersect with what is really reachable (for simulation pruning)
  succ_break = kofola::get_set_intersection(succ_break, glob_reached);

  mstate_col_set result;
  if (succ_break.empty()) { // hit breakpoint
    std::shared_ptr<mstate> ms(new mstate_mh(states, states));
    result.push_back({ms, {0}});
  }
  else { // no breakpoint
    std::shared_ptr<mstate> ms(new mstate_mh(states, succ_break));
    result.push_back({ms, {}});
  }

  return result;
} // get_succ_track() }}}

mstate_set complement_mh::lift_track_to_active(const mstate* src) const
{
  assert(false);
}

mstate_col_set complement_mh::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}

complement_mh::~complement_mh()
{ }

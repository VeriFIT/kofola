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
  std::string res = std::string("[MH(") + ((this->active_)? "A" : "T") + "): ";
  res += "S=" + std::to_string(this->states_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

bool complement_mh::mstate_mh::eq(const mstate& rhs) const
{
  const mstate_mh* rhs_mh = dynamic_cast<const mstate_mh*>(&rhs);
  assert(rhs_mh);
  return (this->active_ == rhs_mh->active_) &&
    (this->states_ == rhs_mh->states_) &&
    (this->breakpoint_ == rhs_mh->breakpoint_);
}

bool complement_mh::mstate_mh::lt(const mstate& rhs) const
{
  assert(false);
}

complement_mh::mstate_mh::~mstate_mh()
{ }

complement_mh::complement_mh(const cmpl_info& info, unsigned scc_index)
  : abstract_complement_alg(info, scc_index)
{ }

mstate_set complement_mh::get_init() const
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.scc_info_.scc_of(orig_init) == this->scc_index_) {
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
  assert(!src_mh->active_);

  std::set<unsigned> states;
  for (unsigned st : glob_reached) {
    if (this->info_.scc_info_.scc_of(st) == this->scc_index_) {
      states.insert(st);
    }
  }

  std::shared_ptr<mstate> ms(new mstate_mh(states, {}));
  mstate_col_set result = {{ms, {}}};
  return result;
} // get_succ_track() }}}

mstate_col_set complement_mh::get_succ_track_to_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
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

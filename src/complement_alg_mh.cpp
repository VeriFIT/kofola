// implementation of the Miyano & Hayashi complementation algorithm for
// inherently weak SCCs

#include "complement_alg_mh.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;


complement_mh::mstate_mh::mstate_mh(
  const std::set<unsigned>&  states,
  const std::set<unsigned>&  breakpoint
  ) :
  states_(states),
  breakpoint_(breakpoint)
{ }

std::string complement_mh::mstate_mh::to_string() const
{
  return "TODO TODO TODO";
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

  std::shared_ptr<mstate> ms(new mstate_mh(init_state, init_state));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

mstate_set complement_mh::get_succ_track(const mstate& src, const bdd& symbol) const
{
  assert(false);
}

mstate_set complement_mh::get_succ_track_to_active(const mstate& src, const bdd& symbol) const
{
  assert(false);
}

mstate_set complement_mh::get_succ_active(const mstate& src, const bdd& symbol) const
{
  assert(false);
}

complement_mh::~complement_mh()
{ }

// implementation of NCSB-based complementation algorithm for deterministic SCCs

#include "complement_alg_ncsb.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

complement_ncsb::mstate_ncsb::mstate_ncsb(
  const std::set<unsigned>&  states,
  const std::set<unsigned>&  safe,
  const std::set<unsigned>&  breakpoint
  ) :
  states_(states),
  safe_(safe),
  breakpoint_(breakpoint)
{ }

std::string complement_ncsb::mstate_ncsb::to_string() const
{
  std::string res = std::string("[NCSB(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->states_);
  res += ", S=" + std::to_string(this->safe_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

bool complement_ncsb::mstate_ncsb::eq(const mstate& rhs) const
{
  assert(false);
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

  std::shared_ptr<mstate> ms(new mstate_ncsb(init_state, {}, {}));
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

  std::shared_ptr<mstate> ms(new mstate_ncsb(succ_states, succ_safe, {}));
  mstate_col_set result = {{ms, {}}};
  return result;
} // get_succ_track() }}}

mstate_col_set complement_ncsb::get_succ_track_to_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}

mstate_col_set complement_ncsb::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{

  // std::set<unsigned> tmp_break = kofola::get_all_successors_in_scc(
  //   this->info_.aut_, this->info_.scc_info_, this->scc_index_, src_ncsb->breakpoint_, symbol);
  //
  // std::set<unsigned> succ_break;
  // std::set_difference(tmp_break.begin(), tmp_break.end(),
  //   succ_safe.begin(), succ_safe.end(),
  //   std::inserter(succ_break, succ_break.end()));
  assert(false);
}

complement_ncsb::~complement_ncsb()
{ }

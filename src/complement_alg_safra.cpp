// implementation of determinization-based complementation algorithm for
// nondeterministic accepting SCCs based on Safra's trees and producing parity
// condition in the similar way as in
// Yong Li, Andrea Turrini, Weizhi Feng,
// Moshe Y. Vardi, Lijun Zhang: Divide-and-Conquer Determinization of Büchi
// Automata Based on SCC Decomposition. CAV (2) 2022: 152-173

#include "complement_alg_safra.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


complement_safra::mstate_safra::mstate_safra()
{ }


std::string complement_safra::mstate_safra::to_string() const
{ // {{{
  std::string res = std::string("[PAR(A):");

  res += "]";
  return res;
} // to_string() }}}


bool complement_safra::mstate_safra::lt(const mstate& rhs) const
{ // {{{
  assert(false);
} // lt() }}}


bool complement_safra::mstate_safra::eq(const mstate& rhs) const
{ // {{{
  assert(false);
} // eq() }}}


complement_safra::mstate_safra::~mstate_safra()
{ }


complement_safra::complement_safra(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index)
{ }


mstate_set complement_safra::get_init() const
{ // {{{
  assert(false);
} // get_init() }}}


mstate_col_set complement_safra::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  assert(false);
} // get_succ_track() }}}


mstate_set complement_safra::lift_track_to_active(const mstate* src) const
{ // {{{
  assert(false);
} // lift_track_to_active() }}}


mstate_col_set complement_safra::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  assert(false);
} // get_succ_active() }}}


spot::acc_cond complement_safra::get_acc_cond() const
{ // {{{
  assert(false);
} // get_acc_cond() }}}


complement_safra::~complement_safra()
{ }

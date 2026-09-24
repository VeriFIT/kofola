// implementation of the Miyano & Hayashi determinization algorithm for
// inherently weak SCCs

#include "determinize_alg_mh.hpp"

using namespace kofola;
using mstate_p = abstract_determinize_alg::mstate_p;
using mstate_col = abstract_determinize_alg::mstate_col;


namespace { // {{{

/// partial macrostate for the given partition
class mstate_det_mh : public abstract_determinize_alg::mstate
{ // {{{
private: // DATA MEMBERS

  /// states of the partition reached so far
  std::set<unsigned> states_;
  /// breakpoint: the runs that are currently checked for staying in the partition
  std::set<unsigned> breakpoint_;

public: // METHODS

  /// constructor
  mstate_det_mh(
    const std::set<unsigned>&  states,
    const std::set<unsigned>&  breakpoint
  ) : states_(states),
    breakpoint_(breakpoint)
  { }

  virtual std::string to_string() const override;
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_det_mh() override { }

  friend class kofola::determinize_mh;
}; // mstate_det_mh }}}


std::string mstate_det_mh::to_string() const
{
  std::string res = std::string("[MH-det: ");
  res += "C=" + std::to_string(this->states_);
  res += ", B=" + std::to_string(this->breakpoint_);
  res += "]";
  return res;
}

bool mstate_det_mh::eq(const mstate& rhs) const
{
  const mstate_det_mh* rhs_mh = dynamic_cast<const mstate_det_mh*>(&rhs);
  assert(rhs_mh);
  return (this->states_ == rhs_mh->states_) &&
    (this->breakpoint_ == rhs_mh->breakpoint_);
}

bool mstate_det_mh::lt(const mstate& rhs) const
{
  const mstate_det_mh* rhs_mh = dynamic_cast<const mstate_det_mh*>(&rhs);
  assert(rhs_mh);

  if (this->states_ != rhs_mh->states_) { return this->states_ < rhs_mh->states_; }
  if (this->breakpoint_ != rhs_mh->breakpoint_) { return this->breakpoint_ < rhs_mh->breakpoint_; }

  return false;   // if all are equal
}

} // anonymous namespace }}}

determinize_mh::determinize_mh(const cmpl_info& info, unsigned part_index)
  : abstract_determinize_alg(info, part_index)
{ }

mstate_p determinize_mh::get_init()
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();

  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state.insert(orig_init);
  }

  // the breakpoint initially checks all the runs in the partition
  return std::make_shared<mstate_det_mh>(init_state, init_state);
} // get_init() }}}

mstate_col determinize_mh::get_succ(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{ // {{{
  const mstate_det_mh* src_mh = dynamic_cast<const mstate_det_mh*>(src);
  assert(src_mh);

  DEBUG_PRINT_LN("Miyano-Hayashi determinization successor");
  DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("src = " + std::to_string(*src_mh));
  DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

  // all states of the partition reached over 'symbol' (runs may also enter the
  // partition from the outside)
  std::set<unsigned> succ_states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == static_cast<int>(this->part_index_)) {
      succ_states.insert(st);
    }
  }

  // the checked runs need to stay within their SCC; a run that leaves it drops
  // out of the breakpoint
  std::set<unsigned> succ_break = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_mh->breakpoint_, symbol);
  succ_break = kofola::get_set_intersection(succ_break, glob_reached);

  if (succ_break.empty()) { // hit breakpoint: no run settled, restart the check
    return {std::make_shared<mstate_det_mh>(succ_states, succ_states), {0}};
  }

  return {std::make_shared<mstate_det_mh>(succ_states, succ_break), {}};
} // get_succ() }}}

determinize_mh::~determinize_mh()
{ }

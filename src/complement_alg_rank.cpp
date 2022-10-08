// implementation of rank-based complementation

#include "complement_alg_rank.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

namespace { // {{{

/// partial macrostate for the given component
class mstate_rank : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  bool active_;
  std::set<unsigned> states_;
  std::set<unsigned> breakpoint_;

public: // METHODS

  /// constructor
  mstate_rank(
    const std::set<unsigned>&  states,
    const std::set<unsigned>&  breakpoint,
    bool                       active
  ) : states_(states),
    breakpoint_(breakpoint),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_rank() override { }

  friend class kofola::complement_rank;
}; // mstate_rank }}}

std::string mstate_rank::to_string() const
{
  assert(false);
}

bool mstate_rank::eq(const mstate& rhs) const
{
  assert(false);
}

bool mstate_rank::lt(const mstate& rhs) const
{
  assert(false);
}


} // anonymous namespace }}}


complement_rank::complement_rank(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index)
{ }


mstate_set complement_rank::get_init() const
{
  assert(false);
}

mstate_col_set complement_rank::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}

mstate_set complement_rank::lift_track_to_active(const mstate* src) const
{
  assert(false);
}

mstate_col_set complement_rank::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}


complement_rank::~complement_rank()
{ }

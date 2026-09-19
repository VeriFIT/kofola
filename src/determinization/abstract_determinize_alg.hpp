// abstract class for partial determinization algorithm and partial macrostate

#pragma once

#include <memory>
#include <ostream>
#include <set>
#include <string>

#include "../algorithms/abstract_complement_alg.hpp"

// SPOT
#include <spot/misc/bddlt.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/sccinfo.hh>


namespace kofola { // {{{

/// abstract class for partial determinization algorithms
///
/// The interface is deliberately much simpler than the one of
/// @p abstract_complement_alg: every partial algorithm is itself deterministic,
/// so there is exactly one initial partial macrostate and exactly one successor
/// for a given macrostate and symbol.  In particular, there is no breakpoint
/// handling and no tracking/active (round robin) distinction.
class abstract_determinize_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given partition
  class mstate
  { // {{{
  public: // METHODS

    /// returns string representation of the partial macrostate
    virtual std::string to_string() const = 0;

    /// equality test
    virtual bool eq(const mstate& rhs) const = 0;

    /// less-than relation
    virtual bool lt(const mstate& rhs) const = 0;

    /// virtual destructor (to allow deletion via pointer)
    virtual ~mstate() { }
  }; // mstate }}}

  /// pointer to a partial macrostate
  using mstate_p = std::shared_ptr<mstate>;
  /// partial macrostate together with the set of colours on the edge
  using mstate_col = std::pair<mstate_p, std::set<unsigned>>;

protected: // DATA MEMBERS

  /// information about the input automaton (shared with complementation)
  const cmpl_info& info_;

  /// index of the partition
  unsigned part_index_;

public: // METHODS

  /// constructor
  abstract_determinize_alg(const cmpl_info& info, unsigned part_index) :
    info_(info),
    part_index_(part_index)
  { }

  /// returns the (unique) initial partial macrostate
  virtual mstate_p get_init() = 0;

  /// returns the (unique) successor of 'src' over 'symbol'
  virtual mstate_col get_succ(
    const std::set<unsigned>&  glob_reached,  // all states reached over symbol
    const mstate*              src,           // partial macrostate
    const bdd&                 symbol) = 0;   // symbol

  /// returns the acceptance condition of the partition (over local colours)
  virtual spot::acc_cond get_acc_cond() = 0;

  /// returns the minimum colour used - allows colour reshuffle for algorithms
  /// (such as Safra-based ones) that only discover their colour range during
  /// the construction
  virtual unsigned get_min_colour() const { return 0; }

  /// virtual destructor (to allow deletion via pointer)
  virtual ~abstract_determinize_alg() { }
}; // abstract_determinize_alg }}}

/// output stream conversion
std::ostream& operator<<(std::ostream& os, const abstract_determinize_alg::mstate& ms);

/// equality operator
bool operator==(
  const abstract_determinize_alg::mstate& lhs,
  const abstract_determinize_alg::mstate& rhs);

/// disequality operator
bool operator!=(
  const abstract_determinize_alg::mstate& lhs,
  const abstract_determinize_alg::mstate& rhs);

/// ordering relation
bool operator<(
  const abstract_determinize_alg::mstate& lhs,
  const abstract_determinize_alg::mstate& rhs);

} // namespace kofola }}}

// abstract class for partial complementation algorithm and partial macrostate

#pragma once

#include <string>
#include <vector>

#include "kofola.hpp"

// SPOT
#include <spot/misc/bddlt.hh>
#include <spot/twa/twagraph.hh>
#include <spot/twaalgos/sccinfo.hh>


namespace kofola { // {{{

/// common information about automaton etc. for complementation
struct cmpl_info
{ // {{{
  /// automaton
  const spot::const_twa_graph_ptr aut_;

  /// information about SCCs
  const spot::scc_info& scc_info_;

  /// direct simulation
  const Simulation& dir_sim_;

  /// types of SCCs
  const std::string scc_types_;

  /// accepting states
  const std::vector<bool>& state_accepting_;

  /// constructor
  cmpl_info(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const Simulation&                 dir_sim,
    const std::string&                scc_types,
    const std::vector<bool>&          state_accepting
    ) :
    aut_(aut),
    scc_info_(scc_info),
    dir_sim_(dir_sim),
    scc_types_(scc_types),
    state_accepting_(state_accepting)
  { }
}; // struct cmpl_info }}}


/// abstract class for partial complementation algorithms
class abstract_complement_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given component
  class mstate
  { // {{{
  public: // METHODS

    /// is the macrostate active?
    virtual bool is_active() const = 0;

    /// returns string representation of the partial macrostate
    virtual std::string to_string() const = 0;

    /// equality test
    virtual bool eq(const mstate& rhs) const = 0;

    /// less-than relation
    virtual bool lt(const mstate& rhs) const = 0;

    /// virtual destructor (to allow deletion via pointer)
    virtual ~mstate() { }
  }; // mstate }}}

  /// set of partial macrostates together with sets of colours on the edge
  using mstate_set = std::vector<std::shared_ptr<mstate>>;
  using mstate_col = std::pair<std::shared_ptr<mstate>, std::set<unsigned>>;
  using mstate_col_set = std::vector<mstate_col>;

protected: // DATA MEMBERS

  /// information for complementation
  const cmpl_info& info_;

  /// index of the component
  unsigned scc_index_;

public: // METHODS

  /// constructor
  abstract_complement_alg(const cmpl_info& info, unsigned scc_index) :
    info_(info),
    scc_index_(scc_index)
  { }

  /// returns initial partial macrostate
  virtual mstate_set get_init() const = 0;

  /// tracking successors
  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,      // all states reached over symbol
    const mstate*              src,               // partial macrostate
    const bdd&                 symbol) const = 0; // symbol

  /// lifts tracking state to active state
  virtual mstate_set lift_track_to_active(const mstate* src) const = 0;

  /// active successors
  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,      // all states reached over symbol
    const mstate*              src,               // partial macrostate
    const bdd&                 symbol) const = 0; // symbol

  /// determines whether the algorithm should be use in round-robin scheme;
  /// in particular:
  ///   true: the algorithm uses get_succ_track(), get_succ_track_to_active(),
  ///         get_succ_active()
  ///   false: the algorithm only uses get_succ_track()
  virtual bool use_round_robin() const = 0;

  /// virtual destructor (to allow deletion via pointer)
  virtual ~abstract_complement_alg() { }
}; // abstract_complement_alg }}}

/// output stream conversion
std::ostream& operator<<(std::ostream& os, const abstract_complement_alg::mstate& ms);

/// equality operator
bool operator==(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs);

/// disequality operator
bool operator!=(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs);

/// ordering relation
bool operator<(
  const abstract_complement_alg::mstate& lhs,
  const abstract_complement_alg::mstate& rhs);

} // namespace kofola }}}

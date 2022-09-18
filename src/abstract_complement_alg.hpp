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

  /// constructor
  cmpl_info(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const Simulation&                 dir_sim
    ) :
    aut_(aut),
    scc_info_(scc_info),
    dir_sim_(dir_sim)
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

    /// returns string representation of the partial macrostate
    virtual std::string to_string() const = 0;

    /// virtual destructor (to allow deletion via pointer)
    virtual ~mstate() { }
  }; // mstate }}}

  /// set of partial macrostates
  using mstate_set = std::vector<std::shared_ptr<mstate>>;

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
  virtual mstate_set get_succ_track(const mstate& src, const bdd& symbol) const = 0;

  /// tracking to active successors
  virtual mstate_set get_succ_track_to_active(const mstate& src, const bdd& symbol) const = 0;

  /// active successors
  virtual mstate_set get_succ_active(const mstate& src, const bdd& symbol) const = 0;

  /// virtual destructor (to allow deletion via pointer)
  virtual ~abstract_complement_alg() { }
}; // abstract_complement_alg }}}
} // namespace kofola }}}

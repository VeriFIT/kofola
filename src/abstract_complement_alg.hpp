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

  /// number of partitions
  const size_t num_partitions_;

  /// types of partitions
  const PartitionToTypeMap& part_to_type_map_;

  /// a map of states to partitions
  const StateToPartitionMap& st_to_part_map_;

  /// vector of states reachable from given state
  const ReachableVector& reachable_vector_;

  /// map of partition to the SCCs it contains
  const PartitionToSCCMap& part_to_scc_map_;

  /// maps SCCs to sets of their predecessors
  const SCCToSCCSetMap& scc_to_pred_sccs_map_;

  /// information about SCCs
  const spot::scc_info& scc_info_;

  /// direct simulation
  const Simulation& dir_sim_;

  /// accepting states
  const std::vector<bool>& state_accepting_;

  /// use shared breakpoint
  const bool shared_breakpoint_;

  /// constructor
  cmpl_info(
    const spot::const_twa_graph_ptr&  aut,
    size_t                            num_partitions,
    const PartitionToTypeMap&         part_to_type_map,
    const StateToPartitionMap&        st_to_part_map,
    const ReachableVector&            reachable_vector,
    const PartitionToSCCMap&          part_to_scc_map,
    const SCCToSCCSetMap&             scc_to_pred_sccs_map,
    const spot::scc_info&             scc_info,
    const Simulation&                 dir_sim,
    const std::vector<bool>&          state_accepting,
    const bool&                       shared_breakpoint) :
    aut_(aut),
    num_partitions_(num_partitions),
    part_to_type_map_(part_to_type_map),
    st_to_part_map_(st_to_part_map),
    reachable_vector_(reachable_vector),
    part_to_scc_map_(part_to_scc_map),
    scc_to_pred_sccs_map_(scc_to_pred_sccs_map),
    scc_info_(scc_info),
    dir_sim_(dir_sim),
    state_accepting_(state_accepting),
    shared_breakpoint_(shared_breakpoint)
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

    /// get breakpoint
    virtual const std::set<unsigned>& get_breakpoint() const = 0;

    /// set breakpoint
    virtual void set_breakpoint(const std::set<unsigned>& breakpoint) = 0;

    /// clear breakpoint
    virtual void clear_breakpoint() { this->set_breakpoint(std::set<unsigned>()); };

    virtual bool subsum_less_early(const mstate& rhs, const std::set<unsigned>&  glob_reached) { return this->eq(rhs); };

    virtual bool subsum_less_early_plus(const mstate& rhs, const std::set<unsigned>&  glob_reached) { return this->eq(rhs); };

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

  /// index of the partition
  unsigned part_index_;

public: // METHODS

  /// constructor
  abstract_complement_alg(const cmpl_info& info, unsigned part_index) :
    info_(info),
    part_index_(part_index)
  { }

  /// returns initial partial macrostate
  virtual mstate_set get_init() = 0;

  /// tracking successors
  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,  // all states reached over symbol
    const mstate*              src,           // partial macrostate
    const bdd&                 symbol) = 0;   // symbol

  /// lifts tracking state to active state
  virtual mstate_set lift_track_to_active(const mstate* src) = 0;

  /// active successors
  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,  // all states reached over symbol
    const mstate*              src,           // partial macrostate
    const bdd&                 symbol,        // symbol
    bool resample = true) = 0;                // resample breakpoint

  /// determines whether the algorithm should be use in round-robin scheme;
  /// in particular:
  ///   true: the algorithm uses get_succ_track(), get_succ_track_to_active(),
  ///         get_succ_active()
  ///   false: the algorithm only uses get_succ_active()
  virtual bool use_round_robin() const = 0;

  /// returns the acceptance condition
  virtual spot::acc_cond get_acc_cond() = 0;

  /// determines whether the algorithm should use shared breakpoint
  virtual bool use_shared_breakpoint() const = 0;

  /// returns the minimum colour used - HACK to allow colour reshuffle for Safra-based algorithm
  virtual unsigned get_min_colour() const = 0;

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

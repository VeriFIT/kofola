// abstract class for partial complementation algorithm and partial macrostate

#pragma once

#include <string>
#include <vector>

namespace kofola { // {{{

class abstract_complement_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given component
  class mstate
  { // {{{
  public:

    /// returns string representation of the partial macrostate
    virtual std::string to_string() const = 0;

    /// pure virtual destructor (to allow deletion via pointer)
    virtual ~mstate() = 0;
  }; // mstate }}}

  /// set of partial macrostates
  using mstate_set = std::vector<const std::unique_ptr<mstate>>;


public: // METHODS

  /// returns initial partial macrostate
  virtual mstate_set get_init() const = 0;

  /// tracking successors
  virtual mstate_set get_succ_track(const mstate& src, const bdd& symbol) const = 0;

  /// tracking to active successors
  virtual mstate_set get_succ_track_to_active(const mstate& src, const bdd& symbol) const = 0;

  /// active successors
  virtual mstate_set get_succ_active(const mstate& src, const bdd& symbol) const = 0;

  /// pure virtual destructor (to allow deletion via pointer)
  virtual ~abstract_complement_alg() = 0;
}; // abstract_complement_alg }}}

} // namespace kofola }}}

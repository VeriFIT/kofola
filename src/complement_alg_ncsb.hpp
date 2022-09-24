// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of NCSB-based complementation algorithm for deterministic SCCs
class complement_ncsb : public abstract_complement_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given component
  class mstate_ncsb : public abstract_complement_alg::mstate
  { // {{{
  private: // DATA MEMBERS

    bool active_ = false;            // true = active ; false = track
    std::set<unsigned> states_;
    std::set<unsigned> safe_;        // safe states (cannot see accepting transition)
    std::set<unsigned> breakpoint_;

  public: // METHODS

    /// constructor
    mstate_ncsb(
      const std::set<unsigned>&  states,
      const std::set<unsigned>&  safe,
      const std::set<unsigned>&  breakpoint);

    virtual std::string to_string() const override;
    virtual bool eq(const mstate& rhs) const override;
    virtual bool lt(const mstate& rhs) const override;
    virtual ~mstate_ncsb() override;

    friend class complement_ncsb;
  }; // mstate_ncsb }}}

public: // METHODS

  /// constructor
  complement_ncsb(const cmpl_info& info, unsigned scc_index);

  virtual mstate_set get_init() const override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual mstate_col_set get_succ_track_to_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual ~complement_ncsb() override;
}; // complement_ncsb }}}
} // namespace kofola }}}

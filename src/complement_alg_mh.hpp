// implementation of the Miyano & Hayashi complementation algorithm for
// inherently weak SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of the Miyano & Hayashi complementation algorithm for
/// inherently weak SCCs
class complement_mh : public abstract_complement_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given component
  class mstate_mh : public abstract_complement_alg::mstate
  { // {{{
  private: // DATA MEMBERS

    std::set<unsigned> states_;
    std::set<unsigned> breakpoint_;

  public: // METHODS

    /// constructor
    mstate_mh(const std::set<unsigned>& states, const std::set<unsigned>& breakpoint);

    virtual std::string to_string() const override;
    virtual ~mstate_mh() override;
  }; // mstate_mh }}}

public: // METHODS

  /// constructor
  complement_mh(const cmpl_info& info, unsigned scc_index);

  virtual mstate_set get_init() const override;
  virtual mstate_set get_succ_track(const mstate& src, const bdd& symbol) const override;
  virtual mstate_set get_succ_track_to_active(const mstate& src, const bdd& symbol) const override;
  virtual mstate_set get_succ_active(const mstate& src, const bdd& symbol) const override;

  virtual ~complement_mh() override;
}; // complement_mh }}}
} // namespace kofola }}}


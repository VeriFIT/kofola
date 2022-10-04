// implementation of determinization-based complementation algorithm for
// nondeterministic accepting SCCs based on Safra's trees and producing parity
// condition in the similar way as in
// Yong Li, Andrea Turrini, Weizhi Feng,
// Moshe Y. Vardi, Lijun Zhang: Divide-and-Conquer Determinization of Büchi
// Automata Based on SCC Decomposition. CAV (2) 2022: 152-173

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of determinization-based complementation algorithm for
/// nondeterministic accepting SCCs based on Safra's trees and producing parity
/// condition
class complement_safra : public abstract_complement_alg
{ // {{{
public: // TYPES

  /// partial macrostate for the given component
  class mstate_safra : public abstract_complement_alg::mstate
  { // {{{
  private: // DATA MEMBERS


  public: // METHODS

    /// constructor
    mstate_safra();

    virtual std::string to_string() const override;
    virtual bool is_active() const override { return true; }
    virtual bool eq(const mstate& rhs) const override;
    virtual bool lt(const mstate& rhs) const override;
    virtual ~mstate_safra() override;

    friend class complement_safra;
  }; // mstate_safra }}}

public: // METHODS

  /// constructor
  complement_safra(const cmpl_info& info, unsigned part_index);

  virtual mstate_set get_init() const override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual mstate_set lift_track_to_active(const mstate* src) const override;

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) const override;

  virtual bool use_round_robin() const override { return false; }

  virtual spot::acc_cond get_acc_cond() const override;

  virtual ~complement_safra() override;

}; // complement_safra }}}
} // namespace kofola }}}

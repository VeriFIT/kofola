// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of NCSB-based complementation algorithm for deterministic SCCs
class complement_sd_tela : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_sd_tela(const cmpl_info& info, unsigned part_index);

  virtual mstate_set get_init() override;

  virtual mstate_col_set get_succ_track(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol) override;

  virtual mstate_set lift_track_to_active(const mstate* src) override;

  virtual mstate_col_set get_succ_active(
    const std::set<unsigned>&  glob_reached,
    const mstate*              src,
    const bdd&                 symbol,
    bool resample = true) override;

  virtual bool use_round_robin() const override { return false; }

  virtual bool use_shared_breakpoint() const override { return this->info_.shared_breakpoint_; }

  virtual spot::acc_cond get_acc_cond() override
  { return spot::acc_cond(1, spot::acc_cond::inf({0})); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~complement_sd_tela() override { };

private:
  CondDNF acc_cond_ {};
}; // complement_sd_tela }}}

namespace sd_tela {

/// partial macrostate for the given component
class mstate_sd_tela : public abstract_complement_alg::mstate
{ // {{{
public: // DATA MEMBERS

  std::set<unsigned> check_ {};       // states for runs that need to be checked
  std::vector<std::set<unsigned>> safe_models_ {};        // safe states for models (cannot accept Fin colors)
  std::set<unsigned> breakpoint_;
  unsigned model_index_ {0};           // index of the model 0 --> universal quantification over runs; 1 --> first model, ....
  unsigned inf_index_ {0};             // index of the INF condition in model_index_-th model
  bool active_ {false};           // true = active ; false = track

public: // METHODS

  /// constructor
  mstate_sd_tela(
    const std::set<unsigned>&  check,
    const std::vector<std::set<unsigned>>&  safe_models,
    const std::set<unsigned>&  breakpoint,
    unsigned                   model_index,
    unsigned                   inf_index,
    bool                       active
  ) : check_(check),
    safe_models_(safe_models),
    breakpoint_(breakpoint),
    model_index_(model_index),
    inf_index_(inf_index),
    active_(active)
  { }

  mstate_sd_tela(
    const std::set<unsigned>&  check,
    size_t  num_safe_models,
    const std::set<unsigned>&  breakpoint,
    unsigned                   model_index,
    unsigned                   inf_index,
    bool                       active
  ) : check_(check),
    safe_models_(std::vector<std::set<unsigned>>(num_safe_models)),
    breakpoint_(breakpoint),
    model_index_(model_index),
    inf_index_(inf_index),
    active_(active)
  { } // mstate_sd_tela() }}}

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_sd_tela() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { this->breakpoint_ = get_set_intersection(breakpoint, this->check_); }

  virtual bool subsum_less_early(const mstate& rhs) override {
    // TODO: implement subsumption for SD-TELA
    return false;
  };

  friend class kofola::complement_sd_tela;
}; // mstate_sd_tela }}}

std::vector<mstate_sd_tela> guess_safe_models(const mstate_sd_tela& init, const std::set<unsigned>& states, unsigned num_models);

} // namespace sd_tela

} // namespace kofola }}}

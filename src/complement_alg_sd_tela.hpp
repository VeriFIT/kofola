// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

using SafeModels = std::vector<std::set<unsigned>>;

namespace sd_tela {

/// partial macrostate for the given component
class mstate_sd_tela : public abstract_complement_alg::mstate
{ // {{{
public: // DATA MEMBERS

  std::set<unsigned> check_ {};       // states for runs that need to be checked
  std::vector<std::set<unsigned>> safe_models_ {};        // safe states for models (cannot accept Fin colors)
  std::set<unsigned> breakpoint_ {};
  unsigned model_index_ {0};           // index of the model 0 --> universal quantification over runs; 1 --> first model, ....
  unsigned inf_index_ {0};             // index of the INF condition in model_index_-th model
  bool active_ {false};           // true = active ; false = track

public: // METHODS

  mstate_sd_tela() = default;

  /// constructor
  mstate_sd_tela(
    const std::set<unsigned>&  check,
    const SafeModels&  safe_models,
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
  
  // intersection with the breakpoint counterpart
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { 
    this->breakpoint_ = get_set_intersection(breakpoint, this->get_lift_breakpoint()); 
  }

  virtual bool subsum_less_early(const mstate& rhs) override {
    (void)rhs; // suppress unused parameter warning
    // TODO: implement subsumption for SD-TELA
    return false;
  };

public:

  /**
   * @brief Returns the lift breakpoint set for the current macrostate.
   *
   * If the macrostate represents universal quantification over runs (model_index_ == 0),
   * returns the set of states to be checked. Otherwise, returns the breakpoint set
   * for the specific safe model indexed by model_index_ - 1.
   *
   * @return The set of states representing the lift breakpoint.
   */
  std::set<unsigned> get_lift_breakpoint() const {
    // universal quantification over runs
    if(this->model_index_ == 0) {
      return this->check_;
    } else {
      // return the breakpoint of the model
      return this->safe_models_[this->model_index_ - 1];
    }
  }

  /**
   * @brief Increments the model and inf indices for the macrostate according to the acceptance condition.
   *
   * If currently at universal quantification (model_index_ == 0), increments model_index_ and resets inf_index_ to 0.
   * Otherwise, if inf_index_ is at the last inf condition for the current model, increments model_index_ and resets inf_index_ to 0.
   * If not at the last inf condition, simply increments inf_index_.
   *
   * @param acc_cond Acceptance condition (CondDNF) used to determine inf condition bounds.
   */
  void increment_indices(const CondDNF& acc_cond) {
    if(this->model_index_ == 0) {
      // increment model index
      this->model_index_++;
      this->inf_index_ = 0; // reset inf index
    } else {
      if(this->inf_index_ >= acc_cond[this->model_index_ - 1].infs.size() - 1) {
        // increment model index
        this->model_index_++;
        this->inf_index_ = 0; // reset inf index
      } else {
        // increment inf index
        this->inf_index_++;
      }
    }
  }
}; // mstate_sd_tela }}}


// FUNCTIONS

std::vector<mstate_sd_tela> guess_safe_models(const mstate_sd_tela& init, const std::set<unsigned>& states, unsigned num_models);

} // namespace sd_tela


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

// AUXILIARY METHODS
public:

  std::vector<std::pair<sd_tela::mstate_sd_tela, std::set<unsigned>>> get_succ_mst_track(
    const std::set<unsigned>& glob_reached,
    const sd_tela::mstate_sd_tela* src_mst,
    const bdd& symbol) const;

  bool contains_transition_color(
    const std::set<unsigned>& states, 
    const bdd& bdd, 
    const spot::acc_cond::mark_t& col) const;

  std::set<unsigned> get_succ_excluding_colors(
    const std::set<unsigned>& states,
    const bdd& bdd,
    const spot::acc_cond::mark_t& col) const;

  std::pair<SafeModels, std::set<unsigned>> get_safe_succ_reach(
    const SafeModels& safe_models, 
    const bdd& symbol) const;


  std::set<unsigned> get_succ_breakpoint_tmp(
    const sd_tela::mstate_sd_tela* src_mst,
    const std::set<unsigned>& succ_safe_reach,
    const bdd& symbol) const;

// DATA MEMBERS
private:
  CondDNF acc_cond_ {};
}; // complement_sd_tela }}}

} // namespace kofola }}}


#include "complement_alg_sd_tela.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // anonymous namespace {{{

/// partial macrostate for the given component
class mstate_sd_tela : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

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
} // anonymous namespace }}}

/**
 * @brief Returns a string representation of the macrostate for debugging and logging purposes.
 * @return String describing the macrostate.
 */
std::string mstate_sd_tela::to_string() const
{
  std::string res = std::string("[SD-TELA(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->check_);
  for (size_t i = 0; i < this->safe_models_.size(); ++i) {
    res += ", S" + std::to_string(i) + "=" + std::to_string(this->safe_models_[i]);
    if(i < this->safe_models_.size() - 1) {
      res += ", ";
    }
  }
  res += ", ModelI=" + std::to_string(this->model_index_);
  res += ", InfI=" + std::to_string(this->inf_index_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

/**
 * @brief Checks equality of this macrostate with another.
 * @param rhs The macrostate to compare with.
 * @return True if all relevant data members are equal, false otherwise.
 */
bool mstate_sd_tela::eq(const mstate& rhs) const
{
  const mstate_sd_tela* rhs_sd_tela = dynamic_cast<const mstate_sd_tela*>(&rhs);
  assert(rhs_sd_tela);
  return (this->active_ == rhs_sd_tela->active_) &&
    (this->model_index_ == rhs_sd_tela->model_index_) &&
    (this->inf_index_ == rhs_sd_tela->inf_index_) &&
    (this->check_ == rhs_sd_tela->check_) &&
    (this->breakpoint_ == rhs_sd_tela->breakpoint_) &&
    (this->safe_models_ == rhs_sd_tela->safe_models_);
}


/**
 * @brief Defines a strict weak ordering for macrostates.
 * @param rhs The macrostate to compare with.
 * @return True if this macrostate is less than the other, used for set/map ordering.
 */
bool mstate_sd_tela::lt(const mstate& rhs) const
{ // {{{
  const mstate_sd_tela* rhs_sd_tela = dynamic_cast<const mstate_sd_tela*>(&rhs);
  assert(rhs_sd_tela);

  if (this->active_ != rhs_sd_tela->active_) { return this->active_ < rhs_sd_tela->active_; }
  if (this->model_index_ != rhs_sd_tela->model_index_) { return this->model_index_ < rhs_sd_tela->model_index_; }
  if (this->inf_index_ != rhs_sd_tela->inf_index_) { return this->inf_index_ < rhs_sd_tela->inf_index_; }
  if (this->check_ != rhs_sd_tela->check_) { return this->check_ < rhs_sd_tela->check_; }
  if (this->breakpoint_ != rhs_sd_tela->breakpoint_) { return this->breakpoint_ < rhs_sd_tela->breakpoint_; }
  if (this->safe_models_ != rhs_sd_tela->safe_models_) { return this->safe_models_ < rhs_sd_tela->safe_models_; }

  return false;   // if all are equal
} // lt() }}}


/**
 * @brief Constructor for the SD-TELA complementation algorithm.
 * @param info Automaton and complementation info.
 * @param part_index Index of the SCC partition for this instance.
 *
 * The partition index identifies the strongly connected component (SCC) of the automaton
 * that this complementation instance operates on.
 */
complement_sd_tela::complement_sd_tela(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
{ }

/**
 * @brief Returns the initial set of macrostates for the SD-TELA complementation algorithm.
 *
 * Only includes the initial state if it belongs to the current SCC partition (identified by part_index).
 * @return Set of initial macrostates.
 */
mstate_set complement_sd_tela::get_init()
{ // {{
  DEBUG_PRINT_LN("init SD-TELA for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state.insert(orig_init);
  }

  std::shared_ptr<mstate> ms(new mstate_sd_tela(init_state, {}, {}, 0, 0, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

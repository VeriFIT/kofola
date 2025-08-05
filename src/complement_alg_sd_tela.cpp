#include "complement_alg_sd_tela.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace kofola {

namespace sd_tela {

/**
 * @brief Generates all possible assignments of the given states to the safe models.
 *
 * For each state in the input set, the function creates new macrostates by assigning the state
 * to each of the available models (from 0 to num_models-1). The result is a vector containing
 * all combinations where each state is assigned to one model, and all states are distributed
 * across the models. Used for exploring all possible safe model configurations.
 *
 * @param init The initial macrostate to start from.
 * @param states The set of states to assign to models.
 * @param num_models The number of models to distribute states into.
 * @return Vector of macrostates with all possible safe model assignments.
 */
std::vector<mstate_sd_tela> guess_safe_models(const mstate_sd_tela& init, const std::set<unsigned>& states, unsigned num_models) {
  std::vector<mstate_sd_tela> safe_models;
  std::deque<std::pair<mstate_sd_tela, unsigned>> queue;
  queue.push_back({init, 0});

  std::vector<unsigned> state_vec(states.begin(), states.end());

  while(!queue.empty()) {
    auto [macrostate, state_idx] = queue.front();
    queue.pop_front();
    if(state_idx >= state_vec.size()) {
      safe_models.push_back(macrostate);
      continue;
    }

    for(unsigned i = 0; i < num_models; i++) {
      mstate_sd_tela new_macrostate {macrostate};
      new_macrostate.safe_models_[i].insert(state_vec[state_idx]);
      queue.push_back({new_macrostate, state_idx + 1});
    }
  }
  return safe_models;
}

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

} // namespace sd_tela

/**
 * @brief Constructor for the SD-TELA complementation algorithm.
 * @param info Automaton and complementation info.
 * @param part_index Index of the SCC partition for this instance.
 *
 * The partition index identifies the strongly connected component (SCC) of the automaton
 * that this complementation instance operates on.
 */
complement_sd_tela::complement_sd_tela(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index) { 
  
  this->acc_cond_ = info.compute_cond_to_verify();
  for(size_t i = 0; i < this->acc_cond_.size(); ++i) {
    this->acc_cond_[i].simplify();
  }
}

/**
 * @brief Computes, for each safe model, the set of successor states in the same SCC reachable via transitions whose condition is implied by the given BDD symbol.
 * Ensures that states already reached by previous models are excluded from subsequent models' successors.
 * Returns a pair:
 *   - Vector of sets: each set contains the safe successors for the corresponding model (excluding previously reached states).
 *   - Set of all states reached by any model (union of all safe successors).
 *
 * @param safe_models Vector of sets, each representing a safe model's states.
 * @param symbol BDD condition representing the transition label.
 * @return Pair of (vector of safe successor sets per model, set of all reached states).
 */
std::pair<SafeModels, std::set<unsigned>> complement_sd_tela::get_safe_succ_reach(const SafeModels& safe_models, const bdd& symbol) const {
  SafeModels safe_models_ret;
  std::set<unsigned> reached_states {};

  // Get all safe successors for each safe model
  for (const auto& safe_model : safe_models) {
    std::set<unsigned> succ = kofola::get_all_successors_in_scc(
      this->info_.aut_, this->info_.scc_info_, safe_model, symbol);
    
    // remove states that are already in previous check set
    // we have made a wrong guess --> we could move those states to different safe set or waint in breakpoint
    succ = kofola::get_set_difference(succ, reached_states);
    safe_models_ret.push_back(succ);
    reached_states.insert(succ.begin(), succ.end());
  }

  return {safe_models_ret, reached_states};
}

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

  std::shared_ptr<mstate> ms(new sd_tela::mstate_sd_tela(init_state, {}, {}, 0, 0, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

/**
  * @brief Computes the tracking successors for a given macrostate and input symbol in SD-TELA complementation.
  *
  * For each safe model, checks if there are transitions labeled by a Fin condition; if any such transition exists,
  * returns an empty result (no valid successors). Otherwise, computes the successors of safe models and the check set
  * over the given symbol, ensuring that states already in safe models are excluded from the check set. Constructs a new
  * macrostate with updated check and safe sets, and returns it as the only element in the result.
  *
  * @param glob_reached Set of all states reached over the symbol.
  * @param src Source macrostate (should be a tracking state).
  * @param symbol BDD condition representing the transition label.
  * @return Set of macrostate-color pairs representing the possible successors.
  */
mstate_col_set complement_sd_tela::get_succ_track(
    const std::set<unsigned>& glob_reached,
    const mstate* src,
    const bdd& symbol)
{
  const sd_tela::mstate_sd_tela* src_mst = dynamic_cast<const sd_tela::mstate_sd_tela*>(src);
  assert(src_mst);
  assert(!src_mst->active_);

  // check if there are NO transitions labeled by Fin condition in 
  // each safe model
  for(size_t i = 0; i < src_mst->safe_models_.size(); ++i) {
    assert(this->acc_cond_[i].fins.size() <= 1);
    for(size_t fin_cond = 0; fin_cond < this->acc_cond_[i].fins.size(); fin_cond++) {
      if (contains_transition_color(src_mst->safe_models_[i], symbol, this->acc_cond_[i].fins[fin_cond])) {
        return {};
      }
    }
  }

  // successors of safe models over a symbol
  auto [succ_safe, safe_reach] = this->get_safe_succ_reach(src_mst->safe_models_, symbol);

  // successors of check set (do not include states that are already in safe models)
  std::set<unsigned> succ_check;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == static_cast<int>(this->part_index_)) {
      if (safe_reach.find(st) == safe_reach.end()) { // if not in safe
        succ_check.insert(st);
      }
    }
  }

  std::shared_ptr<mstate> ms(new sd_tela::mstate_sd_tela(succ_check, succ_safe, {}, src_mst->model_index_, src_mst->inf_index_, false));
  // breakpoint is empty
  mstate_col_set result = {{ms, {}}}; 
  
  return result;
}

/**
 * @brief Lifts a tracking macrostate to an active macrostate in SD-TELA complementation.
 *
 * Converts the given tracking macrostate (tracking) into an active macrostate, preserving the check set,
 * safe models, model and inf indices, and setting the breakpoint to the lift breakpoint as determined by the macrostate.
 * This operation does not increment model or inf counters, assuming they are already set to the currently tracked entities.
 *
 * @param src Source macrostate (must be a tracking state).
 * @return Set containing the new active macrostate.
 */
mstate_set complement_sd_tela::lift_track_to_active(const mstate* src) {
  const sd_tela::mstate_sd_tela* src_mst = dynamic_cast<const sd_tela::mstate_sd_tela*>(src);
  assert(src_mst);
  assert(!src_mst->active_);

  // do not increment model or inf counters 
  // we assume that the counters are set to currently tracked entities (safe models or runs)
  std::shared_ptr<mstate> ms(new sd_tela::mstate_sd_tela(src_mst->check_, src_mst->safe_models_, src_mst->get_lift_breakpoint(), src_mst->model_index_, src_mst->inf_index_, true));
  return {ms};
}

mstate_col_set complement_sd_tela::get_succ_active(
    const std::set<unsigned>& /*glob_reached*/,
    const mstate* /*src*/,
    const bdd& /*symbol*/,
    bool /*resample*/)
{
    return mstate_col_set{};
}

/**
  * @brief Checks if any transition from the given states under the specified BDD condition
  *        contains the given acceptance color.
  *
  * Iterates over all outgoing transitions from the provided states. For transitions that remain
  * within the same SCC and whose condition is implied by the given BDD, checks if the transition's
  * acceptance set contains the specified color. Returns true if at least one such transition exists.
  *
  * @param states Set of source states to check transitions from.
  * @param bdd BDD condition that transitions must satisfy.
  * @param col Acceptance color to look for in transitions.
  * @return True if any transition matches the criteria, false otherwise.
 */
bool complement_sd_tela::contains_transition_color(const std::set<unsigned>& states, const bdd& bdd, const spot::acc_cond::mark_t& col) const {
  for (unsigned s : states) {
    for (const auto &t : this->info_.aut_->out(s)) {
      if (this->info_.scc_info_.scc_of(s) == this->info_.scc_info_.scc_of(t.dst) && bdd_implies(bdd, t.cond)) {
        if (t.acc & col) { return true; }
      }
    }
  }

  return false;
}

} // namespace kofola
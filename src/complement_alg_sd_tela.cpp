#include "complement_alg_sd_tela.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace kofola {

namespace sd_tela {

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

  std::shared_ptr<mstate> ms(new sd_tela::mstate_sd_tela(init_state, {}, {}, 0, 0, false));
  mstate_set result = {ms};
  return result;
} // get_init() }}}

// Trivial implementations for complement_sd_tela methods
mstate_col_set complement_sd_tela::get_succ_track(
    const std::set<unsigned>& /*glob_reached*/,
    const mstate* /*src*/,
    const bdd& /*symbol*/)
{
    return mstate_col_set{};
}

mstate_set complement_sd_tela::lift_track_to_active(const mstate* /*src*/)
{
    return mstate_set{};
}

mstate_col_set complement_sd_tela::get_succ_active(
    const std::set<unsigned>& /*glob_reached*/,
    const mstate* /*src*/,
    const bdd& /*symbol*/,
    bool /*resample*/)
{
    return mstate_col_set{};
}

} // namespace kofola
// implementation of the initial-almost-deterministic complement algorithm

#include "complement_alg_iadacs.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // {{{

/// partial macrostate for the given component
class mstate_init_almost_det : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::set<unsigned> states_;
  bool active_;
  std::set<unsigned> empty_set_;

public: // METHODS

  /// constructor
  mstate_init_almost_det(
    const std::set<unsigned>&  states,
    // const std::set<unsigned>&  breakpoint,
    bool                       active
  ) : states_(states),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_init_almost_det() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override { return this->empty_set_; }
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { (void)breakpoint; }

  friend class kofola::complement_init_almost_det;
}; // mstate_init_almost_det }}}


std::string mstate_init_almost_det::to_string() const
{
  std::string res = std::string("[INIT_ALMOST_DET(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->states_);
  res += "]";
  return res;
}

bool mstate_init_almost_det::eq(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);
  return (this->states_ == rhs_mh->states_);
}

bool mstate_init_almost_det::lt(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);

  if (this->states_ != rhs_mh->states_) { return this->states_ < rhs_mh->states_; }

  return false;   // if all are equal
}

} // anonymous namespace }}}

complement_init_almost_det::complement_init_almost_det(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
{ 
    dualized_acc_cond_ = info.aut_->get_acceptance().complement();
    min_colour_ = info.aut_->acc().all_sets().min_set() - 1;
    DEBUG_PRINT_LN("Created complement_init_almost_det for part " + std::to_string(part_index) + " with acc cond: " + std::to_string(dualized_acc_cond_) + " and min color: " + std::to_string(min_colour_));
}

mstate_set complement_init_almost_det::get_init()
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state.insert(orig_init);
  }

  mstate_set result;
  std::shared_ptr<mstate> ms(new mstate_init_almost_det(init_state, false));
  result.push_back(ms);

  return result;
} // get_init() }}}

mstate_col_set complement_init_almost_det::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{ // {{{
  DEBUG_PRINT_LN("Init det successor");
  DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("src = " + std::to_string(*src));
  DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

  const mstate_init_almost_det* src_iad = dynamic_cast<const mstate_init_almost_det*>(src);
  assert(src_iad);
  assert(!src_iad->active_);
  (void)src_iad; // make release happy

  std::set<unsigned> states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == static_cast<int>(this->part_index_)) {
      states.insert(st);
    }
  }

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(states, false));
  return {{ms, {}}};
} // get_succ_track() }}}

mstate_set complement_init_almost_det::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_init_almost_det* src_iad = dynamic_cast<const mstate_init_almost_det*>(src);
  assert(src_iad);
  assert(!src_iad->active_);

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(src_iad->states_, true));
  return {ms};
} // lift_track_to_active() }}}

mstate_col_set complement_init_almost_det::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol,
  bool                       resample)
{
  (void)resample; // make release happy
  const mstate_init_almost_det* src_iad = dynamic_cast<const mstate_init_almost_det*>(src);
  assert(src_iad);
  assert(src_iad->active_);

  DEBUG_PRINT_LN("tracking successor of: " + std::to_string(*src_iad));
  mstate_init_almost_det tmp(src_iad->states_, false);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_init_almost_det* track_ms = dynamic_cast<const mstate_init_almost_det*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  bool generate_condition = false;
  spot::acc_cond::mark_t acc = {};
  std::set<unsigned> seen_cols;
  for (unsigned s : src_iad->states_) {
    for (const auto &t : this->info_.aut_->out(s)) {
      if (bdd_implies(symbol, t.cond)) { 
        if (this->info_.scc_info_.scc_of(t.dst) == this->info_.scc_info_.scc_of(s)) {
          acc |= t.acc;
          for(auto col: t.acc.sets()) {
            seen_cols.insert(col);
          }
        } 
      }
    }
  }
  
  generate_condition = (acc != spot::acc_cond::mark_t({}));
  mstate_col_set result;

  if (generate_condition) {
    // round-robin not used
    std::shared_ptr<mstate> ms(new mstate_init_almost_det(track_ms->states_, true));
      result.push_back({ms, seen_cols});
  }
  else {
    std::shared_ptr<mstate> ms(new mstate_init_almost_det(track_ms->states_, true));
    result.push_back({ms, {}});
  }

  return result;
}

complement_init_almost_det::~complement_init_almost_det()
{ }

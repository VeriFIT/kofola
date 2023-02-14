// implementation of the Miyano & Hayashi complementation algorithm for
// inherently weak SCCs

#include "complement_alg_init_det.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // {{{

/// partial macrostate for the given component
class mstate_init_det : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::set<unsigned> states_;
  std::set<unsigned> breakpoint_;
  bool active_;

public: // METHODS

  /// constructor
  mstate_init_det(
    const std::set<unsigned>&  states,
    const std::set<unsigned>&  breakpoint,
    bool                       active
  ) : states_(states),
    breakpoint_(breakpoint),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_init_det() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { this->breakpoint_ = get_set_intersection(breakpoint, this->states_); }

  friend class kofola::complement_init_det;
}; // mstate_init_det }}}


std::string mstate_init_det::to_string() const
{
  std::string res = std::string("[INIT_DET(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->states_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

bool mstate_init_det::eq(const mstate& rhs) const
{
  const mstate_init_det* rhs_mh = dynamic_cast<const mstate_init_det*>(&rhs);
  assert(rhs_mh);
  return (this->states_ == rhs_mh->states_) &&
    (this->breakpoint_ == rhs_mh->breakpoint_);
}

bool mstate_init_det::lt(const mstate& rhs) const
{
  const mstate_init_det* rhs_mh = dynamic_cast<const mstate_init_det*>(&rhs);
  assert(rhs_mh);

  if (this->states_ != rhs_mh->states_) { return this->states_ < rhs_mh->states_; }
  if (this->breakpoint_ != rhs_mh->breakpoint_) { return this->breakpoint_ < rhs_mh->breakpoint_; }

  return false;   // if all are equal
}

} // anonymous namespace }}}

complement_init_det::complement_init_det(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
{ }

mstate_set complement_init_det::get_init()
{ // {{{
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_state.insert(orig_init);
  }

  mstate_set result;
  std::shared_ptr<mstate> ms(new mstate_init_det(init_state, {}, false));
  result.push_back(ms);

  return result;
} // get_init() }}}

mstate_col_set complement_init_det::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol)
{ // {{{
  DEBUG_PRINT_LN("Init det successor");
  DEBUG_PRINT_LN("glob_reached = " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("src = " + std::to_string(*src));
  DEBUG_PRINT_LN("symbol = " + std::to_string(symbol));

  const mstate_init_det* src_mh = dynamic_cast<const mstate_init_det*>(src);
  assert(src_mh);
  assert(!src_mh->active_);

  std::set<unsigned> states;
  for (unsigned st : glob_reached) {
    if (this->info_.st_to_part_map_.at(st) == this->part_index_) {
      states.insert(st);
    }
  }

  std::shared_ptr<mstate> ms(new mstate_init_det(states, {}, false));
  return {{ms, {}}};
} // get_succ_track() }}}

mstate_set complement_init_det::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_init_det* src_mh = dynamic_cast<const mstate_init_det*>(src);
  assert(src_mh);
  assert(!src_mh->active_);

  std::shared_ptr<mstate> ms(new mstate_init_det(src_mh->states_, {}, true));
  return {ms};
} // lift_track_to_active() }}}

mstate_col_set complement_init_det::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol,
  bool                       resample)
{
  const mstate_init_det* src_mh = dynamic_cast<const mstate_init_det*>(src);
  assert(src_mh);
  assert(src_mh->active_);

  DEBUG_PRINT_LN("tracking successor of: " + std::to_string(*src_mh));
  mstate_init_det tmp(src_mh->states_, {}, false);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_init_det* track_ms = dynamic_cast<const mstate_init_det*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  std::set<unsigned> succ_break = kofola::get_all_successors_in_scc(
    this->info_.aut_, this->info_.scc_info_, src_mh->breakpoint_, symbol);

  // intersect with what is really reachable (for simulation pruning)
  succ_break = kofola::get_set_intersection(succ_break, glob_reached);

  bool generate_condition = false;
  spot::acc_cond::mark_t acc = {0};
  for (unsigned s : src_mh->states_) {
    for (const auto &t : this->info_.aut_->out(s)) {
      if (bdd_implies(symbol, t.cond)) { 
        if (this->info_.scc_info_.scc_of(t.dst) == this->info_.scc_info_.scc_of(s)) {
          if (t.acc == acc) {
            generate_condition = true;
            break;
          }
        } 
      }
    }
    if (generate_condition) break;
  }

  mstate_col_set result;
  if (generate_condition) {
    if (this->use_round_robin()) {
      std::shared_ptr<mstate> ms(new mstate_init_det(track_ms->states_, {}, false));
      result.push_back({ms, {1}});
    } else { // no round robin
      std::shared_ptr<mstate> ms(new mstate_init_det(track_ms->states_, {}, true));
      result.push_back({ms, {1}});
    }
  }
  else {
    std::shared_ptr<mstate> ms(new mstate_init_det(track_ms->states_, {}, true));
    result.push_back({ms, {}});
  }

  return result;
}

complement_init_det::~complement_init_det()
{ }

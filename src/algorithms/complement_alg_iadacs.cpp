// implementation of the initial-almost-deterministic complementation algorithm

#include "complement_alg_iadacs.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;


namespace { // {{{

/// partial macrostate for the given component
class mstate_init_almost_det : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::vector<unsigned> compressed_mapping_;
  std::vector<unsigned> uncompressed_mapping_;
  bool active_;
  std::set<unsigned> empty_set_;
  unsigned bound_;

public: // METHODS

  /// constructor
  mstate_init_almost_det(
    // const std::set<unsigned>&  breakpoint,
    const std::vector<unsigned>& compressed_mapping,
    const std::vector<unsigned>& uncompressed_mapping,
    bool                       active,
    unsigned                   bound = 0
  ) : 
    compressed_mapping_(compressed_mapping), uncompressed_mapping_(uncompressed_mapping), active_(active), bound_(bound)
  {
  }

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
  res += "MAPPING=" + std::to_string(this->compressed_mapping_) + ", ";
  res += "]";
  return res;
}

bool mstate_init_almost_det::eq(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);
  return (this->compressed_mapping_ == rhs_mh->compressed_mapping_);
}

bool mstate_init_almost_det::lt(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);

  if (this->compressed_mapping_ != rhs_mh->compressed_mapping_) { return this->compressed_mapping_ < rhs_mh->compressed_mapping_; }

  return false;   // if all are equal
}

} // anonymous namespace }}}

determinisation_acc_cond::determinisation_acc_cond(const spot::acc_cond::acc_code& acc_cond, unsigned disjuncts) {
  auto template_code = acc_cond;
  auto max_col = template_code.used_sets().max_set();

  template_code &= spot::acc_cond::acc_code::fin({max_col}); // add Fin() clause for discontinuation event

  acc_code_ = spot::acc_cond::acc_code::f(); // set to neutral element for disjunction
  disj_size_ = template_code.used_sets().count();
  for(unsigned i = 0; i < disjuncts; ++i) {
    acc_code_ |= (template_code << (disj_size_ * i));
    additional_fins_.push_back(max_col + disj_size_ * i);
  }
  DEBUG_PRINT_LN("Created determinisation_acc_cond with acc cond: " + std::to_string(acc_code_));
}

spot::acc_cond determinisation_acc_cond::get_acc_cond() const {
    return spot::acc_cond(acc_code_);
}

unsigned determinisation_acc_cond::get_min_colour() const {
    return acc_code_.used_sets().min_set() - 1;
}

spot::acc_cond::mark_t determinisation_acc_cond::get_all_discontinuation_colours() const {
    spot::acc_cond::mark_t additional_fins = {};
    for(auto fin: additional_fins_) {
        additional_fins.set(fin);
    }

    return additional_fins;
}

complement_init_almost_det::complement_init_almost_det(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
  , runs_bound_(count_part_states(info, part_index))
  , acc_cond_(info.aut_->get_acceptance(), runs_bound_)
{ }

unsigned complement_init_almost_det::count_part_states(const cmpl_info& info, unsigned part_index) {
  unsigned count = 0;
  for (unsigned i = 0; i < info.aut_->num_states(); ++i) {
    if (info.st_to_part_map_.at(i) == static_cast<int>(part_index)) {
      ++count;
    }
  }
  return count;
}

mstate_set complement_init_almost_det::get_init()
{ // {{{
  unsigned bound = this->runs_bound_;
  std::vector<unsigned> init_state(bound, UINT_MAX);

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == static_cast<int>(this->part_index_)) {
    init_state[0] = orig_init;
  }

  mstate_set result;
  std::shared_ptr<mstate> ms(new mstate_init_almost_det(init_state, init_state, false, bound));
  result.push_back(ms);

  return result;
} // get_init() }}}

std::vector<unsigned> compress(const std::vector<unsigned>& mapping) {
  std::vector<unsigned> compressed(mapping.size(), UINT_MAX); // initialize with "infinite"

  for (unsigned i = 0; i < mapping.size(); ++i) {
    if (mapping[i] == UINT_MAX) { 
      continue; // compress
    }
    compressed[i] = mapping[i];
  }

  return compressed;
}

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

  const unsigned bound = src_iad->bound_;
  std::vector<unsigned> g(bound, UINT_MAX);
  std::unordered_set<unsigned> defined_runs;
  unsigned i = 0;

  for (unsigned run : src_iad->compressed_mapping_) {
    if (run == UINT_MAX) { // end of defined runs in the mapping
      break;
    }

    const auto run_scc = this->info_.scc_info_.scc_of(run);
    for (const auto& t : this->info_.aut_->out(run)) {
      if (!bdd_implies(symbol, t.cond)) { // not successor over the symbol
        continue;
      }

      if (this->info_.scc_info_.scc_of(t.dst) != run_scc) { // not successor in the same SCC
        continue;
      }

      if (defined_runs.find(t.dst) != defined_runs.end()) { // already defined in the mapping
        continue;
      }

      g[i] = t.dst;
      defined_runs.insert(t.dst);
      break; // only one successor in the same deterministic SCC
    }
    i++; // always increase i to keep the order of runs in the mapping, even if some runs are disconnected
  }

  const auto& sccs_in_part = this->info_.part_to_scc_map_.at(this->part_index_);
  for (unsigned run : glob_reached) {
    const auto run_scc = this->info_.scc_info_.scc_of(run);

    if (sccs_in_part.find(run_scc) == sccs_in_part.end()) {
      continue;
    }

    if (defined_runs.find(run) != defined_runs.end()) {
      continue;
    }

    g[i++] = run;
    defined_runs.insert(run);
  }

  auto f = compress(g);

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(f, g, false, bound));

  auto disc_cols = acc_cond_.get_all_discontinuation_colours();
  auto colors_iterable = disc_cols.sets();
  std::set<unsigned> colors_set(colors_iterable.begin(), colors_iterable.end());

  return {{ms, colors_set}};
} // get_succ_track() }}}

mstate_set complement_init_almost_det::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_init_almost_det* src_iad = dynamic_cast<const mstate_init_almost_det*>(src);
  assert(src_iad);
  assert(!src_iad->active_);

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(src_iad->compressed_mapping_, src_iad->uncompressed_mapping_, true, src_iad->bound_));
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
  mstate_init_almost_det tmp(src_iad->compressed_mapping_, src_iad->uncompressed_mapping_, false, src_iad->bound_);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_init_almost_det* track_ms = dynamic_cast<const mstate_init_almost_det*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  std::set<unsigned> colors_set;
  unsigned matching_prefix_len;
  for(matching_prefix_len = 0; matching_prefix_len < this->runs_bound_; ++matching_prefix_len) {
    const unsigned tracked_run = track_ms->compressed_mapping_[matching_prefix_len];
    const unsigned tracked_uncompressed = track_ms->uncompressed_mapping_[matching_prefix_len];
    const unsigned src_run = src_iad->compressed_mapping_[matching_prefix_len];

    if (tracked_run != tracked_uncompressed || tracked_run == UINT_MAX) {
      break;
    }

    if (src_run == UINT_MAX) {
      break;
    }

    for (const auto& t : this->info_.aut_->out(src_run)) {
      DEBUG_PRINT_LN("b = " + std::to_string(matching_prefix_len));

      if (!bdd_implies(symbol, t.cond)) {
        continue;
      }

      if (t.dst != tracked_run) {
        continue;
      }

      for (auto col : t.acc.sets()) {
        colors_set.insert(acc_cond_.map_colour(col, matching_prefix_len));
      }
      break;
    }
  }


  for(unsigned i = matching_prefix_len; i < this->runs_bound_; ++i) {
    colors_set.insert(acc_cond_.get_fin_mark(i));
  }

  
  mstate_col_set result;
  std::shared_ptr<mstate> ms(new mstate_init_almost_det(track_ms->compressed_mapping_, track_ms->uncompressed_mapping_, true, track_ms->bound_));

  result.push_back({ms, colors_set});  

  return result;
}

complement_init_almost_det::~complement_init_almost_det()
{ }

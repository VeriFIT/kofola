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

  std::vector<unsigned> f_; // compressed g_
  std::vector<unsigned> g_;
  bool active_;
  std::set<unsigned> empty_set_;
  unsigned bound_;

public: // METHODS

  /// constructor
  mstate_init_almost_det(
    // const std::set<unsigned>&  breakpoint,
    const std::vector<unsigned>& f, // compressed g
    const std::vector<unsigned>& g,
    bool                       active,
    unsigned                   bound = 0
  ) : 
    f_(f), g_(g), active_(active), bound_(bound)
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
  res += "ORDER=" + std::to_string(this->f_) + ", ";
  res += "]";
  return res;
}

bool mstate_init_almost_det::eq(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);
  return (this->f_ == rhs_mh->f_);
}

bool mstate_init_almost_det::lt(const mstate& rhs) const
{
  const mstate_init_almost_det* rhs_mh = dynamic_cast<const mstate_init_almost_det*>(&rhs);
  assert(rhs_mh);

  if (this->f_ != rhs_mh->f_) { return this->f_ < rhs_mh->f_; }

  return false;   // if all are equal
}

} // anonymous namespace }}}

determinisation_acc_cond::determinisation_acc_cond(const spot::acc_cond::acc_code& acc_cond, unsigned disjuncts) {
  auto template_code = acc_cond;
  auto max_col = template_code.used_sets().max_set();
  template_code &= spot::acc_cond::acc_code::fin({max_col});

  acc_code_ &= spot::acc_cond::acc_code::f();

  auto cols_used = template_code.used_sets().count();
  disj_size_ = cols_used;
  for(unsigned i = 0; i < disjuncts; ++i) {
    acc_code_ |= (template_code << (cols_used * i));
    additional_fins_.push_back(max_col + cols_used * i);
  }
  DEBUG_PRINT_LN("Created determinisation_acc_cond with acc cond: " + std::to_string(acc_code_));

  dnf_ = cmpl_info::acc_code_dnf(acc_code_);
}

spot::acc_cond determinisation_acc_cond::get_acc_cond() const {
    return spot::acc_cond(acc_code_);
}

unsigned determinisation_acc_cond::get_min_colour() const {
    return acc_code_.used_sets().min_set() - 1;
}

spot::acc_cond::mark_t determinisation_acc_cond::get_additional_fins_mark() const {
    spot::acc_cond::mark_t additional_fins = {};
    for(const auto& disj : dnf_) {
        for(const auto& fin : disj.fins) {
            additional_fins |= fin;
        }
    }

    return additional_fins;
}

complement_init_almost_det::complement_init_almost_det(const cmpl_info& info, unsigned part_index)
  : abstract_complement_alg(info, part_index)
  , part_states_(count_part_states(info, part_index))
  , acc_cond_(info.aut_->get_acceptance(), part_states_)
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
  unsigned bound = this->part_states_;
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

std::vector<unsigned> compress(const std::vector<unsigned>& g) {
  std::vector<unsigned> compressed(g.size(), UINT_MAX); // initialize with invalid state number
  for (unsigned i = 0; i < g.size(); ++i) {
    if (g[i] == UINT_MAX) { continue; } // compress
    compressed[i] = g[i];
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
  (void)src_iad; // make release happy
  
  auto bound = src_iad->bound_;
  
  std::vector<unsigned> g(bound, UINT_MAX); // initialize with "inf"
  std::unordered_set<unsigned> defined_runs; // to avoid duplicates in g
  
  unsigned i = 0;
  for(auto run: src_iad->f_) {
    if(run == UINT_MAX) { break; } // reached the end of defined runs in f
    for (const auto &t : this->info_.aut_->out(run)) {
      if (bdd_implies(symbol, t.cond)) { 
        if (this->info_.scc_info_.scc_of(t.dst) == this->info_.scc_info_.scc_of(run)) {
          if (defined_runs.find(t.dst) == defined_runs.end()) {
            g[i++] = t.dst;
            defined_runs.insert(t.dst);
          }
        } 
      }
    }
  }


  for(auto run: glob_reached) {
    auto scc_of_run = this->info_.scc_info_.scc_of(run);
    auto sccs_in_part = this->info_.part_to_scc_map_.at(this->part_index_);
    if (sccs_in_part.find(scc_of_run) != sccs_in_part.end()) {
        if (defined_runs.find(run) == defined_runs.end()) {
          g[i++] = run;
          defined_runs.insert(run);
        }
    }
  }

  auto f = compress(g);

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(f, g, false, bound));
  auto fins_mark = acc_cond_.get_additional_fins_mark();
  auto colors_iterable = fins_mark.sets(); 
  std::set<unsigned> colors_set(colors_iterable.begin(), colors_iterable.end());

  return {{ms, colors_set}};
} // get_succ_track() }}}

mstate_set complement_init_almost_det::lift_track_to_active(const mstate* src)
{ // {{{
  const mstate_init_almost_det* src_iad = dynamic_cast<const mstate_init_almost_det*>(src);
  assert(src_iad);
  assert(!src_iad->active_);

  std::shared_ptr<mstate> ms(new mstate_init_almost_det(src_iad->f_, src_iad->g_, true, src_iad->bound_));
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
  mstate_init_almost_det tmp(src_iad->f_, src_iad->g_, false, src_iad->bound_);
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_init_almost_det* track_ms = dynamic_cast<const mstate_init_almost_det*>(track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  std::set<unsigned> colors_set;
  unsigned b;
  for(b = 0; b < this->part_states_; ++b) {
    if (track_ms->f_[b] != track_ms->g_[b] || track_ms->f_[b] == UINT_MAX) { break; }

    if(src_iad->f_[b] == UINT_MAX) { break; }

    for (const auto &t : this->info_.aut_->out(src_iad->f_[b])) {
      DEBUG_PRINT_LN("b = " + std::to_string(b));
      if (bdd_implies(symbol, t.cond)) { 
        if (t.dst == track_ms->f_[b]) {
          for(auto col: t.acc.sets()) {
            colors_set.insert(acc_cond_.map_colour(col, b));
          }
          break;
        } 
      }
    }
  }


  for(unsigned i = b; i < this->part_states_; ++i) {
    colors_set.insert(acc_cond_.get_fin_mark(i));
  }

  
  mstate_col_set result;
  std::shared_ptr<mstate> ms(new mstate_init_almost_det(track_ms->f_, track_ms->g_, true, track_ms->bound_));

  result.push_back({ms, colors_set});  

  return result;
}

complement_init_almost_det::~complement_init_almost_det()
{ }

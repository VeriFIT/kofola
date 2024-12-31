#pragma once

#include "abstract_complement_alg.hpp"
#include "complement_alg_ncsb.hpp"

namespace kofola { // {{{

class complement_ncsb_delay;

namespace { // anonymous namespace {{{

/// partial macrostate for the given component
class mstate_ncsb : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  std::set<unsigned> check_;       // states for runs that need to be checked
  std::set<unsigned> safe_;        // safe states (cannot see accepting transition)
  std::set<unsigned> breakpoint_;
  bool active_;                    // true = active ; false = track

public: // METHODS

  /// constructor
  mstate_ncsb(
    const std::set<unsigned>&  check,
    const std::set<unsigned>&  safe,
    const std::set<unsigned>&  breakpoint,
    bool                       active
  ) : check_(check),
    safe_(safe),
    breakpoint_(breakpoint),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_ncsb() override { }

  virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
  virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { this->breakpoint_ = get_set_intersection(breakpoint, this->check_); }

  friend class kofola::complement_ncsb_delay;
}; // mstate_ncsb }}}

} // anonymous namespace }}}

std::string mstate_ncsb::to_string() const
{
  std::string res = std::string("[NCSB(") + ((this->active_)? "A" : "T") + "): ";
  res += "C=" + std::to_string(this->check_);
  res += ", S=" + std::to_string(this->safe_);
  if (this->active_) {
    res += ", B=" + std::to_string(this->breakpoint_);
  }
  res += "]";
  return res;
}

bool mstate_ncsb::eq(const mstate& rhs) const
{
  const mstate_ncsb* rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);
  assert(rhs_ncsb);
  return (this->active_ == rhs_ncsb->active_) &&
    (this->check_ == rhs_ncsb->check_) &&
    (this->safe_ == rhs_ncsb->safe_) &&
    (this->breakpoint_ == rhs_ncsb->breakpoint_);
}


bool mstate_ncsb::lt(const mstate& rhs) const
{ // {{{
  const mstate_ncsb* rhs_ncsb = dynamic_cast<const mstate_ncsb*>(&rhs);
  assert(rhs_ncsb);

  if (this->active_ != rhs_ncsb->active_) { return this->active_ < rhs_ncsb->active_; }
  if (this->check_ != rhs_ncsb->check_) { return this->check_ < rhs_ncsb->check_; }
  if (this->safe_ != rhs_ncsb->safe_) { return this->safe_ < rhs_ncsb->safe_; }
  if (this->breakpoint_ != rhs_ncsb->breakpoint_) { return this->breakpoint_ < rhs_ncsb->breakpoint_; }

  return false;   // if all are equal
} // lt() }}}

/// implementation of NCSB-based complementation algorithm for deterministic SCCs
class complement_ncsb_delay : public complement_ncsb 
{ // {{{
private:
    struct cmp{
        bool operator()(const mstate_ncsb *const &a,
                         const mstate_ncsb *const &b) const
        {
            return a->lt(*b); 
        }
    };
    bool comparator = [](const mstate_ncsb *const &a,
                         const mstate_ncsb *const &b)
    { return a->lt(*b); };
    std::set<const mstate_ncsb*, cmp> active_mstates_;
    std::map<const mstate_ncsb *, std::set<const mstate_ncsb *, cmp>, cmp> successors_;
    std::map<mstate_ncsb, std::set<mstate_ncsb>> succ_ncsb_;

public: // METHODS

  /// constructor
  complement_ncsb_delay(const cmpl_info& info, unsigned part_index);

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

  virtual ~complement_ncsb_delay() override;

  bool closes_a_cycle(const mstate_ncsb *src, const mstate_ncsb *dst);
}; // complement_ncsb }}}
} // namespace kofola }}}

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// The acceptance condition of the IADACs determinization: one disjunct per run
/// label.
///
/// The input condition is first compacted onto the colours it mentions (colours
/// of the input the formula does not mention do not influence acceptance and
/// are dropped, the remaining ones are renumbered from 0).  With k colours left,
/// the label i owns the block [(k+1)*i, (k+1)*(i+1)): the k compacted colours
/// followed by its discontinuation colour, and its disjunct is the compacted
/// condition over that block conjoined with Fin(discontinuation colour).  The
/// colours thus start at 0, so there is nothing to rebase (get_min_colour() of
/// the algorithms using this class is 0).
class determinisation_acc_cond 
{ // {{{
private:
    spot::acc_cond::acc_code acc_code_;
    std::vector<unsigned> additional_fins_;
    unsigned disj_size_;
    /// the colours of the input the condition mentions
    spot::acc_cond::mark_t used_;
    
public:
    /// throws std::runtime_error if colours_needed() exceeds SPOT_MAX_ACCSETS
    determinisation_acc_cond(const spot::acc_cond::acc_code& acc_cond, unsigned disjuncts);

    /// the number of colours of the condition built for 'acc_cond' and 'disjuncts'
    static unsigned colours_needed(const spot::acc_cond::acc_code& acc_cond, unsigned disjuncts);

    spot::acc_cond get_acc_cond() const;

    spot::acc_cond::mark_t get_all_discontinuation_colours() const;

    unsigned get_fin_mark(unsigned disjunct_index) const {
        assert(disjunct_index < additional_fins_.size());
        return additional_fins_[disjunct_index];
    }

    /// inserts into 'out' the colours of the label 'disjunct_index' that
    /// correspond to the input colours 'acc' (the colours the condition does not
    /// mention are dropped)
    void add_colours(spot::acc_cond::mark_t acc, unsigned disjunct_index, std::set<unsigned>& out) const {
        assert(disjunct_index < additional_fins_.size());
        for (unsigned col : acc.strip(~used_).sets()) {
            out.insert(col + disj_size_ * disjunct_index);
        }
    }
};

class complement_init_almost_det : public abstract_complement_alg
{ // {{{
public: // METHODS

  /// constructor
  complement_init_almost_det(const cmpl_info& info, unsigned part_index, unsigned runs_bound);

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

  // no breakpoint at all
  virtual bool use_shared_breakpoint() const override { return false; }

  virtual spot::acc_cond get_acc_cond() override
  { return acc_cond_.get_acc_cond().get_acceptance().complement(); }

  virtual unsigned get_min_colour() const override { return 0; }

  virtual ~complement_init_almost_det() override;

private:
    unsigned runs_bound_;
    determinisation_acc_cond acc_cond_;
}; // complement_init_almost_det }}}
} // namespace kofola }}}


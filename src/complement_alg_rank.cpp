// implementation of rank-based complementation from Sven Schewe's paper

#include "complement_alg_rank.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

namespace { // {{{

/// representation of all other runs (outside the partition block)
const int BOX = -1;

class ranking : public std::map<int, int>
{
private:
  unsigned max_rank = 0;

public:
  ranking() : std::map<int, int>(){};
  std::string get_name();
  unsigned get_max_rank(){return max_rank;};
  void set_max_rank(unsigned max_rank){this->max_rank = max_rank;};
  void check_tight(std::vector<ranking> rankings);
  bool is_bigger(ranking other);

  bool operator==(const ranking &other) const
  {
    if (this->max_rank != other.max_rank) { return false; }
    for (auto it=this->begin(); it!=this->end(); it++) {
      if (it->second != other.at(it->first)) { return false; }
    }
    return true;
  }

  bool operator<(const ranking &other) const
  {
    if (this->max_rank == other.max_rank) {
        for (auto it=this->begin(); it!=this->end(); it++) {
          if (it->second != other.at(it->first)) {
            return it->second < other.at(it->first);
          }
        }
        return false;
    } else {
      return this->max_rank < other.max_rank;
    }
  }
};


/// partial macrostate for the given component
class mstate_rank : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  bool active_;
  std::set<unsigned> states_;
  std::set<unsigned> breakpoint_;
  ranking f_;
  int i_;

public: // METHODS

  /// constructor
  mstate_rank(
    const std::set<unsigned>&  states,
    const std::set<unsigned>&  breakpoint,
    const ranking&             f,
    int                        i,
    bool                       active
  ) : states_(states),
    breakpoint_(breakpoint),
    f_(f),
    i_(i),
    active_(active)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return this->active_; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_rank() override { }

  friend class kofola::complement_rank;
}; // mstate_rank }}}

std::string mstate_rank::to_string() const
{
  assert(false);
}

bool mstate_rank::eq(const mstate& rhs) const
{
  assert(false);
}

bool mstate_rank::lt(const mstate& rhs) const
{
  assert(false);
}


} // anonymous namespace }}}


complement_rank::complement_rank(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index)
{ }


mstate_set complement_rank::get_init() const
{
  DEBUG_PRINT_LN("init RANK for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_state.insert(orig_init);
  }

  // std::shared_ptr<mstate> ms(new mstate_rank(init_state, {}, {}, false));
  // mstate_set result = {ms};
  // return result;
  //
  // complement_mstate mstate(scc_info_);
  //
  // unsigned orig_init = aut_->get_init_state_number();
  // mstate.curr_reachable_.push_back(orig_init);
  //
  // rank_state tmp;
  // tmp.reachable.insert(BOX);
  // mstate.na_sccs_.push_back(tmp);
  //
  // return mstate;

  assert(false);
}

mstate_col_set complement_rank::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}

mstate_set complement_rank::lift_track_to_active(const mstate* src) const
{
  assert(false);
}

mstate_col_set complement_rank::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}


complement_rank::~complement_rank()
{ }

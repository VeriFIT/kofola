// implementation of rank-based complementation from Sven Schewe's paper

#include "complement_alg_rank.hpp"
#include "rankings.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

namespace { // {{{

/// representation of all other runs (outside the partition block)
const unsigned BOX = UINT_MAX;

// class ranking : public std::map<int, int>
// {
// private:
//   unsigned max_rank = 0;
//
// public:
//   ranking() : std::map<int, int>(){};
//   std::string get_name();
//   unsigned get_max_rank(){return max_rank;};
//   void set_max_rank(unsigned max_rank){this->max_rank = max_rank;};
//   void check_tight(std::vector<ranking> rankings);
//   bool is_bigger(ranking other);
//
//   bool operator==(const ranking &other) const
//   {
//     if (this->max_rank != other.max_rank) { return false; }
//     for (auto it=this->begin(); it!=this->end(); it++) {
//       if (it->second != other.at(it->first)) { return false; }
//     }
//     return true;
//   }
//
//   bool operator<(const ranking &other) const
//   {
//     if (this->max_rank == other.max_rank) {
//         for (auto it=this->begin(); it!=this->end(); it++) {
//           if (it->second != other.at(it->first)) {
//             return it->second < other.at(it->first);
//           }
//         }
//         return false;
//     } else {
//       return this->max_rank < other.max_rank;
//     }
//   }
// };


/// partial macrostate for the given component
class mstate_rank : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  bool active_;
  std::set<unsigned> states_;
  bool is_waiting_;                 // true - waiting, false - tight
  std::set<unsigned> breakpoint_;
  ranking f_;
  int i_;

public: // METHODS

  /// constructor
  mstate_rank(
    const std::set<unsigned>&  states,
    bool                       is_waiting,
    const std::set<unsigned>&  breakpoint,
    const ranking&           f,
    int                        i,
    bool                       active
  ) : states_(states),
    is_waiting_(is_waiting),
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

  friend std::set<unsigned> get_successors_with_box(
    const std::set<unsigned>&  glob_reach,
    const mstate_rank&         rank_state,
    unsigned                   part_index,
    const cmpl_info&           info);
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


/// returns true iff 'pred' is a predecessor of 'state'
bool is_predecessor_of(unsigned pred, unsigned state, const cmpl_info& info)
{
  return kofola::is_in(state, info.reachable_vector_[pred]);
}


std::set<unsigned> get_successors_with_box(
  const std::set<unsigned>&  glob_reach,
  const mstate_rank&         rank_state,
  unsigned                   part_index,
  const cmpl_info&           info)
{ // {{{
  std::set<unsigned> succ;

  std::set<unsigned> state_set;
  if (rank_state.is_waiting_) {
    state_set = rank_state.states_;
  } else {
    for (auto pr : rank_state.f_) {
      state_set.insert(pr.first);
    }
  }

  if (state_set.empty()) { return {}; }

  for (auto state : glob_reach) { // collect all reached states from part block
    if (info.st_to_part_map_.at(state) == part_index) {
      succ.insert(state);
    }

    if (kofola::is_in(BOX, state_set)) { // BOX is present
      for (auto state : glob_reach) { // if the reached state not from this
        // partition can still reach this partition, add BOX
        if (info.st_to_part_map_.at(state) != part_index) {
          unsigned state_scc_index = info.scc_info_.scc_of(state);
          for (unsigned part_scc_index : info.part_to_scc_map_.at(part_index)) {
            if (kofola::is_in(state_scc_index,
                info.scc_to_pred_sccs_map_.at(part_scc_index))) {
              succ.insert(BOX);
              break;
            }
          }
        }
      }
    }
  }

  return succ;
} // get_successors_with_box() }}}


/// retrieves the WAITING part of an automaton
static waiting get_waiting_part(const spot::const_twa_graph_ptr& aut)
{ // {{{
  std::set<std::set<unsigned>> states;
  std::map<std::set<unsigned>, std::set<std::set<unsigned>>> trans;

  // initial state
  unsigned init = aut->get_init_state_number();
  std::set<unsigned> init_mstate;
  init_mstate.insert(init);

  std::stack<std::set<unsigned>> stack;
  stack.push(init_mstate);
  states.insert(init_mstate);

  while (not stack.empty())
  {
    auto state = stack.top();
    stack.pop();
    std::set<std::set<unsigned>> succ;

    // outgoing symbols
    std::vector<bdd> alphabet;

    bdd msupport = bddtrue;
    bdd n_s_compat = bddfalse;

    // const std::set<unsigned> &reach_set = state;
    // compute the occurred variables in the outgoing transitions of ms, stored in msupport
    for (unsigned s : state)
    {
      for (const auto &out : aut->out(s)) { // FIXME: this should be pre-computed for every state
        msupport &= bdd_support(out.cond);
        n_s_compat |= out.cond;
      }
    }

    bdd all = n_s_compat;

    while (all != bddfalse)
    {
      bdd letter = bdd_satoneset(all, msupport, bddfalse);
      all -= letter;
      alphabet.push_back(letter);
    }

    std::vector<std::set<unsigned>> alphabet_map;
    for (unsigned i=0; i<alphabet.size(); i++)
    {
      alphabet_map.push_back(std::set<unsigned>());
    }

    for (unsigned s : state) // every state in a macrostate
    {
      for (const auto &t : aut->out(s))
      {
        for (unsigned i=0; i<alphabet.size(); i++)
        {
          if (bdd_implies(alphabet[i], t.cond))
          {
            alphabet_map[i].insert(t.dst);
          }
        }
      }
    }

    for (auto item : alphabet_map)
    {
      if (states.find(item) == states.end())
      {
        states.insert(item);
        stack.push(item);
      }
      succ.insert(item);
    }

    trans.insert({state, succ});
  }

  return waiting(states, trans);
} // get_waiting_part() }}}


/// returns maximum rankings from a set of rankings
std::vector<ranking> get_max_rankings(const std::vector<ranking>& rankings)
{
  std::vector<ranking> tmp(rankings.begin(), rankings.end());

  for (auto r : rankings)
  {
      unsigned max_rank = r.get_max_rank();
      if (std::any_of(tmp.begin(), tmp.end(), [max_rank, r](ranking r2)
                      { return r != r2 and max_rank == r2.get_max_rank() and r2.is_bigger(r); }))
      {
          // there is bigger ranking
          tmp.erase(std::remove(tmp.begin(), tmp.end(), r), tmp.end());
      }
  }

  return tmp;
}

} // anonymous namespace }}}


complement_rank::complement_rank(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index),
  waiting_(get_waiting_part(info.aut_))
{ // {{{
  // compute rank restrictions
  unsigned states_in_part = 0;
  for (unsigned scc_index : info.part_to_scc_map_.at(part_index)) {
    // TODO: maybe we can consider only SCCs?
    info.scc_info_.states_of(scc_index).size();
  }

  for (const auto& mst : this->waiting_.get_states()) {
    // initialization
    unsigned nonacc = 0;
    for (auto state : mst) {
      if (!info.state_accepting_[state] &&
        info_.st_to_part_map_.at(state) == part_index)
      {
        nonacc++;
      }
    }
    rank_restr_.insert({mst, 2*(nonacc + 1)});
  }

  DEBUG_PRINT_LN("waiting_ for partition " + std::to_string(part_index) + ": " +
    std::to_string(this->waiting_));
  DEBUG_PRINT_LN("rank_restr_ for partition " + std::to_string(part_index) + ": " +
    std::to_string(this->rank_restr_));
} // complement_rank() }}}


mstate_set complement_rank::get_init() const
{ // {{{
  DEBUG_PRINT_LN("init RANK for partition " + std::to_string(this->part_index_));
  std::set<unsigned> init_state;

  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_state.insert(orig_init);
  } else { // TODO: we assert that the partition block is reachable from initial
    // FIXME: the box should not be necessary in cases such as a NAC reachable
    // from an initial state via an acyclic path
    init_state.insert(BOX);
  }

  std::shared_ptr<mstate> ms(new mstate_rank(
    init_state,        // reachable states (S)
    true,              // is it Waiting?
    {},                // breakpoint (O)
    {},                // ranking (f)
    -1,                // index of tracked rank (i)
    false));           // active
  mstate_set result = {ms};
  return result;
} // get_init() }}}


mstate_col_set complement_rank::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{
  assert(false);
}


mstate_set complement_rank::lift_track_to_active(const mstate* src) const
{ // {{{
  const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
  assert(src_rank);
  assert(!src_rank->active_);

  mstate_set result;
  if (src_rank->is_waiting_) { // src is from WAITING
    std::shared_ptr<mstate> src_cpy(new mstate_rank(*src_rank));
    result.push_back(src_cpy);          // one option is to stay in WAITING

    // and let's compute the successors that move to TIGHT
    std::vector<std::tuple<int, int, bool>> r;
    // FIXME is this correct? to consider only states from this partition
    // auto bound = rank_restr_[std::set<unsigned>(mstate.curr_reachable_.begin(), mstate.curr_reachable_.end())];
    auto bound = rank_restr_.at(src_rank->states_);

    for (auto s : src_rank->states_) {
      bool accepting = ((s != BOX)? this->info_.state_accepting_[s] : false);
      r.push_back(std::make_tuple(s, bound, accepting));
    }

    std::vector<ranking> rankings = get_tight_rankings(r);
    rankings = get_max_rankings(rankings);
    for (auto rnking : rankings)
    {
      std::set<unsigned> breakpoint;
      for (auto pr : rnking) { // construct breakpoint
        if (pr.second == 0) { breakpoint.insert(pr.first); }
      }
      std::shared_ptr<mstate> ms(new mstate_rank(
        {},                // reachable states (S)
        false,             // is it Waiting?
        breakpoint,        // breakpoint (O)
        rnking,            // ranking (f)
        0,                 // index of tracked rank (i)
        true));            // active
      result.push_back(ms);
    }
  } else { // src is from TIGHT
    assert(false);
  }

  DEBUG_PRINT_LN("complement_rank::lift returning " + std::to_string(result));

  return result;
} // lift_track_to_active() }}}


mstate_col_set complement_rank::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
  assert(src_rank);
  assert(src_rank->active_);

  assert(false);
#if 0
  mstate_col_set result;

  if (src_rank->is_waiting_)
  { // if the source is from the waiting part
    if (src_rank->states_.size() == 0 ||
        (src_rank->states_.size() == 1 && kofola::is_in(BOX, src_rank->states_))) {
      // in case src does not track any state from the partition block
      rank_state new_state;
      new_state.is_waiting_ = true;
      new_state.reachable = get_successors_with_box(glob_reached, src_rank, symbol, scc_index_[0]);
      complement_mstate tmp(scc_info_);
      tmp.na_sccs_.push_back(new_state);
      result.push_back({tmp, true});
    } else {
      result = get_succ_track_to_active(mstate, symbol);
    }
  }

  else
  { // tight part
      std::vector<std::pair<complement_mstate, bool>> ret = get_succ_track(mstate, symbol);
      if (ret.size() > 0)
      {
          ranking g = ret[0].first.na_sccs_[0].f;

          std::vector<rank_state> eta_3;
          std::vector<rank_state> eta_4;

          // eta 3
          if (src_rank->O.size() > 0)
          {
              rank_state new_state;
              new_state.is_waiting_ = false;
              new_state.f = g;
              new_state.i = src_rank->i;

              rank_state tmp;
              tmp.reachable = src_rank->O;
              std::set<int> O_succ = get_successors_with_box(glob_reached, tmp, symbol, scc_index_[0]);
              std::set<int> g_rev;
              for (auto pr : g)
              {
                  if (pr.second == src_rank->i)
                      g_rev.insert(pr.first);
              }
              std::set_intersection(O_succ.begin(), O_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(new_state.O, new_state.O.begin()));

              eta_3.push_back(new_state);
          }
          else
          {
              rank_state new_state;
              new_state.is_waiting_ = false;
              new_state.f = g;
              new_state.i = (src_rank->i + 2) % (g.get_max_rank() + 1);

              std::set<int> dom_succ = get_successors_with_box(glob_reached, src_rank, symbol, scc_index_[0]);
              std::set<int> g_rev;
              for (auto pr : g)
              {
                  if (pr.second == new_state.i)
                      g_rev.insert(pr.first);
              }
              std::set_intersection(dom_succ.begin(), dom_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(new_state.O, new_state.O.begin()));

              eta_3.push_back(new_state);
          }

          // eta 4
          if (src_rank->i != src_rank->f.get_max_rank() - 1)
          {
              rank_state new_state;

              std::set<int> M;
              rank_state tmp;
              tmp.reachable = src_rank->O;
              std::set<int> O_succ = get_successors_with_box(glob_reached, tmp, symbol, scc_index_[0]);
              std::set<int> g_rev;
              for (auto pr : g)
              {
                  if (pr.second == src_rank->i)
                      g_rev.insert(pr.first);
              }
              std::set_intersection(O_succ.begin(), O_succ.end(), g_rev.begin(), g_rev.end(), std::inserter(M, M.begin()));

              new_state.is_waiting_ = false;
              new_state.i = src_rank->i;
              for (auto s : M)
              {
                  if (s != BOX and is_accepting_[s])
                      new_state.O.insert(s);
              }

              ranking g_prime = g;
              for (auto &pr : g_prime)
              {
                  if (pr.first != BOX and (not is_accepting_[pr.first] and M.find(pr.first) != M.end()))
                      pr.second = pr.second - 1;
              }
              new_state.f = g_prime;

              eta_4.push_back(new_state);
          }

          std::vector<rank_state> U;
          if (eta_3.size() > 0 and eta_3[0].O.size() == 0 and eta_3[0].i == eta_3[0].f.get_max_rank() - 1)
              U.push_back(eta_3[0]);
          if (eta_4.size() > 0 and eta_4[0].O.size() == 0 and eta_4[0].i == eta_4[0].f.get_max_rank() - 1)
              U.push_back(eta_4[0]);

          if (eta_3.size() > 0 and std::find(U.begin(), U.end(), eta_3[0]) == U.end())
          {
              complement_mstate tmp_mstate(scc_info_);
              tmp_mstate.na_sccs_.push_back(eta_3[0]);
              result.push_back({tmp_mstate, false});
          }
          if (eta_4.size() > 0 and std::find(U.begin(), U.end(), eta_4[0]) == U.end())
          {
              complement_mstate tmp_mstate(scc_info_);
              tmp_mstate.na_sccs_.push_back(eta_4[0]);
              result.push_back({tmp_mstate, false});
          }

          for (rank_state s : U)
          {
              if (not one_scc)
              {
                  rank_state tmp;
                  tmp.is_waiting_ = false;
                  tmp.f = s.f;
                  complement_mstate tmp_mstate(scc_info_);
                  tmp_mstate.na_sccs_.push_back(tmp);
                  result.push_back({tmp_mstate, true});
              }
              else
              {
                  complement_mstate tmp(scc_info_);
                  tmp.na_sccs_.push_back(s);
                  result.push_back({tmp, true});
              }
          }
      }
  }

  return result;
#endif
} // get_succ_active() }}}


complement_rank::~complement_rank()
{ }

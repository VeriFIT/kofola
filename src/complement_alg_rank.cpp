// implementation of rank-based complementation from Sven Schewe's paper

#include "complement_alg_rank.hpp"
#include "rankings.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;
using RankRestriction = complement_rank::RankRestriction;

namespace { // {{{

/// representation of all other runs (outside the partition block)
const unsigned BOX = UINT_MAX;


/// partial macrostate for the given component
class mstate_rank : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  bool active_;
  std::set<unsigned> states_;       // {} when !is_waiting_
  bool is_waiting_;                 // true - waiting, false - tight
  std::set<unsigned> breakpoint_;   // {} when !active_ || is_waiting_
  ranking f_;                       // {} when is_waiting_
  int i_;                           // -1 when !active_ || is_waiting_

public: // METHODS

  /// constructor
  mstate_rank(
    const std::set<unsigned>&  states,
    bool                       is_waiting,
    const std::set<unsigned>&  breakpoint,
    const ranking&             f,
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

  bool invariants_hold() const;

  friend class kofola::complement_rank;

  friend std::set<unsigned> get_successors_with_box(
    const std::set<unsigned>&  glob_reach,
    const mstate_rank&         rank_state,
    unsigned                   part_index,
    const cmpl_info&           info);

  friend std::vector<ranking> get_maxrank(
    const std::set<unsigned>&  glob_reach,
    const mstate_rank&         rank_state,
    const RankRestriction&     rank_restr,
    unsigned                   part_index,
    const bdd&                 symbol,
    const cmpl_info&           info);
}; // mstate_rank }}}

bool mstate_rank::invariants_hold() const
{ // {{{
  if (!this->is_waiting_ && !this->states_.empty()) {
    return false;
  }
  if ((!this->active_ || this->is_waiting_ ) && !this->breakpoint_.empty()) {
    return false;
  }
  if (this->is_waiting_ && !this->f_.empty()) {
    return false;
  }
  if ((!this->active_ || this->is_waiting_ ) && this->i_ != -1) {
    return false;
  }

  return true;
} // invariants_hold() }}}

std::string mstate_rank::to_string() const
{ // {{{
  std::string res = std::string("[RANK(") + ((this->active_)? "A" : "T") + "): ";
  res += "S=";
  if (this->is_waiting_) { // WAITING
    res += std::to_string(this->states_);
  } else { // TIGHT
    std::set<int> dom_of_f;
    for (const auto& pair : this->f_) {
      dom_of_f.insert(pair.first);
    }
    res += std::to_string(dom_of_f);
    if (this->active_) {
      res += ", O=" + std::to_string(this->breakpoint_);
    }
    res += ", f=" + std::to_string(this->f_);
    if (this->active_) {
      res += ", i=" + std::to_string(this->i_);
    }
  }
  res += "]";
  return res;
} // to_string() }}}


bool mstate_rank::eq(const mstate& rhs) const
{
  const mstate_rank* rhs_rank = dynamic_cast<const mstate_rank*>(&rhs);
  assert(rhs_rank);

  return this->active_ == rhs_rank->active_ &&
    this->is_waiting_ == rhs_rank->is_waiting_ &&
    this->states_ == rhs_rank->states_ &&
    this->breakpoint_ == rhs_rank->breakpoint_ &&
    this->f_ == rhs_rank->f_ &&
    this->i_ == rhs_rank->i_;
}

bool mstate_rank::lt(const mstate& rhs) const
{
  const mstate_rank* rhs_rank = dynamic_cast<const mstate_rank*>(&rhs);
  assert(rhs_rank);

  if (this->active_ == rhs_rank->active_) {
    if (this->is_waiting_ == rhs_rank->is_waiting_) {
      if (this->i_ == rhs_rank->i_) {
        if (this->states_ == rhs_rank->states_) {
          if (this->breakpoint_ == rhs_rank->breakpoint_) {
            return this->f_ < rhs_rank->f_;
          } else {
            return this->breakpoint_ < rhs_rank->breakpoint_;
          }
        } else {
          return this->states_ < rhs_rank->states_;
        }
      } else {
        return this->i_ < rhs_rank->i_;
      }
    } else {
      return this->is_waiting_ < rhs_rank->is_waiting_;
    }
  } else {
    return this->active_ < rhs_rank->active_;
  }
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
std::vector<ranking> get_max(const std::vector<ranking>& rankings)
{ // {{{
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
} // get_max() }}}


bool compare_ranks(
  const std::tuple<int, int, bool>& first,
  const std::tuple<int, int, bool>& second)
{
  return (std::get<1>(first) < std::get<1>(second));
}


std::vector<ranking> cart_product(
  std::vector<ranking> rankings,
  std::tuple<int, int, bool> state)
{ // {{{
  int state_name = std::get<0>(state);
  int max_rank = std::get<1>(state);
  bool accepting = std::get<2>(state);

  std::vector<ranking> result;

  for (int i=0; i<max_rank; (accepting ? i+=2 : i++)) {
    std::vector<ranking> new_rankings(rankings.begin(), rankings.end());
    for (auto r : new_rankings) {
      r[state_name] = i;
      if (r.get_max_rank() < i) {
        r.set_max_rank(i);
      }
      result.push_back(r);
    }
  }

  return result;
} // cart_product() }}}


std::set<int> get_all_successors_acc(
  const spot::const_twa_graph_ptr&  aut,
  const spot::scc_info&             scc_info,
  const std::set<unsigned>&         current_states,
  const bdd&                        symbol,
  unsigned                          part_index)
{ // {{{
  std::set<int> successors;
  spot::acc_cond::mark_t acc = {0};

  for (unsigned s : current_states) {
    for (const auto &t : aut->out(s)) {
      if (!bdd_implies(symbol, t.cond)) { continue; }

      if (t.acc == acc && scc_info.scc_of(t.dst) == part_index) {
        successors.insert((int)t.dst);
      }
    }
  }

  return successors;
} // get_all_successors_acc() }}}


void check_tight(std::vector<ranking>& rankings)
{ // {{{
  std::vector<ranking> result;
  for (auto r : rankings) {
    bool skip = false;
    if (r.get_max_rank() % 2 == 1) {
      // check box
      int max_rank = -1;
      if (r.find(BOX) != r.end()) {
        if (r.get_max_rank() != r[BOX]) {
          skip = true;
        }
        max_rank = r[BOX];
      }

      if (not skip) {
        std::set<int> ranks;
        for (auto pr : r) {
            if (pr.first != BOX and pr.second == max_rank) {
              skip = true;
              break;
            }
            ranks.insert(pr.second);
        }

        if (not skip) {
          for (int i=1; i<r.get_max_rank(); i+=2) {
            if (ranks.find(i) == ranks.end()) {
              skip = true;
              break;
            }
          }

          if (not skip) {
            result.push_back(r);
          }
        }
      }
    }
  }

  rankings = result;
} // check_tight() }}}


std::vector<ranking> get_tight_rankings(
  const std::vector<std::tuple<int, int, bool>>& mp)
{ // {{{
  std::vector<ranking> rankings;

  // get max rank
  if (mp.size() > 0)
  {
    auto max = std::max_element(mp.begin(), mp.end(), compare_ranks);
    int max_rank = std::get<1>(*max);
    if (max_rank > 2*mp.size()) {
      max_rank = 2*mp.size();
    }

    // odd rank
    if (max_rank == 0) { return rankings; }
    if (max_rank % 2 == 0) {
      max_rank--;
    }

    for (const auto& state : mp) {
      if (rankings.size() == 0) {
        for (int i=0; i<std::get<1>(state); (std::get<2>(state) ? i+=2 : i++)) {
          ranking r;
          r[std::get<0>(state)] = i;
          r.set_max_rank(i);
          rankings.push_back(r);
        }
      } else {
        rankings = cart_product(rankings, state);
      }
    }

    check_tight(rankings);
  }

  return rankings;
} // get_tight_rankings() }}}


std::vector<ranking> get_succ_rankings(
  const ranking&                                  r,
  const std::vector<std::tuple<int, int, bool>>&  restr,
  const std::set<unsigned>&                       glob_reached,
  const bdd&                                      symbol,
  unsigned                                        part_index,
  const cmpl_info&                                info)
{ // {{{
  std::vector<ranking> rankings = get_tight_rankings(restr);
  std::set<ranking> rankings_set(rankings.begin(), rankings.end());

  DEBUG_PRINT_LN("tight rankings: " + std::to_string(rankings))

  for (ranking r2 : rankings) {
    if (r2.get_max_rank() != r.get_max_rank()) {
      rankings_set.erase(r2);
      continue;
    }

    bool skip = false;
    for (auto pr : r) {
      unsigned state = pr.first;

      mstate_rank tmp(
        {state},       // reachable states (S)
        true,          // is it Waiting?
        {},            // breakpoint (O)
        {},            // ranking (f)
        -1,            // index of tracked rank (i)
        false);        // active

      std::set<unsigned> succ = get_successors_with_box(
        glob_reached, tmp, part_index, info);

      for (auto s : succ) {
        if (r2.at(s) > r.at(state)) {
          skip = true;
          rankings_set.erase(r2);
          break;
        }
      }

      if (skip) { break; }

      if (state != BOX) {
        std::set<int> succ = get_all_successors_acc(
          info.aut_, info.scc_info_, {state}, symbol, part_index);

        unsigned rank = (r.at(state) % 2 == 0 ? r.at(state) : r.at(state) - 1);
        for (auto s : succ) {
          if (r2.at(s) > rank) {
            skip = true;
            rankings_set.erase(r2);
            break;
          }
        }

        if (skip) { break; }
      }
    }
  }

  return std::vector<ranking>(rankings_set.begin(), rankings_set.end());
} // get_succ_rankings() }}}


std::vector<ranking> get_maxrank(
  const std::set<unsigned>&  glob_reached,
  const mstate_rank&         rank_state,
  const RankRestriction&     rank_restr,
  unsigned                   part_index,
  const bdd&                 symbol,
  const cmpl_info&           info)
{ // {{{
  std::vector<std::tuple<int, int, bool>> restr;

  std::set<unsigned> domain;
  for (auto pr : rank_state.f_) {
    domain.insert(pr.first);
  }

  auto succ_domain = get_successors_with_box(
    glob_reached, rank_state, part_index, info);

  auto bound = rank_restr.at(domain);
  for (auto s : succ_domain) {
    restr.push_back({s, bound, (s == BOX) ? false : info.state_accepting_[s]});
  }

  std::vector<ranking> succ_rankings = get_succ_rankings(
    rank_state.f_, restr, glob_reached, symbol, part_index, info);

  succ_rankings = get_max(succ_rankings);
  DEBUG_PRINT_LN("obtained succ rankings = " + std::to_string(succ_rankings));

  return succ_rankings;
} // get_maxrank() }}}


} // anonymous namespace }}}


complement_rank::complement_rank(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index),
  waiting_(get_waiting_part(info.aut_))
{ // {{{
  // compute rank restrictions
  unsigned states_in_part = 0;
  for (unsigned scc_index : this->info_.part_to_scc_map_.at(part_index)) {
    // TODO: maybe we can consider only SCCs?
    this->info_.scc_info_.states_of(scc_index).size();
  }

  for (const auto& mst : this->waiting_.get_states()) {
    // initialization
    unsigned nonacc = 0;
    for (auto state : mst) {
      if (!this->info_.state_accepting_[state] &&
        this->info_.st_to_part_map_.at(state) == part_index)
      {
        nonacc++;
      }
    }
    this->rank_restr_.insert({mst, 2*(nonacc + 1)});
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
{ // {{{
  const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
  assert(src_rank);
  assert(!src_rank->active_);
  assert(src_rank->invariants_hold());

  if (src_rank->is_waiting_) { // WAITING
    std::set<unsigned> succs = get_successors_with_box(glob_reached, *src_rank,
      this->part_index_, this->info_);

    std::shared_ptr<mstate> ms(new mstate_rank(
      succs,             // reachable states (S)
      true,              // is it Waiting?
      {},                // breakpoint (O)
      {},                // ranking (f)
      -1,                // index of tracked rank (i)
      false));           // active
    mstate_col_set result = {{ms, {}}};
    return result;
  } else { // TIGHT
    std::vector<ranking> maxrank = get_maxrank(glob_reached, *src_rank,
       this->rank_restr_, this->part_index_, symbol, this->info_);

    if (maxrank.size() == 0) { return {}; }
    assert(maxrank.size() == 1);

    std::shared_ptr<mstate> ms(new mstate_rank(
      {},                // reachable states (S)
      false,             // is it Waiting?
      {},                // breakpoint (O)
      maxrank[0],        // ranking (f)
      -1,                // index of tracked rank (i)
      false));           // active
    mstate_col_set result = {{ms, {}}};
    return result;
  }
} // get_succ_track() }}}


mstate_set complement_rank::lift_track_to_active(const mstate* src) const
{ // {{{
  const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
  assert(src_rank);
  assert(!src_rank->active_);
  assert(src_rank->invariants_hold());

  mstate_set result;
  if (src_rank->is_waiting_) { // src is from WAITING
    std::shared_ptr<mstate_rank> src_cpy(new mstate_rank(*src_rank));
    src_cpy->active_ = true;
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
    rankings = get_max(rankings);
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
  assert(src_rank->invariants_hold());

  DEBUG_PRINT_LN("tracking successor of: " + std::to_string(*src_rank));
  mstate_rank tmp(
    src_rank->states_,       // reachable states (S)
    src_rank->is_waiting_,   // is it Waiting?
    {},                      // breakpoint (O)
    src_rank->f_,            // ranking (f)
    -1,                      // index of tracked rank (i)
    false);                  // active
  mstate_col_set track_succ = this->get_succ_track(glob_reached, &tmp, symbol);

  if (track_succ.size() == 0) { return {};}
  assert(track_succ.size() == 1);

  const mstate_rank* track_ms = dynamic_cast<const mstate_rank*>(
    track_succ[0].first.get());
  assert(track_ms);

  DEBUG_PRINT_LN("obtained track ms: " + std::to_string(*track_ms));

  if (src_rank->is_waiting_) { // WAITING
    if (src_rank->states_.size() == 0 ||
        (src_rank->states_.size() == 1 && kofola::is_in(BOX, src_rank->states_))) {
      // in case src does not track any state from the partition block

      std::shared_ptr<mstate> ms(new mstate_rank(*track_ms));
      return {{ms, {0}}};
    } else {
      mstate_set lifted = this->lift_track_to_active(track_ms);
      mstate_col_set result;
      for (const auto& ms : lifted) {
        result.push_back({ms, {}});
      }
      return result;
    }
  } else { // TIGHT
    ranking g = track_ms->f_;

    std::vector<mstate_rank> eta_3;
    std::vector<mstate_rank> eta_4;

    // eta 3
    if (src_rank->breakpoint_.size() > 0) { // breakpoint is nonempty
      // TODO: would it make sense to detect emptiness of breakpoint after
      // making the transition?
      mstate_rank tmp(
        src_rank->breakpoint_,   // reachable states (S)
        true,                    // is it Waiting?
        {},                      // breakpoint (O)
        {},                      // ranking (f)
        -1,                      // index of tracked rank (i)
        false);                  // active
      std::set<unsigned> O_succ = get_successors_with_box(glob_reached, tmp,
        this->part_index_, this->info_);
      std::set<unsigned> g_rev;
      for (auto pr : g) {
        if (pr.second == src_rank->i_) {
          g_rev.insert(pr.first);
        }
      }
      std::set<unsigned> new_breakpoint = kofola::get_set_intersection(O_succ, g_rev);

      mstate_rank new_state(
        {},               // reachable states (S)
        false,            // is it Waiting?
        new_breakpoint,   // breakpoint (O)
        g,                // ranking (f)
        src_rank->i_,     // index of tracked rank (i)
        true);            // active

      eta_3.push_back(new_state);
    } else { // breakpoint is empty
      int new_i = (src_rank->i_ + 2) % (g.get_max_rank() + 1);
      std::set<unsigned> dom_succ = get_successors_with_box(glob_reached, *src_rank,
        this->part_index_, this->info_);
      std::set<unsigned> g_rev;
      for (auto pr : g) {
        if (pr.second == new_i) {
          g_rev.insert(pr.first);
        }
      }

      std::set<unsigned> new_breakpoint =
        kofola::get_set_intersection(dom_succ, g_rev);

      mstate_rank new_state(
        {},               // reachable states (S)
        false,            // is it Waiting?
        new_breakpoint,   // breakpoint (O)
        g,                // ranking (f)
        new_i,            // index of tracked rank (i)
        true);            // active

      eta_3.push_back(new_state);
    }

    // eta 4
    if (src_rank->i_ != src_rank->f_.get_max_rank() - 1) {
      mstate_rank tmp(
        src_rank->breakpoint_,   // reachable states (S)
        true,                    // is it Waiting?
        {},                      // breakpoint (O)
        {},                      // ranking (f)
        -1,                      // index of tracked rank (i)
        false);                  // active
      std::set<unsigned> O_succ = get_successors_with_box(glob_reached, tmp,
        this->part_index_, this->info_);
      std::set<unsigned> g_rev;
      for (auto pr : g) {
        if (pr.second == src_rank->i_) {
          g_rev.insert(pr.first);
        }
      }
      std::set<unsigned> M = kofola::get_set_intersection(O_succ, g_rev);

      std::set<unsigned> new_breakpoint;
      for (auto s : M) {
        if (s != BOX && this->info_.state_accepting_[s]) {
          new_breakpoint.insert(s);
        }
      }

      ranking g_prime = g;
      for (auto &pr : g_prime) {
        if (pr.first != BOX && !this->info_.state_accepting_[pr.first] &&
            kofola::is_in(pr.first, M)) {
          pr.second = pr.second - 1;
        }
      }

      mstate_rank new_state(
        {},               // reachable states (S)
        false,            // is it Waiting?
        new_breakpoint,   // breakpoint (O)
        g_prime,          // ranking (f)
        src_rank->i_,     // index of tracked rank (i)
        true);            // active
      eta_4.push_back(new_state);
    }

    std::vector<mstate_rank> U;
    if (eta_3.size() > 0 and
        eta_3[0].breakpoint_.size() == 0 and
        eta_3[0].i_ == eta_3[0].f_.get_max_rank() - 1) {
      U.push_back(eta_3[0]);
    }
    if (eta_4.size() > 0 and
        eta_4[0].breakpoint_.size() == 0 and
        eta_4[0].i_ == eta_4[0].f_.get_max_rank() - 1) {
      U.push_back(eta_4[0]);
    }

    mstate_col_set result;
    if (eta_3.size() > 0 and std::find(U.begin(), U.end(), eta_3[0]) == U.end()) {
      std::shared_ptr<mstate> tmp_mstate(new mstate_rank(eta_3[0]));
      result.push_back({tmp_mstate, {}});
    }

    if (eta_4.size() > 0 and std::find(U.begin(), U.end(), eta_4[0]) == U.end()) {
      std::shared_ptr<mstate> tmp_mstate(new mstate_rank(eta_4[0]));
      result.push_back({tmp_mstate, {}});
    }

    for (const mstate_rank& s : U) { // accepting transitions
      // switch to track
      std::shared_ptr<mstate> new_state(new mstate_rank(
        {},        // reachable states (S)
        false,     // is it Waiting?
        {},        // breakpoint (O)
        s.f_,      // ranking (f)
        -1,        // index of tracked rank (i)
        false));   // active
      result.push_back({new_state, {0}});
    }

    return result;
  }
} // get_succ_active() }}}


complement_rank::~complement_rank()
{ }

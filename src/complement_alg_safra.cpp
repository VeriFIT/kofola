// implementation of determinization-based complementation algorithm for
// nondeterministic accepting SCCs based on Safra's trees and producing parity
// condition in the similar way as in
// Yong Li, Andrea Turrini, Weizhi Feng,
// Moshe Y. Vardi, Lijun Zhang: Divide-and-Conquer Determinization of Büchi
// Automata Based on SCC Decomposition. CAV (2) 2022: 152-173

#include "complement_alg_safra.hpp"
#include "safra_tree.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;
using kofola::safra::safra_tree;

namespace { // {{{
/// partial macrostate for the given component
class mstate_safra : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

  /// the corresponding Safra tree
  safra_tree st_;

public: // METHODS

  /// constructor
  explicit mstate_safra(const safra_tree& st) : st_(st)
  { }

  virtual std::string to_string() const override;
  virtual bool is_active() const override { return true; }
  virtual bool eq(const mstate& rhs) const override;
  virtual bool lt(const mstate& rhs) const override;
  virtual ~mstate_safra() override { };

  friend class kofola::complement_safra;
}; // mstate_safra }}}


// Returns true if lhs has a smaller nesting pattern than rhs
// If lhs and rhs are the same, return false.
// compare backwards
bool nesting_cmp(const std::vector<int> &lhs,
                 const std::vector<int> &rhs)
{
  unsigned m = std::min(lhs.size(), rhs.size());
  auto lit = lhs.rbegin();
  auto rit = rhs.rbegin();
  for (unsigned i = 0; i != m; ++i)
  {
    if (*lit != *rit)
      return *lit < *rit;
    lit++;
    rit++;
  }
  return lhs.size() > rhs.size();
}


// Backward search for obtaining the nesting pattern
// The obtained nesting pattern is in reverse order
bool compare_braces(const std::vector<int>& braces, int a, int b)
{ // {{{
  std::vector<int> a_pattern;
  std::vector<int> b_pattern;
  a_pattern.reserve(a + 1);
  b_pattern.reserve(b + 1);
  unsigned size_a = 0;
  unsigned size_b = 0;
  while (a != b)
  {
    if (a > b)
    {
      a_pattern.emplace_back(a);
      // go to the parent
      a = braces[a];
      size_a ++;
    }
    else
    {
      b_pattern.emplace_back(b);
      // go to the parent
      b = braces[b];
      size_b ++;
    }
  }
  return nesting_cmp(a_pattern, b_pattern);
} // compare_braces }}}


bool is_simulation_bigger(unsigned j, unsigned i, const Simulation& sim)
{ // {{{
  for (const auto& st_st_pair : sim) {
    if (st_st_pair.first == i && st_st_pair.second == j) { return true; }
  }
  return false;
} // is_simulation_bigger() }}}


void simulation_reduce(safra_tree& next, const cmpl_info& info)
{ // {{{
  auto it1 = next.labels_.begin();
  while (it1 != next.labels_.end()) {
    auto old_it1 = it1++;
    for (auto it2 = next.labels_.begin(); it2 != next.labels_.end(); ++it2) {
      if (old_it1 == it2)
        continue;
      unsigned i = old_it1->first;
      unsigned j = it2->first;
      // j simulates i?
      if (!is_simulation_bigger(j, i, info.dir_sim_)) {
        continue;
      }
      int brace_i = old_it1->second;
      int brace_j = it2->second;
      bool remove = false;
      // need to compare there nesting pattern
      // TODO: should we check partition instead?
      unsigned scc_i = info.scc_info_.scc_of(i);
      // std::cout << "State " << i << " brace: " << brace_i <<
      // std::endl; std::cout << "State " << j << " brace: " <<
      // brace_j << std::endl;
      // print_pattern_vec(braces, braces.size());
      if (compare_braces(next.braces_, brace_j, brace_i)) {
        remove = true;
      }
      if (remove) {
        it1 = next.labels_.erase(old_it1);
        break;
      }
    }
  }
} // simulation_reduce() }}}


unsigned determine_color (safra_tree& next)
{ // {{{
  // std::cout << "Compute the color for next: " << next.to_string() << std::endl;
  int min_dcc = INT_MAX;
  int min_acc = INT_MAX;

  int topbrace = next.braces_.size();
  constexpr char is_empty = 1;
  constexpr char is_green = 2;
  std::vector<char> empty_green;
  // initially both empty and green for a brace
  empty_green.assign(next.braces_.size(), is_empty | is_green);

  for (const auto &n : next.labels_)
    if (n.second >= 0) // not top level, top level will not have red or good events?
    {
      int brace = n.second;
      // Step A4: The brace has a state on its own, so it is not a green event
      empty_green[brace] &= ~is_green;
      while (brace >= 0 && (empty_green[brace] & is_empty))
      {
        // Once there is brace associated with a state, it is certainly not empty
        empty_green[brace] &= ~is_empty;
        // backward search until top brace, to its parent
        brace = next.braces_[brace];
      }
    }

  // Step A4 Remove brackets within green pairs
  // for each bracket, find its highest green ancestor
  // 0 cannot be in a green pair, its highest green ancestor is itself
  // Also find red and green signals to emit
  // And compute the number of braces to remove for renumbering
  std::vector<int> highest_green_ancestor;
  highest_green_ancestor.assign(next.braces_.size(), 0);

  std::vector<unsigned> decr_by;
  decr_by.assign(next.braces_.size(), 0);
  unsigned decr = 0;

  for (int b = 0; b < next.braces_.size(); ++b)
  {
    // At first, set it to iself
    highest_green_ancestor[b] = b;
    const int &ancestor = next.braces_[b];
    // Note that ancestor < b
    if (ancestor >= 0 && (highest_green_ancestor[ancestor] != ancestor || (empty_green[ancestor] & is_green)))
    {
      // if we do not go further up to the tree or ancester is green
      highest_green_ancestor[b] = highest_green_ancestor[ancestor];
      empty_green[b] |= is_empty; // mark brace for removal
    }

    if (empty_green[b] & is_empty)
    {
      // Step A5 renumber braces
      ++decr;

      // Any brace above topbrace was added while constructing
      // this successor, so it should not emit any red.
      if (b < topbrace)
        // Step A3 emit red
        min_dcc = std::min(min_dcc, b);
    }
    else if (empty_green[b] & is_green)
    {
      assert(b < topbrace);
      // Step A4 emit green
      min_acc = std::min(min_acc, b);
    }

    decr_by[b] = decr;
  }
  // std::cout << "min_dcc = " << min_dcc << ", min_acc = " << min_acc << std::endl;
  // nondet_labellings.emplace_back(min_dcc, min_acc);
  // drease the values
  // Update nodes with new braces numbers
  std::vector<int> newbs = std::vector<int>(next.braces_.size() - decr, -1);
  for (auto& n : next.labels_)
    {
      // if the brace is not -1
      if (n.second >= 0)
        {
          // highest ancester that is green
          unsigned i = highest_green_ancestor[n.second];
          int j = next.braces_[i] >=0 // except 0, every ancester should be this way
                    ? next.braces_[i] - decr_by[next.braces_[i]]
                    : -1;
          n.second = i - decr_by[i];
          // succ.ordered_states_[n.first] = n.second;
          newbs[n.second] = j;
        }
    }
  // std::cout << "done color" << std::endl;
  next.braces_ = newbs;
  int parity;
  if (min_dcc == INT_MAX && min_acc != INT_MAX)
  {
    // only good events
    parity = 2 * (min_acc + 1);
  }
  else if (min_dcc != INT_MAX && min_acc == INT_MAX)
  {
    // only bad events
    parity = 2 * min_dcc + 1;
  }
  else if (min_dcc != INT_MAX && min_acc != INT_MAX)
  {
    // have both good and bad events
    parity = std::min(2 * min_dcc + 1, 2 * min_acc + 2);
  }
  else // both are
  {
    parity = (-1);
  }
  // std::cout << "Color: " << parity << std::endl;
  return parity;
} // determine_color () }}}


spot::acc_cond::acc_code
get_acc_formula(unsigned base, bool odd, unsigned num_colors) {
  // should have the minimal and maximal colors
  // return spot::acc_cond::acc_code::inf({base});
  // assert((num_colors & 1) == odd);
  spot::acc_cond::acc_code res = spot::acc_cond::acc_code::f();

  //    acc-name: parity min even 5
  //    Acceptance: 5 Inf(0) | (Fin(1) & (Inf(2) | (Fin(3) & Inf(4))))
  // built from right to left
  int start = num_colors - 1;
  int inc = -1;
  int end = -1;
  for (int i = start; i != end; i += inc)
  {
    if ((i & 1) == odd)
      res |= spot::acc_cond::acc_code::inf({(unsigned)(i + base)});
    else
      res &= spot::acc_cond::acc_code::fin({(unsigned)(i + base)});
  }

  return res;
}

} // anonymous namespace }}}


std::string mstate_safra::to_string() const
{ // {{{
  std::string res = std::string("[SAFRA(A): ");
  res += this->st_.to_string();
  res += "]";
  return res;
} // to_string() }}}


bool mstate_safra::lt(const mstate& rhs) const
{ // {{{
  const mstate_safra* rhs_safra = dynamic_cast<const mstate_safra*>(&rhs);
  assert(rhs_safra);
  return this->st_ < rhs_safra->st_;
} // lt() }}}


bool mstate_safra::eq(const mstate& rhs) const
{ // {{{
  const mstate_safra* rhs_safra = dynamic_cast<const mstate_safra*>(&rhs);
  assert(rhs_safra);
  return this->st_ == rhs_safra->st_;
} // eq() }}}


complement_safra::complement_safra(const cmpl_info& info, unsigned part_index) :
  abstract_complement_alg(info, part_index)
{ }


mstate_set complement_safra::get_init() const
{ // {{{
  std::set<unsigned> init_states;
  unsigned orig_init = this->info_.aut_->get_init_state_number();
  if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
    init_states.insert(orig_init);
  }

  mstate_set result;
  safra_tree stree;
  for (unsigned st : init_states) {
    stree.labels_.push_back({st, 0});
    stree.braces_.push_back(-1);
  }

  std::shared_ptr<mstate> ms(new mstate_safra(stree));
  result.push_back(ms);
  return result;
} // get_init() }}}


mstate_col_set complement_safra::get_succ_track(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  assert(false);
} // get_succ_track() }}}


mstate_set complement_safra::lift_track_to_active(const mstate* src) const
{ // {{{
  const mstate_safra* src_safra = dynamic_cast<const mstate_safra*>(src);
  assert(src_safra);

  std::shared_ptr<mstate> ms(new mstate_safra(src_safra->st_));
  return {ms};
} // lift_track_to_active() }}}


mstate_col_set complement_safra::get_succ_active(
  const std::set<unsigned>&  glob_reached,
  const mstate*              src,
  const bdd&                 symbol) const
{ // {{{
  const mstate_safra* src_safra = dynamic_cast<const mstate_safra*>(src);
  assert(src_safra);

  // all states in the scc_index should be on the same order
  const int min_new_brace = src_safra->st_.braces_.size();
  DEBUG_PRINT_LN("src: " + src_safra->to_string());
  DEBUG_PRINT_LN("glob_reached: " + std::to_string(glob_reached));
  DEBUG_PRINT_LN("part_index_: " + std::to_string(this->part_index_));

  // std::map<unsigned, int> curr_nodes;
  std::vector<int> braces = src_safra->st_.braces_;
  std::map<unsigned, int> succ_nodes;
  std::set<unsigned> succs;

  // first deal with all states already in the SCCs
  for (const auto &node : src_safra->st_.labels_) {
    const unsigned state = node.first;
    for (const auto &tr : this->info_.aut_->out(state)) {
      const unsigned dst = tr.dst;
      if (!kofola::is_in(dst, glob_reached)) { continue; }

      DEBUG_PRINT_LN("curr s: " + std::to_string(state) +
        " label: " + std::to_string(node.second));
      DEBUG_PRINT_LN("succ t: " + std::to_string(dst));

      // OL: changed from scc_index (maybe should be checked)
      const unsigned succ_part = this->info_.st_to_part_map_.at(dst);

      // the brace of this state
      int newb = node.second;
      // Only care about the states in the current SCCs
      if (this->part_index_ == succ_part) {
        succs.insert(dst);

        DEBUG_PRINT_LN("tr.acc = " + std::to_string(tr.acc));
        DEBUG_PRINT_LN("this->info_.st_to_part_map_.at(state) = " +
          std::to_string(this->info_.st_to_part_map_.at(state)));
        DEBUG_PRINT_LN("succ_part = " + std::to_string(succ_part));
        if (tr.acc || this->info_.st_to_part_map_.at(state) != succ_part)
        { // accepting edges or leaving the partition (TODO: should be SCC?)
          // Step A1: Accepting edges generate new braces
          newb = braces.size();
          // put current brace node.second as the parent of newb
          braces.emplace_back(node.second);
        }
        auto i = succ_nodes.emplace(dst, newb);
        DEBUG_PRINT_LN("succ_nodes: " + std::to_string(succ_nodes));
        if (!i.second) // dst already exists
        {
          // Step A2: Only keep the smallest nesting pattern.
          if (compare_braces(braces, newb, i.first->second)) {
            DEBUG_PRINT_LN("compare_braces(" + std::to_string(braces) + ", " +
                std::to_string(newb) + ", " + std::to_string(i.first->second) +
                ") returned true");
            // newb is smaller
            i.first->second = newb;
            // std::cout << "label of " << dst << " updated to " << newb << std::endl;
          } else {
            DEBUG_PRINT_LN("compare_braces(" + std::to_string(braces) + ", " +
                std::to_string(newb) + ", " + std::to_string(i.first->second) +
                ") returned false");
            // the newb is not smaller than current one
            // new brace was created but is not needed
            if (newb != node.second)
              braces.pop_back();
          }
        }
      }
    }
  }

  // newly incoming states
  std::set<unsigned> reach_diff = kofola::get_set_difference(glob_reached, succs);
  DEBUG_PRINT_LN("newly incoming states: " + std::to_string(reach_diff));

  // std::cout << "After computation of nondet inside " << i << " size = " <<
  // next_nondetstates.size() << "\n";
  // New incoming states Top level is 0, if we enter the SCC, we need more
  // braces Order each entry states since each run can have accepting runs
  for (unsigned dst : reach_diff) {
    // put them all in top brace 0
    int newb = braces.size();
    // Step A1
    auto i = succ_nodes.emplace(dst, newb);
    DEBUG_PRINT_LN("succ_nodes: " + std::to_string(succ_nodes));
    // If the state has not been added
    // means dst is already existing
    if (i.second || i.first->second == -1) {
      // If the state has not been added
      braces.push_back(-1); // top level, so parent is -1
      i.first->second = newb;
    }
  }

  DEBUG_PRINT_LN("succ_nodes: " + std::to_string(succ_nodes));
  // now store the results to succ
  safra_tree next;
  for (auto &node : succ_nodes) {
    next.labels_.emplace_back(std::make_pair(node.first, node.second));
    // std::cout << "succ state = " << node.first << " brace = " << node.second << std::endl;
  }
  // replace the braces
  next.braces_ = braces;
  // use simulation relation to delete states
  if (this->info_.options_.dir_sim) {
    simulation_reduce(next, this->info_);
  }

  // now compute the colour
  unsigned colour = determine_color(next);
  std::shared_ptr<mstate> ms(new mstate_safra(next));

  DEBUG_PRINT_LN("Done computing color for trans to " + ms->to_string() + ": "
    + std::to_string(colour));
  this->min_colour = std::min(this->min_colour, static_cast<int>(colour));
  this->max_colour = std::max(this->max_colour, static_cast<int>(colour));

  return {{ms, {colour}}};
} // get_succ_active() }}}


spot::acc_cond complement_safra::get_acc_cond() const
{ // {{{
  if (this->max_colour < 0) { // no colour was generated
    assert(false);
  }

  DEBUG_PRINT_LN("max_color = " + std::to_string(this->max_colour) +
    ", min_color = " + std::to_string(this->min_colour));

  int local_max_colour = (((this->max_colour % 2) == 1)? this->max_colour + 1 : this->max_colour);
  bool is_min_odd = ((this->min_colour % 2) == 1);

  // TODO: OL: I'm not sure what to do with this
  // for (auto& pair : trans2colors_) {
  //   auto& d = res_->edge_data(pair.first);
  //   auto & vec = pair.second;
  //   // res_->edge_vector()
  //   for (unsigned ndx = 0; ndx < vec.size(); ndx ++) {
  //     // std::cout << ndx << "-th NAC: " << d.acc << std::endl;
  //     if (vec.at(ndx) < 0)
  //     {
  //       // std::cout << "-1 color is " << vec[ndx] << std::endl;
  //       // maximal color must be even color
  //       d.acc.set(nac_color_range[ndx].second - nac_color_range[ndx].first);
  //     }
  //     else
  //     {
  //       d.acc.set(((unsigned)(vec[ndx] - nac_color_range[ndx].first)));
  //       // std::cout << ">0 color is " << vec[ndx] << std::endl;
  //     }
  //     // std::cout << ndx << "-th NAC: " << d.acc << std::endl;
  //   }
  // }

  return ::get_acc_formula(0, is_min_odd, local_max_colour - this->min_colour + 1);
} // get_acc_cond() }}}


complement_safra::~complement_safra()
{ }

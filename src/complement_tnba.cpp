// Copyright (C) 2017-2019 Laboratoire de Recherche et Développement
// de l'Epita.
// Copyright (C) 2022  The COLA Authors
//
// COLA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// COLA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

//#include "optimizer.hpp"
#include "kofola.hpp"
#include "complement_mstate.hpp"
#include "complement_class.hpp"
#include "simulation.hpp"
#include "types.hpp"
//#include "struct.hpp"
#include "decomposer.hpp"
#include "rankings.hpp"
#include "complement_mh.hpp"
#include "complement_ncsb.hpp"
#include "complement_rank.hpp"
#include "abstract_complement_alg.hpp"

#include <deque>
#include <map>
#include <set>
#include <stack>

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twa/acc.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/cleanacc.hh>
#include <spot/twaalgos/postproc.hh>
#include <spot/misc/bddlt.hh>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/misc/version.hh>
#include <spot/twa/acc.hh>

// Complementation of Buchi automara based on SCC decomposition
// We classify three types of SCCs in the input NBA:
// 1. inherently weak SCCs (IWCs): every cycle in the SCC will not visit accepting transitions or every cycle visits an accepting transition
// 2. deterministic accepting SCCs (DACs): states in the SCC have at most one successor remain in the same SCC for a letter
// 3. nondeterministic accepting SCCs (NACs): has an accepting transition and nondeterministic

namespace cola
{
  // complementation Buchi automata
  class tnba_complement
  {
  private:
    // The source automaton.
    const spot::const_twa_graph_ptr aut_;

    // Direct simulation on source automaton.
    std::vector<std::pair<unsigned, unsigned>> dir_sim_;

    // SCCs information of the source automaton.
    spot::scc_info &si_;

    // Complement decomposition options
    compl_decomp_options decomp_options_;

    // Number of states in the input automaton.
    unsigned nb_states_;

    // state_simulator
    state_simulator simulator_;

    // delayed simulation
    delayed_simulation delayed_simulator_;

    // The parity automata being built.
    spot::twa_graph_ptr res_;

    // the number of indices
    unsigned sets_ = 0;

    unsigned num_colors_;

    spot::option_map &om_;

    // use ambiguous
    bool use_unambiguous_;

    bool use_scc_;

    // use stutter
    bool use_stutter_;

    bool use_simulation_;

    // Association between labelling states and state numbers of the
    // DPA.
    std::unordered_map<complement_mstate, unsigned, complement_mstate_hash> rank2n_;

    // States to process.
    std::deque<std::pair<complement_mstate, unsigned>> todo_;

    // Support for each state of the source automaton.
    std::vector<bdd> support_;

    // Propositions compatible with all transitions of a state.
    std::vector<bdd> compat_;

    // is accepting for states
    std::vector<bool> is_accepting_;

    // Whether a SCC is deterministic or not
    std::string scc_types_;

    // State names for graphviz display
    std::vector<std::string> *names_;

    // the index of each weak SCCs
    std::vector<unsigned> weaksccs_;
    // the index of each deterministic accepting SCCs
    std::vector<unsigned> acc_detsccs_;
    // the index of each deterministic accepting SCCs
    std::vector<unsigned> acc_nondetsccs_;

    // Show Rank states in state name to help debug
    bool show_names_;

    std::map<std::pair<std::set<unsigned>, std::set<unsigned>>, unsigned> rank_bounds_; // TODO

    std::string
    get_det_string(const std::vector<state_rank> &states)
    {
      std::string res = "[";
      bool first_state = true;
      for (unsigned p = 0; p < states.size(); p++)
      {
        if (!first_state)
          res += " < ";
        first_state = false;
        res += std::to_string(states[p].first);
      }
      res += "]";
      return res;
    }

    std::string
    get_name(const complement_mstate &ms)
    {
      std::string name;
      name += "(";
      name += get_set_string(std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end()));
      name += "),(";
      for (auto partial : ms.iw_sccs_)
      {
        const std::set<unsigned> tmp(partial.begin(), partial.end());
        name += get_set_string(tmp);
      }
      name += ",";
      const std::set<unsigned> breakset(ms.iw_break_set_.begin(), ms.iw_break_set_.end());
      name += get_set_string(breakset);
      name += "),(";
      for (auto partial : ms.acc_detsccs_)
      {
        const std::set<unsigned> tmp(partial.first.begin(), partial.first.end());
        name += get_set_string(tmp);
        name += "+";
        const std::set<unsigned> aux(partial.second.begin(), partial.second.end());
        name += get_set_string(aux);
      }
      name += ",";
      const std::set<unsigned> det_breakset(ms.det_break_set_.begin(), ms.det_break_set_.end());
      name += get_set_string(det_breakset);
      name += "),(";
      for (auto partial : ms.na_sccs_)
      {
        name += get_set_string_box(partial.reachable);
        name += "+";
        name += partial.f.get_name();
        name += "+";
        name += get_set_string_box(partial.O);
        name += "+";
        name += std::to_string(partial.i);
        name += ",";
      }
      name += "),(";
      name += std::to_string(ms.active_index_);
      name += ")";
      return name;
    }

    // From a Rank state, looks for a duplicate in the map before
    // creating a new state if needed.
    unsigned
    new_state(complement_mstate &s)
    {
      complement_mstate dup(s);
      auto p = rank2n_.emplace(dup, 0);
      if (p.second) // This is a new state
      {
        p.first->second = res_->new_state();
        if (show_names_)
        {
          names_->push_back(get_name(p.first->first));
        }
        todo_.emplace_back(dup, p.first->second);
      }
      return p.first->second;
    }

    bool exists(complement_mstate &s)
    {
      return rank2n_.end() != rank2n_.find(s);
    }

    std::set<unsigned> get_all_successors(std::set<unsigned> current_states, bdd symbol)
    {
      std::set<unsigned> successors;

      for (unsigned s : current_states)
      {
        for (const auto &t : aut_->out(s))
        {
          if (!bdd_implies(symbol, t.cond))
            continue;

          successors.insert(t.dst);
        }
      }

      return successors;
    }

    spot::twa_graph_ptr
    postprocess(spot::twa_graph_ptr aut)
    {
      spot::scc_info da(aut, spot::scc_info_options::ALL);
      // set of states -> the forest of reachability in the states.
      mstate_equiv_map set2scc;
      // record the representative of every SCC
      for (auto p = rank2n_.begin(); p != rank2n_.end(); p++)
      {
        const state_set set = p->first.get_reach_set();
        // first the set of reached states
        auto val = set2scc.emplace(set, state_set());
        val.first->second.insert(p->second);
      }
      mstate_merger merger(aut, set2scc, da, om_);
      spot::twa_graph_ptr res = merger.run();
      if (om_.get(VERBOSE_LEVEL) >= 1)
        std::cout << "The number of states reduced by mstate_merger: "
                  << (aut->num_states() - res->num_states()) << " {out of "
                  << aut->num_states() << "}" << std::endl;
      return res;
    }

  public:
    tnba_complement(const spot::const_twa_graph_ptr &aut, spot::scc_info &si, spot::option_map &om, std::vector<bdd> &implications, compl_decomp_options &decomp_options)
        : aut_(aut),
          om_(om),
          decomp_options_(decomp_options),
          use_simulation_(om.get(USE_SIMULATION) > 0),
          use_scc_(om.get(USE_SCC_INFO) > 0),
          use_stutter_(om.get(USE_STUTTER) > 0),
          use_unambiguous_(om.get(USE_UNAMBIGUITY) > 0),
          si_(si),
          nb_states_(aut->num_states()),
          support_(nb_states_),
          compat_(nb_states_),
          is_accepting_(aut->num_states(), false),
          simulator_(aut, si, implications, om.get(USE_SIMULATION) > 0),
          delayed_simulator_(aut, om),
          show_names_(om.get(VERBOSE_LEVEL) >= 1)
    {

      if (om.get(VERBOSE_LEVEL) >= 2)
      {
        simulator_.output_simulation();
      }
      res_ = spot::make_twa_graph(aut->get_dict());
      res_->copy_ap_of(aut);
      res_->prop_copy(aut,
                      {
                          false,        // state based
                          false,        // inherently_weak
                          false, false, // deterministic
                          true,         // complete
                          false         // stutter inv
                      });
      // Generate bdd supports and compatible options for each state.
      // Also check if all its transitions are accepting.
      for (unsigned i = 0; i < nb_states_; ++i)
      {
        bdd res_support = bddtrue;
        bdd res_compat = bddfalse;
        bool accepting = true;
        bool has_transitions = false;
        for (const auto &out : aut->out(i))
        {
          has_transitions = true;
          res_support &= bdd_support(out.cond);
          res_compat |= out.cond;
          if (!out.acc)
            accepting = false;
        }
        support_[i] = res_support;
        compat_[i] = res_compat;
        is_accepting_[i] = accepting && has_transitions;
      }
      // obtain the types of each SCC
      scc_types_ = get_scc_types(si_);
      // find out the DACs and NACs
      unsigned nonacc_weak = 0;
      for (unsigned i = 0; i < scc_types_.size(); i++)
      {
        if (is_accepting_weakscc(scc_types_, i))
        {
          weaksccs_.push_back(i);
        }
        else if (is_accepting_detscc(scc_types_, i))
        {
          acc_detsccs_.push_back(i);
        }
        else if (is_accepting_nondetscc(scc_types_, i))
        {
          // accepting nondeterministic scc
          acc_nondetsccs_.emplace_back(i);
        }
      }

      // std::cerr << "IWA: " << weaksccs_.size() << ", DET: " << acc_detsccs_.size() << ", NAC: " << acc_nondetsccs_.size() << std::endl;
    }

    unsigned
    get_num_states()
    {
      return this->nb_states_;
    }

    spot::scc_info &
    get_scc_info()
    {
      return this->si_;
    }

    void compute_simulation()
    {
      // compute simulation
      std::vector<bdd> implications;
      auto aut_tmp = aut_;
      spot::simulation(aut_, &implications, -1);

      // get vector of simulated states
      std::vector<std::vector<char>> implies(
          implications.size(),
          std::vector<char>(implications.size(), 0));
      {
        for (unsigned i = 0; i != implications.size(); ++i)
        {
          if (!si_.reachable_state(i))
            continue;
          unsigned scc_of_i = si_.scc_of(i);
          for (unsigned j = 0; j != implications.size(); ++j)
          {
            // reachable states
            if (!si_.reachable_state(j))
              continue;
            // j simulates i and j cannot reach i
            bool i_implies_j = bdd_implies(implications[i], implications[j]);
            if (i_implies_j)
              dir_sim_.push_back({i, j});
          }
        }
      }
    }

    void get_initial_index(complement_mstate &init_state, int &active_index)
    {
      if (decomp_options_.merge_iwa and is_weakscc(scc_types_, active_index))
      {
        init_state.set_active_index(this->weaksccs_[0]);
      }
      else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, active_index))
      {
        init_state.set_active_index(this->acc_detsccs_[0]);
      }
      else
      {
        if (is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)))
        {
          // initial state in nonaccepting scc
          unsigned tmp_index = active_index;
          do
          {
            active_index = (active_index + 1) % si_.scc_count();

            if (active_index == tmp_index)
              break;
          } while (is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)));
          init_state.set_active_index(active_index);
        }
        else
        {
          init_state.set_active_index(active_index);
        }
      }
    }

    std::set<int> reachable_vertices(std::vector<std::vector<int>> list, std::set<int> from)
    {
      std::set<int> all(from);
      std::stack<int> stack;
      int item;
      for (int it : all)
        stack.push(it);

      while (stack.size() > 0)
      {
        item = stack.top();
        stack.pop();
        for (int dst : list[item])
        {
          if (all.find(dst) == all.end())
          {
            stack.push(dst);
            all.insert(dst);
          }
        }
      }
      return all;
    }

    std::vector<std::set<int>> get_reachable_vector()
    {
      std::vector<std::set<int>> reachable_vector;

      std::vector<std::set<int>> list_set(aut_->num_states());
      std::vector<std::vector<int>> list_vector(aut_->num_states());

      for (unsigned s = 0; s < aut_->num_states(); s++)
      {
        reachable_vector.push_back(std::set<int>());
        list_set[s] = std::set<int>();

        // iterate over all transitions from s
        for (const auto &t : aut_->out(s))
        {
          list_set[s].insert(t.dst);
        }

        list_vector[s] = std::vector<int>(list_set[s].begin(), list_set[s].end());
      }

      for (int s = 0; s < aut_->num_states(); s++)
      {
        std::set<int> tmp({s});
        reachable_vector[s] = reachable_vertices(list_vector, tmp);
      }

      return reachable_vector;
    }

    void get_initial_state(complement_mstate &init_state, int &active_index, unsigned &orig_init, std::vector<std::vector<unsigned>> &iw_sccs, std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> &acc_detsccs, std::vector<rank_state> &na_sccs)
    {
      // weak SCCs
      for (unsigned index : weaksccs_)
      {
        if (index != active_index or active_index != si_.scc_of(orig_init))
          iw_sccs.push_back(std::vector<unsigned>());
        else
          iw_sccs.push_back(std::vector<unsigned>(1, orig_init));
      }
      if (decomp_options_.merge_iwa)
      {
        iw_sccs.clear();
        if (is_weakscc(scc_types_, active_index))
          iw_sccs.push_back(std::vector<unsigned>(1, orig_init));
        else
          iw_sccs.push_back(std::vector<unsigned>());
      }

      // det SCCs
      if (not decomp_options_.merge_det)
      {
        for (unsigned index : acc_detsccs_)
        {
          if (index != active_index or si_.scc_of(orig_init) != active_index)
            acc_detsccs.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
          else
          {
            acc_detsccs.push_back({std::vector<unsigned>(1, orig_init), std::vector<unsigned>()});
          }
        }
      }
      else
      {
        if (is_accepting_detscc(scc_types_, active_index))
        {
          acc_detsccs.push_back({std::vector<unsigned>(1, orig_init), std::vector<unsigned>()});
        }
        else
        {
          acc_detsccs.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
        }
      }

      // nondet accepting SCCs
      for (unsigned index : acc_nondetsccs_)
      {
        if (index != active_index or si_.scc_of(orig_init) != active_index)
        {
          rank_state tmp;
          tmp.reachable.insert(-1);
          na_sccs.push_back(tmp);
        }
        else
        {
          rank_state tmp;
          tmp.reachable.insert(orig_init);
          na_sccs.push_back(tmp);
        }
      }
      init_state.na_sccs_ = na_sccs;

      if ((not decomp_options_.merge_iwa) or (not(is_weakscc(scc_types_, active_index) and not is_accepting_weakscc(scc_types_, active_index))))
        init_state.set_iw_sccs(iw_sccs);
      else
      {
        std::vector<std::vector<unsigned>> tmp;
        tmp.push_back(std::vector<unsigned>());
        init_state.set_iw_sccs(tmp);
      }

      init_state.set_acc_detsccs(acc_detsccs);
      auto acc_detsccs_orig = acc_detsccs;
      init_state.curr_reachable_.push_back(orig_init);

      // get break set for active scc
      if (is_weakscc(scc_types_, active_index) and active_index == si_.scc_of(orig_init))
      {
        if (not decomp_options_.merge_iwa or is_accepting_weakscc(scc_types_, active_index))
          init_state.set_iw_break_set(std::vector<unsigned>(1, orig_init));
        else
          init_state.set_iw_break_set(std::vector<unsigned>());
        init_state.det_break_set_ = std::vector<unsigned>();
      }
      else if (is_accepting_detscc(scc_types_, active_index) and active_index == si_.scc_of(orig_init))
      {
        init_state.det_break_set_ = std::vector<unsigned>(1, orig_init);
        init_state.iw_break_set_ = std::vector<unsigned>();
      }
      else
      {
        init_state.set_iw_break_set(std::vector<unsigned>());
        init_state.det_break_set_ = std::vector<unsigned>();
      }

      if (si_.scc_of(orig_init) != active_index)
        init_state.set_iw_break_set(std::vector<unsigned>());
    }

    spot::twa_graph_ptr
    run()
    {
      if (decomp_options_.iw_sim)
      {
        compute_simulation();
      }

      if (show_names_)
      {
        names_ = new std::vector<std::string>();
        res_->set_named_prop("state-names", names_);
      }

      if (this->weaksccs_.size() == 0)
        decomp_options_.merge_iwa = false;
      if (this->acc_detsccs_.size() == 0)
        decomp_options_.merge_det = false;
      if (this->acc_detsccs_.size() == 0 and this->acc_nondetsccs_.size() == 0)
        decomp_options_.tgba = false; // no TGBA for IW SCCs only

      // complementation algorithm
      // auto acc = res_->set_buchi();
      if (decomp_options_.tgba)
        res_->set_generalized_buchi(2);
      else
        res_->set_generalized_buchi(1);

      // spot::print_hoa(std::cerr, aut_);
      // std::cerr << std::endl << std::endl;

      // initial macrostate
      auto scc_info = get_scc_info();
      complement_mstate init_state(scc_info, acc_detsccs_.size());
      unsigned orig_init = aut_->get_init_state_number();
      int active_index = scc_info.scc_of(orig_init);
      get_initial_index(init_state, active_index);

      std::vector<complement_mstate> all_states;
      std::vector<std::vector<unsigned>> iw_sccs;
      std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs;
      std::vector<rank_state> na_sccs;
      bool acc_edge = false;

      get_initial_state(init_state, active_index, orig_init, iw_sccs, acc_detsccs, na_sccs);

      // std::cerr << "Initial: " << get_name(init_state) << std::endl;
      auto init = new_state(init_state);
      res_->set_init_state(init);

      all_states.push_back(init_state);

      // mh_complement mh(aut_, scc_info, scc_types_, decomp_options_, dir_sim_);
      std::vector<std::set<int>> reachable_vector = get_reachable_vector();

      // rank_complement rank_compl(aut_, scc_info, scc_types_, decomp_options_, dir_sim_, reachable_vector, is_accepting_);

      bool sink_state = false;
      bool is_empty = aut_->is_empty();

      while (!todo_.empty())
      {
        auto top = todo_.front();
        todo_.pop_front();
        complement_mstate ms = top.first;

        // no successors for sink state
        if (ms.active_index_ == -1)
          continue;

        // std::cerr << std::endl
        //           << "State: " << get_name(ms) << std::endl;
        active_index = ms.active_index_;

        // skip nonaccepting sccs
        if (active_index >= 0 and is_weakscc(scc_types_, active_index) and (not is_accepting_weakscc(scc_types_, active_index)) and not is_empty)
        {
          ms.active_index_ = (ms.active_index_ + 1) % si_.scc_count();
          todo_.emplace_back(ms, top.second);
          continue;
        }

        // reachable states
        std::set<unsigned> reachable = std::set<unsigned>(ms.curr_reachable_.begin(), ms.curr_reachable_.end());

        // Compute support of all available states.
        bdd msupport = bddtrue;
        bdd n_s_compat = bddfalse;
        const std::set<unsigned> &reach_set = ms.get_reach_set();
        // compute the occurred variables in the outgoing transitions of ms, stored in msupport
        for (unsigned s : reach_set)
        {
          msupport &= support_[s];
          n_s_compat |= compat_[s];
        }

        bdd all = n_s_compat;
        if (all != bddtrue)
        {
          // direct the rest to sink state
          complement_mstate succ(si_, acc_detsccs_.size());
          succ.active_index_ = -1;
          auto sink = new_state(succ);
          // empty state use 0 as well as the weak ones
          res_->new_edge(top.second, sink, !all);
          if (not sink_state)
          {
            res_->new_edge(sink, sink, !all, {0});
            res_->new_edge(sink, sink, all, {0});
            if (decomp_options_.tgba)
            {
              res_->new_edge(sink, sink, !all, {1});
              res_->new_edge(sink, sink, all, {1});
            }
            sink_state = true;
          }
        }

        while (all != bddfalse)
        {
          bdd letter = bdd_satoneset(all, msupport, bddfalse);
          all -= letter;
          // std::cerr << "Current symbol: " << letter << std::endl;

          std::set<unsigned> all_succ = get_all_successors(reachable, letter);

          bool active_type = true;
          bool active_type2 = true;
          bool no_succ = false;
          bool active_iw = true;

          // na succ
          std::vector<rank_state> na_succ(ms.na_sccs_.size());
          std::vector<std::pair<rank_state, bool>> succ_na;

          std::vector<complement_mstate> new_succ;
          complement_mstate new_succ1(scc_info, acc_detsccs_.size());
          new_succ1.na_sccs_ = na_succ;
          complement_mstate new_succ2(scc_info, acc_detsccs_.size());
          new_succ2.na_sccs_ = na_succ;
          new_succ.push_back(new_succ1);
          new_succ.push_back(new_succ2);

          std::vector<bool> acc_succ;

          // iw succ
          std::vector<std::vector<unsigned>> iw_succ(this->weaksccs_.size());
          new_succ[0].iw_sccs_ = iw_succ;
          new_succ[1].iw_sccs_ = iw_succ;

          if (decomp_options_.merge_iwa /*or iw_succ.size() == 0*/)
          {
            iw_succ.clear();
            iw_succ.push_back(std::vector<unsigned>());
          }

          // det succ
          std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_det_succ;
          if (not decomp_options_.merge_det)
          {
            for (unsigned i = 0; i < this->acc_detsccs_.size(); i++)
            {
              acc_det_succ.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
            }
          }
          else
          {
            acc_det_succ.push_back({std::vector<unsigned>(), std::vector<unsigned>()});
          }

          // scc indices
          std::vector<unsigned> indices;
          if (not decomp_options_.merge_iwa)
            indices.insert(indices.end(), this->weaksccs_.begin(), this->weaksccs_.end());
          else if (this->weaksccs_.size() > 0)
            indices.push_back(this->weaksccs_[0]);
          if (not decomp_options_.merge_det)
            indices.insert(indices.end(), this->acc_detsccs_.begin(), this->acc_detsccs_.end());
          else if (this->acc_detsccs_.size() > 0)
            indices.push_back(this->acc_detsccs_[0]);
          if (this->acc_nondetsccs_.size() > 0)
            indices.insert(indices.end(), this->acc_nondetsccs_.begin(), this->acc_nondetsccs_.end());

          // index of value active_index
          auto it = std::find(indices.begin(), indices.end(), active_index);
          unsigned true_index = std::distance(indices.begin(), it);
          unsigned orig_index = true_index;

          std::vector<complement_mstate> succ_det;

          if (ms.iw_break_set_.size() == 0 and ms.det_break_set_.size() == 0)
            active_type = false;

          bool iwa_done = false;
          bool det_done = false;

          std::vector<std::vector<std::pair<complement_mstate, bool>>> succ;

          for (unsigned i = 0; i < indices.size(); i++)
          {
            true_index = (orig_index + i) % indices.size();

            std::vector<unsigned> index;
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, indices[true_index]))
              index.insert(index.end(), weaksccs_.begin(), weaksccs_.end());
            else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, indices[true_index]))
              index.insert(index.end(), acc_detsccs_.begin(), acc_detsccs_.end());
            else
              index.push_back(indices[true_index]);

            // merge iwa
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, index[0]) and iwa_done)
              continue;
            // merge det
            if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, index[0]) and det_done)
              continue;

            // reachable states in this scc
            std::set<unsigned> reach_track;
            std::set<unsigned> scc_states;
            if (decomp_options_.merge_iwa and is_weakscc(scc_types_, index[0]))
            {
              for (auto i : weaksccs_)
              {
                scc_states.insert(scc_info.states_of(i).begin(), scc_info.states_of(i).end());
              }
            }
            else if (decomp_options_.merge_det and is_accepting_detscc(scc_types_, index[0]))
            {
              for (auto i : acc_detsccs_)
              {
                scc_states.insert(scc_info.states_of(i).begin(), scc_info.states_of(i).end());
              }
            }
            else
            {
              scc_states.insert(scc_info.states_of(index[0]).begin(), scc_info.states_of(index[0]).end());
            }
            std::set_intersection(scc_states.begin(), scc_states.end(), reachable.begin(), reachable.end(), std::inserter(reach_track, reach_track.begin()));

            // successors in this scc
            std::set<unsigned> succ_in_scc;
            std::set_intersection(scc_states.begin(), scc_states.end(), all_succ.begin(), all_succ.end(), std::inserter(succ_in_scc, succ_in_scc.begin()));

            bool active_scc = not(std::find(index.begin(), index.end(), active_index) == index.end() and (not decomp_options_.tgba or not is_weakscc(scc_types_, index[0])));
            bool next_to_active = (true_index == (orig_index + 1) % indices.size());

            if (is_weakscc(scc_types_, index[0]))
            {
              mh_compl mhc(aut_, index, scc_info, ms, decomp_options_, letter, true_index, dir_sim_, reachable_vector, is_accepting_);

              if (active_scc)
                succ.push_back(mhc.get_succ_active());

              else if (next_to_active)
              {
                succ.push_back(mhc.get_succ_track_to_active());
                succ.push_back(mhc.get_succ_track());
              }

              else
                succ.push_back(mhc.get_succ_track());
            }
            else if (is_accepting_detscc(scc_types_, index[0]))
            {
              ncsb_compl ncsb(aut_, index, scc_info, ms, decomp_options_, letter, true_index - ms.iw_sccs_.size(), dir_sim_, reachable_vector, is_accepting_);

              if (active_scc)
              {
                succ.push_back(ncsb.get_succ_active());
              }

              else if (next_to_active)
              {
                succ.push_back(ncsb.get_succ_track_to_active());
                succ.push_back(ncsb.get_succ_track());
              }

              else
              {
                succ.push_back(ncsb.get_succ_track());
              }
            }
            else
            {
              rank_comp rank(aut_, index, scc_info, ms, decomp_options_, letter, true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size(), dir_sim_, reachable_vector, is_accepting_);

              if (active_scc)
              {
                succ.push_back(rank.get_succ_active());
              }

              else if (next_to_active)
              {
                succ.push_back(rank.get_succ_track_to_active());
                succ.push_back(rank.get_succ_track());
              }

              else
              {
                succ.push_back(rank.get_succ_track());
              }
            }
          }

          // combine states
          std::vector<std::pair<complement_mstate, bool>> successors;
          // cartesian product
          unsigned k = 0;
          true_index = orig_index;
          for (auto mstate : succ)
          {
            if (k == 0)
            {
              // active component
              for (auto &state : mstate)
              {
                bool iw = not state.first.iw_sccs_.empty();
                bool det = not state.first.acc_detsccs_.empty();
                state.first.iw_sccs_.resize(ms.iw_sccs_.size());
                state.first.acc_detsccs_.resize(ms.acc_detsccs_.size());
                state.first.na_sccs_.resize(ms.na_sccs_.size());
                if (iw)
                {
                  if (true_index != 0)
                  {
                    state.first.iw_sccs_[true_index] = state.first.iw_sccs_[0];
                    state.first.iw_sccs_[0].clear();
                  }
                }
                else if (det)
                {
                  if (true_index - ms.iw_sccs_.size() != 0)
                  {
                    state.first.acc_detsccs_[true_index - ms.iw_sccs_.size()] = state.first.acc_detsccs_[0];
                    state.first.acc_detsccs_[0].first.clear();
                    state.first.acc_detsccs_[0].second.clear();
                  }
                }
                else
                {
                  if (true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size() != 0)
                  {
                    state.first.na_sccs_[true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size()] = state.first.na_sccs_[0];
                    state.first.na_sccs_[0] = rank_state();
                  }
                }
                successors.push_back(state);
              }
            }
            else if (k == 1)
            {
              // active + 1 component - track to active
              true_index = (true_index + 1) % indices.size();
              std::vector<std::pair<complement_mstate, bool>> new_succ;
              for (auto &succ : successors)
              {
                if (succ.second)
                {
                  // track to active
                  for (auto &state : mstate)
                  {
                    complement_mstate tmp(succ.first);
                    if (not state.first.iw_sccs_.empty())
                    {
                      tmp.iw_sccs_[true_index] = state.first.iw_sccs_[0];
                      tmp.iw_break_set_ = state.first.iw_break_set_;
                    }
                    else if (not state.first.acc_detsccs_.empty())
                    {
                      tmp.acc_detsccs_[true_index - ms.iw_sccs_.size()] = state.first.acc_detsccs_[0];
                      tmp.det_break_set_ = state.first.det_break_set_;
                    }
                    else
                    {
                      tmp.na_sccs_[true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size()] = state.first.na_sccs_[0];
                    }
                    tmp.active_index_ = ms.active_index_;
                    unsigned i = 0;
                    do
                    {
                      tmp.active_index_ = indices[(true_index + i) % indices.size()];
                      i++;
                    } while (is_weakscc(scc_types_, tmp.active_index_) and (not is_accepting_weakscc(scc_types_, tmp.active_index_)));
                    new_succ.push_back({tmp, true});
                  }
                }
                else
                  new_succ.push_back(succ);
              }
              successors = new_succ;
            }
            else if (k == 2)
            {
              // active + 1 component - track
              std::vector<std::pair<complement_mstate, bool>> new_succ;
              for (auto &succ : successors)
              {
                if (not succ.second)
                {
                  // track
                  for (auto &state : mstate)
                  {
                    complement_mstate tmp(succ.first);
                    if (not state.first.iw_sccs_.empty())
                    {
                      tmp.iw_sccs_[true_index] = state.first.iw_sccs_[0];
                    }
                    else if (not state.first.acc_detsccs_.empty())
                    {
                      tmp.acc_detsccs_[true_index - ms.iw_sccs_.size()] = state.first.acc_detsccs_[0];
                    }
                    else
                    {
                      tmp.na_sccs_[true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size()] = state.first.na_sccs_[0];
                    }
                    tmp.active_index_ = ms.active_index_;
                    new_succ.push_back({tmp, false});
                  }
                }
                else
                  new_succ.push_back(succ);
              }
              successors = new_succ;
            }
            else
            {
              true_index = (true_index + 1) % indices.size();
              // other components - track
              std::vector<std::pair<complement_mstate, bool>> new_succ;
              for (auto &succ : successors)
              {
                // track
                for (auto &state : mstate)
                {
                  complement_mstate tmp(succ.first);
                  if (not state.first.iw_sccs_.empty())
                  {
                    tmp.iw_sccs_[true_index] = state.first.iw_sccs_[0];
                  }
                  else if (not state.first.acc_detsccs_.empty())
                  {
                    tmp.acc_detsccs_[true_index - ms.iw_sccs_.size()] = state.first.acc_detsccs_[0];
                  }
                  else
                  {
                    tmp.na_sccs_[true_index - ms.iw_sccs_.size() - ms.acc_detsccs_.size()] = state.first.na_sccs_[0];
                  }
                  new_succ.push_back({tmp, succ.second});
                }
              }
              successors = new_succ;
            }
            k++;
          }

          for (unsigned i = 0; i < successors.size(); i++)
          {
            successors[i].first.curr_reachable_ = std::vector<unsigned>(all_succ.begin(), all_succ.end());

            if (std::find(all_states.begin(), all_states.end(), successors[i].first) == all_states.end())
            {
              all_states.push_back(successors[i].first);
              auto s = new_state(successors[i].first);
            }

            // std::cerr << "New succ: " << get_name(successors[i].first) << std::endl;
            auto p = rank2n_.emplace(successors[i].first, 0);
            if (not successors[i].second)
            {
              res_->new_edge(top.second, p.first->second, letter);
              // std::cerr << "Nonaccepting" << std::endl;
            }
            else
            {
              res_->new_edge(top.second, p.first->second, letter, {0});
              // std::cerr << "Accepting" << std::endl;
            }
          }

          //       if (decomp_options_.iw_sim)
          //       {
          //         // simulation on currently reachable states
          //         std::set<unsigned> aux_reach(all_succ.begin(), all_succ.end());
          //         std::set<unsigned> new_reach;
          //         for (auto state : all_succ)
          //         {
          //           if (not is_weakscc(scc_types_, scc_info.scc_of(state)))
          //           {
          //             aux_reach.erase(state);
          //             new_reach.insert(state);
          //           }
          //         }

          //         for (auto pr : dir_sim_)
          //         {
          //           if (pr.first != pr.second and aux_reach.find(pr.first) != aux_reach.end() and aux_reach.find(pr.second) != aux_reach.end() and reachable_vector[pr.second].find(pr.first) == reachable_vector[pr.second].end())
          //             aux_reach.erase(pr.first);
          //         }

          //         aux_reach.insert(new_reach.begin(), new_reach.end());

          //         new_succ[i].curr_reachable_ = std::vector<unsigned>(aux_reach.begin(), aux_reach.end());
          //       }
          //       else
          //       {
          //         new_succ[i].curr_reachable_ = std::vector<unsigned>(all_succ.begin(), all_succ.end());
          //       }

          //       if (is_weakscc(scc_types_, new_succ[i].active_index_) and not decomp_options_.tgba)
          //       {
          //         new_succ[i].det_break_set_.clear();
          //       }

          //       // det sim
          //       if (is_accepting_detscc(scc_types_, new_succ[i].active_index_) and decomp_options_.merge_det and decomp_options_.det_sim)
          //       {
          //         // remove smaller states from S

          //         // all reachable states
          //         std::set<unsigned> new_S;
          //         for (auto item : new_succ[i].acc_detsccs_)
          //         {
          //           new_S.insert(item.first.begin(), item.first.end());
          //           new_S.insert(item.second.begin(), item.second.end());
          //         }

          //         for (auto pr : dir_sim_)
          //         {
          //           if (pr.first != pr.second and new_S.find(pr.first) != new_S.end() and new_S.find(pr.second) != new_S.end())
          //           {
          //             // reachability check
          //             if (reachable_vector[pr.first].find(pr.second) != reachable_vector[pr.first].end() and reachable_vector[pr.second].find(pr.first) == reachable_vector[pr.second].end())
          //             {
          //               // both states in S -> we can remove the smaller one from S
          //               new_S.erase(pr.first);
          //             }
          //           }
          //         }

          //         // erase state if not in new_S
          //         for (auto &item : new_succ[i].acc_detsccs_)
          //         {
          //           std::vector<unsigned> result;
          //           std::set_intersection(item.first.begin(), item.first.end(), new_S.begin(), new_S.end(), std::back_inserter(result));
          //           item.first = result;

          //           std::vector<unsigned> result2;
          //           std::set_intersection(item.second.begin(), item.second.end(), new_S.begin(), new_S.end(), std::back_inserter(result2));
          //           item.second = result2;
          //         }
          //         // erase state from B if not in new_S
          //         std::vector<unsigned> result;
          //         std::set_intersection(new_succ[i].det_break_set_.begin(), new_succ[i].det_break_set_.end(), new_S.begin(), new_S.end(), std::back_inserter(result));
          //         new_succ[i].det_break_set_ = result;
          //       }

          //       // std::cerr << "New succ: " << get_name(new_succ[i]) << std::endl;
          //       if (std::find(all_states.begin(), all_states.end(), new_succ[i]) == all_states.end())
          //       {
          //         all_states.push_back(new_succ[i]);
          //         auto s = new_state(new_succ[i]);
          //       }

          //       auto p = rank2n_.emplace(new_succ[i], 0);
          //       if (active_type and not acc_edge and active_iw)
          //       {
          //         res_->new_edge(top.second, p.first->second, letter);
          //         // std::cerr << "Nonaccepting" << std::endl;
          //       }
          //       else
          //       {
          //         res_->new_edge(top.second, p.first->second, letter, {0});
          //         // std::cerr << "Accepting" << std::endl;
          //       }

          //       if (decomp_options_.tgba and ms.iw_break_set_.size() == 0)
          //         res_->new_edge(top.second, p.first->second, letter, {1});
          //     }
          //}
        }
      }

      // spot::print_hoa(std::cerr, res_);
      // std::cerr << std::endl;

      if (this->acc_detsccs_.size() == 0 and this->acc_nondetsccs_.size() == 0)
        res_ = postprocess(res_);
      return res_;
    }

    spot::twa_graph_ptr
    run_new()
    {

      // TODO: SCC preprocessing
      auto& scc_info = get_scc_info();
      const auto scc_types = get_scc_types(scc_info);

      using abs_cmpl_alg_p = std::unique_ptr<kofola::abstract_complement_alg>;
      using abs_cmpl_ms_p = std::unique_ptr<kofola::abstract_complement_alg::mstate>;
      using vec_algorithms = std::vector<abs_cmpl_alg_p>;
      using vec_macrostates = std::vector<abs_cmpl_ms_p>;

      vec_algorithms algos;

      for (size_t i = 0; i < scc_info.scc_count(); ++i)
      { // determine which algorithms to run on each of the SCCs
        abs_cmpl_alg_p algo;
        if (is_accepting_weakscc(scc_types, i)) {
          // weaksccs_.push_back(i);
        }
        else if (is_accepting_detscc(scc_types, i))
        {
          // acc_detsccs_.push_back(i);
        }
        else if (is_accepting_nondetscc(scc_types, i))
        {
          // acc_nondetsccs_.emplace_back(i);
        }
      }

      assert(false);
    }
  };

  bool
  all_trans_acc(const spot::twa_graph_ptr &aut, unsigned current_state, unsigned scc, spot::scc_info si)
  {
    auto current_scc = si.scc_of(current_state);

    for (auto &t : aut->out(current_state))
    {
      if (si.scc_of(t.dst) == current_scc)
      {
        if (not t.acc)
          return false;
      }
    }

    return true;
  }

  spot::twa_graph_ptr
  saturation(const spot::twa_graph_ptr &aut, spot::scc_info si)
  {
    bool change;
    for (unsigned i = 0; i < si.scc_count(); i++)
    {
      do
      {
        change = false;
        for (auto state : si.states_of(i))
        {
          if (all_trans_acc(aut, state, i, si))
          {
            for (auto s : si.states_of(i))
            {
              for (auto &t : aut->out(s))
              {
                if (t.dst == state)
                {
                  if (not t.acc)
                  {
                    t.acc = spot::acc_cond::mark_t{0};
                    change = true;
                  }
                }
              }
            }
          }
        }
      } while (change);
    }

    return aut;
  }

  spot::twa_graph_ptr
  complement_tnba(const spot::twa_graph_ptr &aut, spot::option_map &om, compl_decomp_options decomp_options)
  {
    const int trans_pruning = om.get(NUM_TRANS_PRUNING);
    // now we compute the simulator
    spot::twa_graph_ptr aut_reduced;
    std::vector<bdd> implications;
    spot::twa_graph_ptr aut_tmp = nullptr;
    if (om.get(USE_SIMULATION) > 0)
    {
      aut_tmp = spot::scc_filter(aut);
      auto aut2 = spot::simulation(aut_tmp, &implications, trans_pruning);
      aut_tmp = aut2;
    }
    if (aut_tmp)
      aut_reduced = aut_tmp;
    else
      aut_reduced = aut;

    spot::scc_info scc(aut_reduced, spot::scc_info_options::ALL);

    if (decomp_options.scc_compl)
    {
      // saturation
      if (decomp_options.sat)
      {
        aut_reduced = saturation(aut_reduced, scc);
        spot::scc_info scc_sat(aut_reduced, spot::scc_info_options::ALL);
        scc = scc_sat;
      }

      // decompose source automaton
      cola::decomposer decomp(aut_reduced, om);
      auto decomposed = decomp.run(true, decomp_options.merge_iwa, decomp_options.merge_det);

      if (decomposed.size() > 0)
      {
        std::vector<spot::twa_graph_ptr> part_res;

        spot::postprocessor p;

        for (auto aut : decomposed)
        {
          if (decomp_options.scc_compl_high)
            p.set_level(spot::postprocessor::High);
          else
            p.set_level(spot::postprocessor::Low);
          // complement each automaton
          auto aut_preprocessed = p.run(aut);
          spot::scc_info part_scc(aut_preprocessed, spot::scc_info_options::ALL);

          auto comp = cola::tnba_complement(aut_preprocessed, part_scc, om, implications, decomp_options);
          auto dec_aut = comp.run();
          // postprocessing for each automaton
          part_res.push_back(p.run(dec_aut));
        }

        // intersection of all complements
        spot::twa_graph_ptr result;
        bool first = true;
        for (auto aut : part_res)
        {
          if (first)
          {
            result = aut;
            first = false;
            continue;
          }
          result = spot::product(result, aut);
        }

        return result;
      }
    }

    spot::const_twa_graph_ptr aut_to_compl;
    aut_to_compl = aut_reduced;

    auto comp = cola::tnba_complement(aut_to_compl, scc, om, implications, decomp_options);
    return comp.run_new();
  }
}

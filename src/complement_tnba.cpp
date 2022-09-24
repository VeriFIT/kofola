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
#include "complement_alg_mh.hpp"
#include "complement_alg_ncsb.hpp"

#include <deque>
#include <map>
#include <set>
#include <stack>

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/product.hh>
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
    kofola::Simulation dir_sim_;

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

          std::set<unsigned> all_succ = kofola::get_all_successors(this->aut_, reachable, letter);

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

    // ######################################################################
    // NEW INTERFACE
    // ######################################################################

    using abs_cmpl_alg_p = std::unique_ptr<kofola::abstract_complement_alg>;
    using vec_algorithms = std::vector<abs_cmpl_alg_p>;

    using abs_cmpl_ms_p = std::shared_ptr<kofola::abstract_complement_alg::mstate>;
    using vec_macrostates = std::vector<abs_cmpl_ms_p>;

    /// the uberstate - combination of all partial macrostates
    class uberstate
    { // {{{
    private:  // DATA MEMBERS

      /// all reached states
      std::set<unsigned> reached_states_;
      vec_macrostates part_macrostates_;

    public:  // METHODS

      /// constructor
      uberstate(const std::set<unsigned>& reached_states,
                const vec_macrostates& part_macrostates) :
        reached_states_(reached_states),
        part_macrostates_(part_macrostates)
      { }

      /// move constructor
      uberstate(uberstate&& us) = default;

      /// deleted copy constructor and assignment operator
      uberstate(const uberstate& us) = delete;
      uberstate& operator=(const uberstate& us) = delete;

      /// converts to string
      std::string to_string() const
      { // {{{
        std::string result;
        result += "<" + std::to_string(this->reached_states_);
        result += " | ";

        for (size_t i = 0; i < this->part_macrostates_.size(); ++i) {
          result += "c" + std::to_string(i) + ": " +
            this->part_macrostates_[i]->to_string();
          if (this->part_macrostates_.size() != i + 1) {
            result += ", ";
          }
        }

        result += ">";
        return result;
      } // to_string() }}}

      /// output stream operator
      friend std::ostream& operator<<(std::ostream& os, const uberstate& us)
      {
        os << us.to_string();
        return os;
      }

      /// returns the set of all reached states
      const std::set<unsigned>& get_reach_set() const
      { return this->reached_states_; }

      /// returns the partial macrostates
      const vec_macrostates& get_part_macrostates() const
      { return this->part_macrostates_; }

      /// total ordering operator to allow use in std::set and std::map
      bool operator<(const uberstate& rhs) const
      { // {{{
        assert(this->part_macrostates_.size() == rhs.part_macrostates_.size());

        // let's start by comparing reached_states_
        // FIXME: inefficient - in C++20, we would use the flying saucer
        // operator <=>
        if (this->reached_states_ != rhs.reached_states_) {
          return this->reached_states_ < rhs.reached_states_;
        }

        // then, compare the partial macrostates lexicographically
        const size_t length = this->part_macrostates_.size();
        for (size_t i = 0; i < length; ++i) {
          if (*(this->part_macrostates_[i]) != *(rhs.part_macrostates_[i])) {
            return *(this->part_macrostates_[i]) < *(rhs.part_macrostates_[i]);
          }
        }

        // they are equal
        return false;
      } // operator< }}}
    }; // uberstate }}}

    /// target of a transition (including colours)
    using state_col = std::pair<unsigned, std::set<unsigned>>;

    /// return type of get_succ (set of pairs <state, set of colours>)
    using vec_state_col = std::vector<state_col>;

    /// functor for comparison of uberstate pointers
    struct uberstate_ptr_less_ftor
    {
      bool operator()(const uberstate* lhs, const uberstate* rhs) const
      { return *lhs < *rhs; }
    };


    // Here we have a bidirectional map between uberstates and state
    // identifiers (unsigned).  The uberstates are physically stored only at
    // 'num_to_uberstate_map_', the reason being that they contain vectors of
    // unique_ptr (no copy is therefore allowed).

    /// maps uberstates to state numbers
    std::map<const uberstate*, unsigned, uberstate_ptr_less_ftor> uberstate_to_num_map_;
    /// maps state numbers to uberstates
    std::vector<uberstate> num_to_uberstate_map_;
    /// counter of states (to be assigned to uberstates) - 0 is reserved for sink
    unsigned cnt_state_ = 0;
    const unsigned SINK_STATE = 0;

    /// accessor into the uberstate table
    unsigned uberstate_to_num(const uberstate& us) const
    { // {{{
      auto it = uberstate_to_num_map_.find(&us);
      if (uberstate_to_num_map_.end() != it) {
        return it->second;
      } else {
        assert(false);
      }
    } // uberstate_to_num() }}}

    /// translates state number to uberstate
    const uberstate& num_to_uberstate(unsigned num) const
    { // {{{
      assert(0 != num);
      assert(num <= num_to_uberstate_map_.size());
      return num_to_uberstate_map_[num - 1];   // offset because of sink
    } // num_to_uberstate() }}}

    /// inserts an uberstate (by moving) and returns its assigned number (if
    /// not present), or just returns the number of an equal uberstate (if
    /// present)
    unsigned insert_uberstate(uberstate&& us)
    { // {{{
      DEBUG_PRINT_LN("inserting uberstate " + us.to_string());
      for (const auto& elem : this->uberstate_to_num_map_) {
        DEBUG_PRINT_LN("  uberstate map element: " + std::to_string(*(elem.first)) +
          " -> " + std::to_string(elem.second));
      }
      DEBUG_PRINT_LN("num_to_uberstate_map_ = " + std::to_string(num_to_uberstate_map_));
      auto it = uberstate_to_num_map_.find(&us);
      if (uberstate_to_num_map_.end() == it) { // not found
        num_to_uberstate_map_.emplace_back(std::move(us));
        assert(num_to_uberstate_map_.size() == cnt_state_ + 1);  // invariant
        const uberstate& us_new = num_to_uberstate_map_[cnt_state_];
        auto jt_bool_pair = uberstate_to_num_map_.insert({&us_new, cnt_state_ + 1});    // increment by one (fixed sink)
        assert(jt_bool_pair.second);    // insertion happened
        ++cnt_state_;
        DEBUG_PRINT_LN("inserted as " + std::to_string(jt_bool_pair.first->second));
        return jt_bool_pair.first->second;
      } else { // found
        DEBUG_PRINT_LN("found as " + std::to_string(it->second));
        return it->second;
      }
    } // insert_uberstate() }}}


    /// computes the Cartesian product of a vector of sets (no repetitions
    /// assumed in the inputs)
    template<class A>
    std::vector<std::vector<A>> compute_cartesian_prod(
      const std::vector<std::vector<A>> vec_of_sets)
    { // {{{
      const size_t length = vec_of_sets.size();
      std::vector<std::vector<A>> result;

      // this vector will iterate over all possible tuples of indices
      std::vector<size_t> indices(length, 0);

      while (true) {
        std::vector<A> vec;
        for (size_t i = 0; i < length; ++i) {
          assert(indices[i] < vec_of_sets[i].size());
          vec.push_back(vec_of_sets[i][indices[i]]);
        }

        assert(vec.size() == length);
        result.push_back(std::move(vec));

        // generate the next vector of indices, if possible
        bool generated = false;
        for (size_t j = 0; j < length; ++j) {
          ++(indices[j]);
          if (indices[j] < vec_of_sets[j].size()) { // indices is set
            generated = true;
            break;
          }
          else { // we need to move into the next index
            indices[j] = 0;
          }
        }

        if (!generated) { break; }
      }

      return result;
    } // compute_cartesian_prod() }}}

    /// removes duplicit values (warning: can change order!)
    template <class T>
    void remove_duplicit(T& t)
    { // {{{
      std::sort(t.begin(), t.end());
      t.erase(std::unique(t.begin(), t.end()), t.end());
    } // remove_duplicit }}}


    /// gets all successors of an uberstate wrt a vector of algorithms and a
    /// symbol
    vec_state_col get_succ_uberstates(
      const vec_algorithms&  algos,
      const uberstate&       src,
      const bdd&             symbol)
    { // {{{
      using mstate_col = kofola::abstract_complement_alg::mstate_col;
      using mstate_col_set = kofola::abstract_complement_alg::mstate_col_set;
      assert(algos.size() == src.get_part_macrostates().size());

      std::set<unsigned> all_succ = kofola::get_all_successors(
          this->aut_, src.get_reach_set(), symbol);
      const vec_macrostates& prev_part_macro = src.get_part_macrostates();
      // this container collects all sets of pairs of macrostates and colours,
      // later, we will turn it into the Cartesian product
      std::vector<mstate_col_set> succ_part_macro_col;
      for (size_t i = 0; i < algos.size(); ++i) {
        const kofola::abstract_complement_alg::mstate* ms = prev_part_macro[i].get();
        mstate_col_set mcs = algos[i]->get_succ_track(all_succ, ms, symbol);
        succ_part_macro_col.emplace_back(std::move(mcs));
      }

      DEBUG_PRINT_LN("generated partial macrostates + colours: " +
        std::to_string(succ_part_macro_col));

      // compute the Cartesian product of the partial macrostates (+ colours)
      std::vector<std::vector<mstate_col>> cp =
        compute_cartesian_prod(succ_part_macro_col);

      // generate uberstates
      vec_state_col result;
      for (const auto& vec : cp) {
        vec_macrostates vm;
        std::set<unsigned> cols;
        for (const auto& ms_col : vec) {
          vm.push_back(ms_col.first);
          // FIXME: renumber colours
          cols.insert(ms_col.second.begin(), ms_col.second.end());
        }

        DEBUG_PRINT_LN("inserting")
        unsigned us_num = insert_uberstate(uberstate(all_succ, vm));
        DEBUG_PRINT_LN("inserted")
        result.emplace_back(us_num, std::move(cols));
      }

      DEBUG_PRINT_LN("computed successors: " + std::to_string(result));

      remove_duplicit(result);

      return result;
    } // get_succ_uberstates() }}}

    /// gets all initial uberstates wrt a vector of algorithms
    std::vector<unsigned> get_initial_uberstates(const vec_algorithms& alg_vec)
    { // {{{
      const size_t length = alg_vec.size();
      std::set<unsigned> initial_states = {aut_->get_init_state_number()};

      using mstate_set = kofola::abstract_complement_alg::mstate_set;
      std::vector<mstate_set> vec_mstate_sets;
      for (auto& alg : alg_vec) { // get outputs of all procedures
        mstate_set init_mstates = alg->get_init();
        assert(1 == init_mstates.size());   // FIXME: we don't support more init states yet!
        if (init_mstates.empty()) { return {};}   // one empty set terminates
        vec_mstate_sets.emplace_back(std::move(init_mstates));
      }

      DEBUG_PRINT_LN("obtained vector of sets of partial macrostates");

      // compute the cartesian product from the vector of sets of macrostates
      std::vector<vec_macrostates> cp = compute_cartesian_prod(vec_mstate_sets);

      DEBUG_PRINT_LN("generated Cartesian product: " + std::to_string(cp));

      std::vector<unsigned> result;
      for (const auto& vec : cp) {
        unsigned us_num = insert_uberstate(uberstate(initial_states, vec));
        result.push_back(us_num);
      }

      // sort and remove duplicates
      remove_duplicit(result);

      return result;
    } // get_initial_uberstates() }}}


    /// selects the algorithms to run on the SCCs
    vec_algorithms select_algorithms(const kofola::cmpl_info& info) const
    { // {{{
      vec_algorithms result;

      for (size_t i = 0; i < info.scc_info_.scc_count(); ++i)
      { // determine which algorithms to run on each of the SCCs
        abs_cmpl_alg_p alg;
        if (!is_accepting_scc(info.scc_types_, i)) { // nonaccepting SCCs can be skipped
          continue;
        }
        else if (is_accepting_weakscc(info.scc_types_, i)) {
          alg = std::make_unique<kofola::complement_mh>(info, i);
        }
        else if (is_accepting_detscc(info.scc_types_, i)) {
          alg = std::make_unique<kofola::complement_ncsb>(info, i);
        }
        else if (is_accepting_nondetscc(info.scc_types_, i)) {
          assert(false);
        }
        else {
          assert(false);
        }

        result.push_back(std::move(alg));
        assert(result.size() == i + 1);
      }

      return result;
    } // select_algorithms() }}}

    /// new modular complementation procedure
    spot::twa_graph_ptr
    run_new()
    { // {{{

      // FIXME: check that what we receive is a normal TBA

      // TODO: SCC preprocessing

      auto& scc_info = get_scc_info();
      const auto scc_types = get_scc_types(scc_info);

      // collect information for complementation
      kofola::cmpl_info info(this->aut_, scc_info, this->dir_sim_, scc_types);

      DEBUG_PRINT_LN("selecting algorithms");

      // creates a vector of algorithms, for every SCC of aut one
      vec_algorithms alg_vec = select_algorithms(info);

      DEBUG_PRINT_LN("algorithms selected");

      DEBUG_PRINT_LN("creating a sink state");

      // our structure for the automaton (TODO: hash table might be better)
      std::map<unsigned, std::vector<std::pair<bdd, vec_state_col>>> compl_states;

      // create a sink state (its transitions)
      // FIXME: should have _all_ colours
      compl_states.insert({SINK_STATE, {{bddtrue, {{SINK_STATE, {0}}}}}});

      // get initial uberstates
      auto init_vec{this->get_initial_uberstates(alg_vec)};
      std::stack<unsigned> todo;
      for (unsigned state : init_vec) {
        compl_states.insert({state, {}});
        todo.push(state);
      }

      DEBUG_PRINT_LN("after get initial");

      DEBUG_PRINT_LN("initial todo: " + std::to_string(todo));

      while (!todo.empty()) { // the main loop
        // get next uberstate
        unsigned us_num = todo.top();
        todo.pop();
        const uberstate& us = num_to_uberstate(us_num);

        // get the post of 'us'
        auto it = compl_states.find(us_num);
        assert(compl_states.end() != it);
        std::vector<std::pair<bdd, vec_state_col>>& us_post = it->second;

        DEBUG_PRINT_LN("processing " + std::to_string(us_num) + ": " + us.to_string());

        // compute support of all available states
        // TODO: this should be cached for each reach_set
        bdd msupport = bddtrue;
        bdd n_s_compat = bddfalse;
        const std::set<unsigned> &reach_set = us.get_reach_set();

        // compute the occurred variables in the outgoing transitions of ms, stored in msupport
        for (unsigned s : reach_set)
        {
          msupport &= support_[s];
          n_s_compat |= compat_[s];
        }

        // direct non-support symbols to sink
        bdd all = n_s_compat;
        if (all != bddtrue) {
          vec_state_col succs = {{SINK_STATE, {}}};
          us_post.emplace_back(std::make_pair(!all, std::move(succs)));
        }

        // iterate over all symbols
        while (all != bddfalse) {
          bdd letter = bdd_satoneset(all, msupport, bddfalse);
          all -= letter;

          vec_state_col succs = get_succ_uberstates(alg_vec, us, letter);
          us_post.emplace_back(std::make_pair(letter, succs));

          for (const auto& state_cols : succs) {
            const unsigned& succ_state = state_cols.first;

            auto it_bool_pair = compl_states.insert({succ_state, {}});
            if (it_bool_pair.second) { // the successor state is new
              todo.push(succ_state);
            }
          }
        }
      }

      DEBUG_PRINT_LN(std::to_string(compl_states));

      // convert the result into a spot automaton
      // FIXME: we should be directly constructing spot aut
      spot::twa_graph_ptr result = spot::make_twa_graph(this->aut_->get_dict());
      for (const auto& st_trans_pair : compl_states) {
        const unsigned& src = st_trans_pair.first;
        unsigned spot_state = result->new_state();
        assert(spot_state == src);

        for (const auto& bdd_vec_tgt_pair : st_trans_pair.second) {
          const bdd& symbol = bdd_vec_tgt_pair.first;
          for (const auto& tgt_col_pair : bdd_vec_tgt_pair.second) {
            const unsigned& tgt = tgt_col_pair.first;
            const std::set<unsigned>& cols = tgt_col_pair.second;
            spot::acc_cond::mark_t spot_cols(cols.begin(), cols.end());
            result->new_edge(src, tgt, symbol, spot_cols);
          }
        }
      }

      spot::print_hoa(std::cerr, result);
      std::cerr << "\n\n\n\n";
      DEBUG_PRINT_LN("FIXME: handle initial states!");

      return result;
    } // run_new() }}}
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

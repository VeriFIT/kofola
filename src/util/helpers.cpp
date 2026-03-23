// Copyright (c) 2017-2020  The Seminator Authors
// Copyright (c) 2021  The COLA Authors
//
// This file is a part of COLA, a tool for determinization
// of omega automata.
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

#include "helpers.hpp"

#include <vector>
#include <sstream>

#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/isweakscc.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/minimize.hh>
#include <spot/misc/optionmap.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/twa/bddprint.hh>
#include <spot/twaalgos/word.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twa/twagraph.hh>

using namespace kofola;

// verbosity of logging
unsigned kofola::LOG_VERBOSITY = 0;

// program options
options kofola::OPTIONS;

namespace helpers
{
  bool
  is_elevator_automaton(const spot::const_twa_graph_ptr &aut)
  {
    spot::scc_info si(aut);
    unsigned nc = si.scc_count();
    for (unsigned scc = 0; scc < nc; ++scc)
    {
      if (is_deterministic_scc(scc, si) || spot::is_inherently_weak_scc(si, scc))
      {
          continue;
      }
      return false;
    }
    return true;
  }

  bool
  is_elevator_automaton(const spot::scc_info &scc, std::string& scc_str)
  {
    for (unsigned sc = 0; sc < scc.scc_count(); ++sc)
    {
      if ((scc_str[sc]&SCC_INSIDE_DET_TYPE) > 0
      || (scc_str[sc]&SCC_WEAK_TYPE) > 0)
      {
          continue;
      }
      return false;
    }
    return true;
  }

  bool
  is_weak_automaton(const spot::const_twa_graph_ptr &aut)
  {
    spot::scc_info si(aut);
    unsigned nc = si.scc_count();
    for (unsigned scc = 0; scc < nc; ++scc)
    {
      if (spot::is_inherently_weak_scc(si, scc))
      {
          continue;
      }
      return false;
    }
    return true;
  }

  bool
  is_weak_automaton(const spot::scc_info &scc, std::string& scc_str)
  {
    for (unsigned sc = 0; sc < scc.scc_count(); ++sc)
    {
      if (scc_str[sc]&SCC_WEAK_TYPE)
      {
          continue;
      }
      return false;
    }
    return true;
  }

  // NOTE: copied from spot/twaalgos/deterministic.cc in SPOT
  //res[i + scccount*j] = 1 iff SCC i is reachable from SCC j
  std::vector<bool>
  find_scc_paths(const spot::scc_info &scc)
  {
    unsigned long scccount = scc.scc_count();
    std::vector<bool> res(scccount * scccount, 0);
    for (unsigned i = 0; i < scccount; ++i)
      {
        // reach itself
        res[i + scccount * i] = true;
      }
    for (unsigned i = 0; i < scccount; ++i)
    {
      unsigned ibase = i * scccount;
      for (unsigned d : scc.succ(i))
      {
        // we necessarily have d < i because of the way SCCs are
        // numbered, so we can build the transitive closure by
        // just ORing any SCC reachable from d.
        unsigned dbase = d * scccount;
        // j reach d (i can reach d, so res[d + i * scccount] = 1)
        for (unsigned j = 0; j < scccount; ++j)
        {
          // j is reachable from i if j is reachable from d
          res[ibase + j] = res[ibase + j] || res[dbase + j];
        }
      }
    }
    return res;
  }

  std::vector<bool>
  get_accepting_reachable_sccs(const spot::scc_info &si)
  {
    unsigned nscc = si.scc_count();
    assert(nscc);
    std::vector<bool> reachable_from_acc(nscc);
    std::vector<bool> res(nscc);
    do // iterator of SCCs in reverse topological order
      {
        --nscc;
        // larger nscc is closer to initial state?
        if (si.is_accepting_scc(nscc) || reachable_from_acc[nscc])
          {
            for (unsigned succ: si.succ(nscc))
              reachable_from_acc[succ] = true;
            res[nscc] = true;
          }
      }
    while (nscc);
    return res;
  }

  bool
  is_limit_deterministic_automaton(const spot::scc_info &si, std::string& scc_str)
  {
    unsigned nscc = si.scc_count();
    assert(nscc);
    std::vector<bool> reachable_from_acc(nscc);
    do // iterator of SCCs in reverse topological order
      {
        --nscc;
        // larger nscc is closer to initial state?
        if ((scc_str[nscc] & SCC_ACC) > 0 || reachable_from_acc[nscc])
          {
            // need to check all outgoing transitions of states in the SCC
            if ((scc_str[nscc] & SCC_DET_TYPE) == 0)
            {
              return false;
            }
            for (unsigned succ: si.succ(nscc))
              reachable_from_acc[succ] = true;
          }
      }
    while (nscc);
    return true;
  }

  std::string
  get_scc_types(const spot::scc_info &si)
  {
    spot::scc_info si_copy = si;
    unsigned nc = si.scc_count();
    std::string res(nc, 0);
    for (unsigned sc = 0; sc < nc; ++sc)
    {
      char type = 0;
      type |= is_deterministic_scc(sc, si) ? SCC_INSIDE_DET_TYPE : 0; // only care about the states inside SCC
      type |= is_deterministic_scc(sc, si, DeterminismScope::ALL) ? SCC_DET_TYPE : 0; // must also be deterministic for all transitions after accepting
      type |= is_deterministic_scc(sc, si, DeterminismScope::BORDER_NONDET) ? SCC_DET_BORDER_NONDET_TYPE : 0;
      type |=  spot::is_inherently_weak_scc(si_copy, sc) ? SCC_WEAK_TYPE : 0;
      type |= si.is_accepting_scc(sc) ? SCC_ACC : 0;
      // other type is 0
      res[sc] = type;
    }

    // Compute predecessors map from SCC successors
    std::vector<std::set<unsigned>> preds(nc);
    for (unsigned sc = 0; sc < nc; ++sc) {
      for (unsigned succ : si.succ(sc)) {
        preds[succ].insert(sc);
      }
    }

    // Fixpoint computation for almost initial deterministic components
    std::vector<bool> is_almost_initial_det(nc, false);
    bool changed = true;
    while (changed) {
      changed = false;
      for (unsigned sc = 0; sc < nc; ++sc) {
        if ((res[sc] & SCC_DET_BORDER_NONDET_TYPE) && !is_almost_initial_det[sc]) {
          bool all_preds_det = true;
          for (unsigned pred : preds[sc]) {
            if (!(res[pred] & SCC_DET_BORDER_NONDET_TYPE) || !is_almost_initial_det[pred]) {
              all_preds_det = false;
              break;
            }
          }
          if (preds[sc].empty() || all_preds_det) {
            is_almost_initial_det[sc] = true;
            res[sc] |= SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE;
            changed = true;
          }
        }
      }
    }

    return res;
  }

  void
  print_scc_types(const std::string& scc_types, const spot::scc_info &scc)
  {
    std::vector<bool> reach_sccs = get_accepting_reachable_sccs(const_cast<spot::scc_info&>(scc));
    for (unsigned i = 0; i < scc.scc_count(); i ++)
    {
      std::cout << "Scc " << i;
      if (scc_types[i] & SCC_WEAK_TYPE)
      {
        std::cout << " weak";
      }
      if (scc_types[i] & SCC_INSIDE_DET_TYPE)
      {
        std::cout << " inside-det";
      }
      if (scc_types[i] & SCC_DET_TYPE)
      {
        std::cout << " det";
      }
      if (scc_types[i] & SCC_INITIAL_DET_TYPE)
      {
        std::cout << " initial-det";
      }
      if (scc_types[i] & SCC_DET_BORDER_NONDET_TYPE)
      {
        std::cout << " det-border-nondet";
      }
      if (scc_types[i] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE)
      {
        std::cout << " initial-almost-det";
      }
      if (scc_types[i] & SCC_ACC)
      {
        std::cout << " accepting";
      }
      std::cout << " " << reach_sccs[i]<< std::endl;
    }
  }

  bool
  is_deterministic_scc(unsigned scc, const spot::scc_info& si,
                     DeterminismScope scope)
  {
    for (unsigned src: si.states_of(scc))
    {
      bdd available = bddtrue;
      bdd border = bddfalse;
      for (auto& t: si.get_aut()->out(src))
      {
        if (scope == DeterminismScope::INSIDE_ONLY && (si.scc_of(t.dst) != scc))
          continue;

        // deterministic inside; nondeterministic on the border (to other SCCs)
        if (scope == DeterminismScope::BORDER_NONDET && (si.scc_of(t.dst) != scc)) {
          border |= t.cond;
          continue;
        }

        if (!bdd_implies(t.cond, available))
          return false;
        else
          available -= t.cond;
      }
      if (scope == DeterminismScope::BORDER_NONDET && !bdd_implies(border, available)) {
        return false;
      }
    }
    return true;
  }

  bool
  is_accepting_scc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_ACC) > 0;
  }

  bool
  is_accepting_detscc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_WEAK_TYPE) == 0 && (scc_types[scc] & SCC_INSIDE_DET_TYPE) > 0 && (scc_types[scc] & SCC_ACC) > 0;
  }

  bool is_accepting_initial_detscc(const std::string& scc_types, unsigned scc)
  {
    return  (scc_types[scc] & SCC_ACC) > 0 && (scc_types[scc] & SCC_INITIAL_DET_TYPE) > 0;
  }

  bool is_accepting_initial_almost_detscc(const std::string& scc_types, unsigned scc) {
    return (scc_types[scc] & SCC_ACC) > 0 && (scc_types[scc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) > 0;
  }

  bool
  is_accepting_weakscc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_WEAK_TYPE) > 0 && (scc_types[scc] & SCC_ACC) > 0 && (scc_types[scc] & SCC_INITIAL_DET_TYPE) == 0;
  }

  bool
  is_weakscc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_WEAK_TYPE) > 0;
  }

  bool
  is_accepting_nondetscc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_WEAK_TYPE) == 0 && (scc_types[scc] & SCC_INSIDE_DET_TYPE) == 0 && (scc_types[scc] & SCC_ACC) > 0;
  }
}

namespace kofola
{
  bool set_contains_accepting_state(
    const std::set<unsigned>&  input,
    const std::vector<bool>&   vec_acceptance)
  {
    return input.end() != std::find_if(input.begin(), input.end(),
        [=](unsigned x) { return vec_acceptance[x]; });
  }

  std::ostream& operator<<(std::ostream& os, const PartitionType& parttype)
  {
    switch (parttype) {
      case PartitionType::INHERENTLY_WEAK: return os << "Inherently weak";
      case PartitionType::DETERMINISTIC: return os << "Deterministic";
      case PartitionType::STRONGLY_DETERMINISTIC: return os << "Strongly deterministic";
      case PartitionType::NONDETERMINISTIC: return os << "Nondeterministic";
      case PartitionType::INITIAL_ALMOST_DETERMINISTIC: return os << "Initial almost deterministic";
      default: throw std::runtime_error("Undefined partition type");
    }
  }
}

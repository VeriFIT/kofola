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

namespace
{
  // Which notion of elevator automaton to check: the classical one
  // only accepts SCCs that are deterministic or inherently weak; the
  // Emerson-Lei one (ELEA) additionally accepts SCCs that are
  // generalized co-Buchi.
  enum class elevator_kind { classic, emerson_lei };

  // Spot's determine_unknown_acceptance() and is_inherently_weak_scc()
  // swap the automaton's acceptance temporarily and restore it from the
  // formula only, so num_sets() can drop below the colours the edges
  // still carry.  Puts the original condition back when destroyed.
  class acc_restorer
  {
    spot::twa_graph_ptr aut_;
    spot::acc_cond acc_;
  public:
    explicit acc_restorer(const spot::const_twa_graph_ptr &aut)
      : aut_(std::const_pointer_cast<spot::twa_graph>(aut)), acc_(aut->acc())
    {}
    ~acc_restorer() { aut_->set_acceptance(acc_); }
  };

  // Resolve SCCs whose acceptance status is ambiguous under mixed
  // Fin/Inf conditions; otherwise is_inherently_weak_scc() can
  // under-report a fully rejecting SCC as not weak (its
  // is_rejecting_scc() fast path relies on this being resolved first).
  // Spot does not support this on alternating automata, so those are
  // left unresolved.
  void
  resolve_unknown_acceptance(spot::scc_info &si)
  {
    if (si.get_aut()->is_existential())
    {
        si.determine_unknown_acceptance();
    }
  }

  // Shared scaffolding for is_elevator_automaton() and
  // is_emerson_lei_elevator_automaton().
  bool
  is_elevator_automaton_aux(const spot::const_twa_graph_ptr &aut,
                            elevator_kind kind)
  {
    // Universal branching is not handled by is_deterministic_scc(),
    // so alternating automata are never considered elevator automata
    // (same convention as spot::is_deterministic()).
    if (!aut->is_existential())
    {
        return false;
    }

    acc_restorer restore(aut);
    spot::scc_info si(aut);
    resolve_unknown_acceptance(si);
    unsigned nc = si.scc_count();
    for (unsigned scc = 0; scc < nc; ++scc)
    {
      bool ok = is_deterministic_scc(scc, si)
        || spot::is_inherently_weak_scc(si, scc);
      if (!ok && kind == elevator_kind::emerson_lei)
      {
          ok = is_generalized_co_buchi_scc(scc, si);
      }
      if (!ok)
      {
          return false;
      }
    }
    return true;
  }

  // Same scaffolding as above, but working off a precomputed per-SCC
  // type bitmask (see get_scc_types()) instead of computing
  // properties directly from an scc_info.
  bool
  is_elevator_automaton_aux(const spot::scc_info &scc, std::string& scc_str,
                            elevator_kind kind)
  {
    // Same convention as the const_twa_graph_ptr overload above.
    if (!scc.get_aut()->is_existential())
    {
        return false;
    }

    for (unsigned sc = 0; sc < scc.scc_count(); ++sc)
    {
      char type = scc_str[sc];
      bool ok = (type & SCC_INSIDE_DET_TYPE) > 0 || (type & SCC_WEAK_TYPE) > 0;
      if (!ok && kind == elevator_kind::emerson_lei)
      {
          ok = (type & SCC_GEN_CO_BUCHI_TYPE) > 0;
      }
      if (!ok)
      {
          return false;
      }
    }
    return true;
  }
} // anonymous namespace

  bool
  is_elevator_automaton(const spot::const_twa_graph_ptr &aut)
  {
    return is_elevator_automaton_aux(aut, elevator_kind::classic);
  }

  bool
  is_elevator_automaton(const spot::scc_info &scc, std::string& scc_str)
  {
    return is_elevator_automaton_aux(scc, scc_str, elevator_kind::classic);
  }

  bool
  is_generalized_co_buchi_scc(unsigned scc, const spot::scc_info& si)
  {
    spot::acc_cond::mark_t sets = si.acc_sets_of(scc);
    spot::acc_cond acc = si.get_aut()->acc().restrict_to(sets);
    acc = acc.remove(si.common_sets_of(scc), false);
    // Neither restrict_to() nor remove() lowers num_sets(), but
    // is_generalized_co_buchi() expects Fin over all sets: drop
    // the unused ones and renumber the rest.
    acc = acc.strip(acc.all_sets() - acc.get_acceptance().used_sets(),
                    false);
    return acc.is_generalized_co_buchi();
  }

  bool
  is_emerson_lei_elevator_automaton(const spot::const_twa_graph_ptr &aut)
  {
    return is_elevator_automaton_aux(aut, elevator_kind::emerson_lei);
  }

  bool
  is_emerson_lei_elevator_automaton(const spot::scc_info &scc, std::string& scc_str)
  {
    return is_elevator_automaton_aux(scc, scc_str, elevator_kind::emerson_lei);
  }

  bool
  is_weak_automaton(const spot::const_twa_graph_ptr &aut)
  {
    acc_restorer restore(aut);
    spot::scc_info si(aut);
    resolve_unknown_acceptance(si);
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
    acc_restorer restore(si.get_aut());
    spot::scc_info si_copy = si;
    resolve_unknown_acceptance(si_copy);
    unsigned nc = si.scc_count();

    // spot::is_inherently_weak_scc() temporarily replaces the acceptance of the
    // automaton and restores it from the formula alone, which shrinks
    // num_sets() to the largest colour the formula mentions.  Transitions that
    // carry a colour the formula does not mention would then refer to a
    // non-existent set, so put the original condition back afterwards.
    auto aut = std::const_pointer_cast<spot::twa_graph>(si.get_aut());
    const spot::acc_cond orig_acc = aut->acc();
    std::string res(nc, 0);
    for (unsigned sc = 0; sc < nc; ++sc)
    {
      char type = 0;
      type |= is_deterministic_scc(sc, si) ? SCC_INSIDE_DET_TYPE : 0; // only care about the states inside SCC
      type |= is_deterministic_scc(sc, si, DeterminismScope::ALL) ? SCC_DET_TYPE : 0; // must also be deterministic for all transitions after accepting
      type |= is_deterministic_scc(sc, si, DeterminismScope::BORDER_NONDET) ? SCC_DET_BORDER_NONDET_TYPE : 0;
      type |=  spot::is_inherently_weak_scc(si_copy, sc) ? SCC_WEAK_TYPE : 0;
      type |= is_generalized_co_buchi_scc(sc, si) ? SCC_GEN_CO_BUCHI_TYPE : 0;
      type |= si_copy.is_accepting_scc(sc) ? SCC_ACC : 0;
      // other type is 0
      res[sc] = type;
    }
    aut->set_acceptance(orig_acc);

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
    // For a universal-branching edge, t.dst does not hold a plain
    // state number (it encodes an index into a destination-set
    // vector instead), so si.scc_of(t.dst) below would read out of
    // bounds.  Alternating automata are simply never considered
    // deterministic here.
    if (!si.get_aut()->is_existential())
    {
        return false;
    }

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

  bool is_accepting_initial_almost_detscc(const std::string& scc_types, unsigned scc) {
    return (scc_types[scc] & SCC_ACC) > 0 && (scc_types[scc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) > 0;
  }

  bool
  is_accepting_weakscc(const std::string& scc_types, unsigned scc)
  {
    return (scc_types[scc] & SCC_WEAK_TYPE) > 0 && (scc_types[scc] & SCC_ACC) > 0 && (scc_types[scc] & SCC_INITIAL_ALMOST_DETERMINISTIC_TYPE) == 0;
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

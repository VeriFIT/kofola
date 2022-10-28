// Copyright (C) 2022  The Kofola Authors
//
// Kofola is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Kofola is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <stack>
#include <string>

#include <spot/twaalgos/hoa.hh>
#include <spot/misc/optionmap.hh>
#include <spot/twaalgos/sccinfo.hh>

static const char SCC_WEAK_TYPE = 1;
static const char SCC_INSIDE_DET_TYPE = 2;
static const char SCC_DET_TYPE = 4;
static const char SCC_ACC = 8;

namespace kofola
{ // {{{

using string_to_string_dict = std::map<std::string, std::string>;

/// structure for command line options
struct options
{ // {{{
	std::vector<std::string> filenames;   ///< input files
	std::string operation;                ///< operation to perform with the inputs
	std::string output_type;              ///< desired automaton on the output
	string_to_string_dict params;         ///< generic parameters
}; // options }}}

} // kofola }}}




// for states ranking/labelling
static const int RANK_N = -1; // nondeterministic states
static const int RANK_M = -2; // missing states

enum automaton_type
{
  NONDETERMINISTIC = 0,
  INHERENTLY_WEAK  = 1,
  ELEVATOR         = 2,
  LIMIT_DETERMINISTIC = 4
};

// Complement decomposition options
struct compl_decomp_options
{
  bool merge_iwa = false;
  bool merge_det = false;
  bool tgba = false;
  bool tba = false;
  bool raw = false;
  bool iw_sim = false;
  bool det_sim = false;
  bool scc_compl = false;
  bool scc_compl_high = false;
  bool dir_sim = true;
  bool sat = true;
  bool dataflow = false;
  bool rank_for_nacs = false;
  bool low_red_interm = false;
};

/// macro for debug outputs
#define PRINT_VERBOSE_LVL(lvl, title, x) {\
  if (kofola::LOG_VERBOSITY >= lvl) {\
    std::cerr << title << ": " << x << "\n";\
  }\
}

#define PRINT_VERBOSE_LVL_LN(lvl, title, x) {\
  PRINT_VERBOSE_LVL(lvl, title, __FILE__ << ":" << __func__ << ":" << __LINE__ << ": " << x)\
}

// #define DEBUG_PRINT(x) { std::cerr << "debug: " << x << "\n"; }
#define DEBUG_PRINT(x) { PRINT_VERBOSE_LVL(2, "debug", x); }
#define DEBUG_PRINT_LN(x) { PRINT_VERBOSE_LVL_LN(2, "debug", x); }
#define WARN_PRINT(x) { PRINT_VERBOSE_LVL(1, "warning", x); }

namespace std
{ // {{{

// DECLARATIONS
template <class A> std::string to_string(const A& value);
template <class A> std::string to_string(const std::set<A>& st);
template <class A> std::string to_string(const std::vector<A>& vec);
template <class A> std::string to_string(const std::list<A>& vec);
template <class A> std::string to_string(const std::stack<A>& stck);
template <class A> std::string to_string(const std::function<A>& func);
template <class A> std::string to_string(const std::shared_ptr<A>& ptr);
template <class A, class B> std::string to_string(const std::pair<A, B>& p);
template <class A, class B> std::string to_string(const std::map<A, B>& mp);
template <class A, class B> std::string to_string(const std::unordered_map<A, B>& unmap);
template <class A, class B> std::string to_string(const std::unordered_multimap<A, B>& unmmap);
} // }}}

namespace kofola
{ // {{{
  /// log verbosity
  extern unsigned LOG_VERBOSITY;

  /// type for representing simulation (NB: this is not very efficient... should be changed)
  using Simulation = std::vector<std::pair<unsigned, unsigned>>;

  /// type for mapping states to partitions where they are (-1 represents trivial SCCs)
  /// TODO: a simple vector should suffice (-1 means invalid partition)
  using StateToPartitionMap = std::unordered_map<unsigned, int>;
  /// similar as above but for mapping SCCs to partitions
  using SCCToPartitionMap = std::unordered_map<unsigned, int>;
  /// map of partitions to sets of SCCs they are composed of
  using PartitionToSCCMap = std::unordered_map<unsigned, std::set<unsigned>>;

  /// types of partitions
  enum class PartitionType
  {
    INHERENTLY_WEAK,
    DETERMINISTIC,
    STRONGLY_DETERMINISTIC,
    NONDETERMINISTIC
  };

  /// output stream overloaded operator
  std::ostream& operator<<(std::ostream& os, const PartitionType& parttype);

  /// the type for mapping partition numbers to their types
  using PartitionToTypeMap = std::map<size_t, PartitionType>;

  /// the type for representing a set of states reachable from a given state
  using ReachableVector = std::vector<std::set<unsigned>>;

  /// type for mapping SCCs to a set of SCCs (e.g., representing all
  /// predecessors of the given SCC)
  using SCCToSCCSetMap = std::unordered_map<unsigned, std::set<unsigned>>;

  /// declaration of printer
  template<class Tuple, size_t N>
  struct TuplePrinter;

  /// get all successors of a given set of states over a given symbol
  template <class T>
  std::set<unsigned> get_all_successors(
    const spot::const_twa_graph_ptr&  aut,
    const T&                          current_states,
    const bdd&                        symbol)
  { // {{{
    std::set<unsigned> successors;

    for (unsigned s : current_states) {
      for (const auto &t : aut->out(s)) {
        if (bdd_implies(symbol, t.cond)) { successors.insert(t.dst); }
      }
    }

    DEBUG_PRINT_LN("all successors of " + std::to_string(current_states) +
        " over " + std::to_string(symbol) +
        " = " + std::to_string(successors));
    return successors;
  } // get_all_successors() }}}


  /// get all successors of a given set over a symbol that are in the partition 'part_num'
  template <class T>
  std::set<unsigned> get_all_successors_in_part(
    const spot::const_twa_graph_ptr&  aut,
    const StateToPartitionMap&        st_to_part_map,
    unsigned                          part_num,
    const T&                          current_states,
    const bdd&                        symbol)
  { // {{{
    std::set<unsigned> successors = get_all_successors(aut, current_states, symbol);
    std::set<unsigned> result;
    std::copy_if(successors.begin(), successors.end(), std::inserter(result, result.end()),
        [=](unsigned x){ return part_num == st_to_part_map.at(x); });

    return result;
  } // get_all_successors_in_part() }}}

  /// get all successors of a given set over a symbol that are in the same SCC as the source state
  template <class T>
  std::set<unsigned> get_all_successors_in_scc(
    const spot::const_twa_graph_ptr&  aut,
    const spot::scc_info&             scc_info,
    const T&                          current_states,
    const bdd&                        symbol)
  { // {{{
    std::set<unsigned> successors;

    for (unsigned s : current_states) {
      for (const auto &t : aut->out(s)) {
        if (scc_info.scc_of(s) == scc_info.scc_of(t.dst) && bdd_implies(symbol, t.cond)) {
          successors.insert(t.dst); }
      }
    }

    return successors;
  } // get_all_successors_in_scc() }}}

  /// checks whether an element is in a container with find()
  template <class T, class C>
  inline bool is_in(const T& elem, const C& container)
  { return container.find(elem) != container.end(); }

  /// checks whether an element in a map-container has given value
  template <class T, class V, class C>
  inline bool has_value(const T& elem, const V& value, const C& container)
  {
		auto it = container.find(elem);
		return it != container.end() && it->second == value;
	}


  /// computes the difference of two sets
  template <class T>
  std::set<T> get_set_difference(const std::set<T>& lhs, const std::set<T>& rhs)
  { // {{{
    std::set<T> result;
    std::set_difference(lhs.begin(), lhs.end(),
      rhs.begin(), rhs.end(),
      std::inserter(result, result.end()));

    return result;
  } // get_set_difference() }}}

  /// computes the union of two sets
  template <class T>
  std::set<T> get_set_union(const std::set<T>& lhs, const std::set<T>& rhs)
  { // {{{
    std::set<T> result;
    std::set_union(lhs.begin(), lhs.end(),
      rhs.begin(), rhs.end(),
      std::inserter(result, result.end()));

    return result;
  } // get_set_union() }}}

  /// computes the intersection of two sets
  template <class T>
  std::set<T> get_set_intersection(const std::set<T>& lhs, const std::set<T>& rhs)
  { // {{{
    std::set<T> result;
    std::set_intersection(lhs.begin(), lhs.end(),
      rhs.begin(), rhs.end(),
      std::inserter(result, result.end()));

    return result;
  } // get_set_intersection() }}}

  /// checks whether a set contains at least one accepting state
  bool set_contains_accepting_state(
    const std::set<unsigned>&  input,               // input set
    const std::vector<bool>&   vec_acceptance);     // vectoring denoting accepting states

} // namespace kofola }}}


namespace cola
{

  /// \brief Complementation
  ///
  /// The automaton \a aut should be an elevator automaton for now.
  /// Output a generalized Buchi automaton
  spot::twa_graph_ptr
  complement_tnba(const spot::twa_graph_ptr &aut, spot::option_map &om, compl_decomp_options decomp_options);

  // ============================ helper functions ===================================

  /// \brief Testing whether the input is an elevator automata in which every scc is either deterministic
  /// or inherently weak (i.e., the states/transitions are either all accepting or nonaccepting)
  ///
  /// Output a bool value
  bool
  is_elevator_automaton(const spot::scc_info &scc, std::string& scc_str);

  bool
  is_elevator_automaton(const spot::const_twa_graph_ptr &aut);

  bool
  is_weak_automaton(const spot::const_twa_graph_ptr &aut);

  bool
  is_weak_automaton(const spot::scc_info &scc, std::string& scc_str);

  bool
  is_limit_deterministic_automaton(const spot::scc_info &scc, std::string& scc_str);

  /// \brief Output the set of states
  ///
  std::string
  get_set_string(const std::set<unsigned> &set);

  std::string
  get_set_string_box(const std::set<int> &set);

  /// \brief Compute the reachability of the SCCs
  ///
  ///
  /// Output a vector res such that res[i + scccount*j] = 1 iff SCC i is reachable from SCC j
  std::vector<bool>
  find_scc_paths(const spot::scc_info &scc);
  /// Output a vector res such that res[i + (j+1)*j/2] = 1 iff SCC i is reachable from SCC
  /// Must ensure that j >= i
  std::vector<bool>
  find_scc_paths_(const spot::scc_info &scc);

  /// \brief Output an automaton to a file
  void output_file(spot::const_twa_graph_ptr aut, const char *file);

  std::vector<bool>
  get_deterministic_sccs(const spot::scc_info &scc);

  std::vector<bool>
  get_accepting_reachable_sccs(const spot::scc_info &scc);

  std::string
  get_scc_types(const spot::scc_info &scc);
  // /// \brief Output an automaton to a file
  // std::vector<bool>
  // is_reachable_weak_sccs(const spot::scc_info &scc, state_simulator& sim);
  void
  print_scc_types(const std::string& scc_types, const spot::scc_info &scc);

  // Check the equivalence of the constructed dpa and the input nba
  void
  check_equivalence(spot::const_twa_graph_ptr nba, spot::twa_graph_ptr dpa);

  bool
  is_accepting_scc(const std::string& scc_types, unsigned scc);

  bool
  is_accepting_detscc(const std::string& scc_types, unsigned scc);

  bool
  is_accepting_weakscc(const std::string& scc_types, unsigned scc);

  bool
  is_weakscc(const std::string& scc_types, unsigned scc);

  bool
  is_accepting_nondetscc(const std::string& scc_types, unsigned scc);

  bool
  is_deterministic_scc(unsigned scc, const spot::scc_info& si,
                     bool inside_only = true);

}

// some things are missing in std
namespace std
{ // {{{

// DEFINITIONS
/** Character to string */
inline std::string to_string(char ch)
{ // {{{
	std::string str;
	str += ch;
	return str;
} // to_string(char) }}}

/** String to string */
inline std::string to_string(const std::string& str) { return str; }

/** Vector to string */
template <class A>
std::string to_string(const std::vector<A>& vec)
{ // {{{
	std::string result = "[";
	bool first = true;
	for (const auto& elem : vec)
	{
		if (!first) { result += ", "; }
		first = false;
		result += std::to_string(elem);
	}
	result += "]";

	return result;
} // to_string(std::vector) }}}

/** List to string */
template <class A>
std::string to_string(const std::list<A>& vec)
{ // {{{
	std::string result = "[";
	bool first = true;
	for (auto elem : vec)
	{
		if (!first) { result += ", "; }
		first = false;
		result += std::to_string(elem);
	}
	result += "]";

	return result;
} // to_string(std::list) }}}

// TODO: the following functions are similar

/** unordered_map to string */
template <class A, class B>
std::string to_string(const std::unordered_map<A, B>& unmap)
{ // {{{
	std::string result = "{";
	bool first = true;
	for (auto key_val_pair : unmap)
	{
		if (!first) { result += ", "; }
		first = false;
		result +=
			std::to_string(key_val_pair.first) +
			" -> " +
			std::to_string(key_val_pair.second);
	}
	result += "}";

	return result;
} // to_string(std::unordered_map) }}}

/** map to string */
template <class A, class B>
std::string to_string(const std::map<A, B>& mp)
{ // {{{
	std::string result = "{";
	bool first = true;
	for (auto key_val_pair : mp)
	{
		if (!first) { result += ", "; }
		first = false;
		result +=
			std::to_string(key_val_pair.first) +
			" -> " +
			std::to_string(key_val_pair.second);
	}
	result += "}";

	return result;
} // to_string(std::map) }}}

/** unordered_multimap to string */
template <class A, class B>
std::string to_string(const std::unordered_multimap<A, B>& unmap)
{ // {{{
	std::string result = "{";
	bool first = true;
	for (auto key_val_pair : unmap)
	{
		if (!first) { result += ", "; }
		first = false;
		result +=
			std::to_string(key_val_pair.first) +
			" -> " +
			std::to_string(key_val_pair.second);
	}
	result += "}";

	return result;
} // to_string(std::unordered_multimap) }}}


/** set to string */
template <class A>
std::string to_string(const std::set<A>& st)
{ // {{{
	std::string result = "{";
	bool first = true;
	for (auto elem : st)
	{
		if (!first) { result += ", "; }
		first = false;
		result += std::to_string(elem);
	}
	result += "}";

	return result;
} // to_string(std::set) }}}

/** stack to string */
template <class A>
std::string to_string(const std::stack<A>& stck)
{ // {{{
	std::stack<A> copy = stck;
	std::vector<A> vec;
	while (!copy.empty()) {
		vec.push_back(copy.top());
		copy.pop();
	}
	std::reverse(vec.begin(), vec.end());
	return std::to_string(vec);
} // to_string(std::stack) }}}

/** function to string */
template <class A>
std::string to_string(const std::function<A>& fun)
{ // {{{
	return std::to_string(static_cast<const void*>(&fun));
} // to_string(std::function) }}}

/** shared_ptr to string */
template <class A>
std::string to_string(const std::shared_ptr<A>& ptr)
{ // {{{
  return "@(" + std::to_string(*ptr) + ")";
} // to_string(std::shared_ptr) }}}

/** tuple to string */
template <class... Ts>
std::string to_string(const std::tuple<Ts...>& tup)
{ // {{{
	std::string str = "<";
  str += kofola::TuplePrinter<decltype(tup), sizeof...(Ts)>::print(tup);
	str += ">";

	return str;
} // to_string(std::tuple) }}}

template <class A, class B>
std::string to_string(const std::pair<A, B>& p)
{ // {{{
	return std::to_string(std::tuple<A, B>(p.first, p.second));
} // to_string(std::pair) }}}

/** arbitrary type with the << operator */
template <class A>
std::string to_string(const A& value)
{ // {{{
	std::ostringstream os;
  os << value;
  return os.str();
} // to_string(T) }}}

} // namespace std }}}

namespace kofola
{ // {{{
  // Taken from
  //   http://en.cppreference.com/w/cpp/utility/tuple/tuple_cat
  template<class Tuple, size_t N>
  struct TuplePrinter
  {
    static std::string print(const Tuple& t)
    {
      std::string res = TuplePrinter<Tuple, N-1>::print(t);
      return res + ", " + std::to_string(std::get<N-1>(t));
    }
  };

  template<class Tuple>
  struct TuplePrinter<Tuple, 1>
  {
    static std::string print(const Tuple& t)
    {
      return std::to_string(std::get<0>(t));
    }
  };
} // namespace kofola }}}

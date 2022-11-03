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

// spot
#include <spot/twa/twa.hh>
#include <spot/twaalgos/sccinfo.hh>

namespace kofola
{ // {{{

/// saturates an automaton with accepting marks
spot::twa_graph_ptr saturate(
	const spot::const_twa_graph_ptr&  aut,
	const spot::scc_info&             si);


// Copied from https://stackoverflow.com/questions/216823/how-to-trim-an-stdstring

/// trim from left
inline std::string str_trim_left(std::string s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
		[](unsigned char ch) { return !std::isspace(ch); }
	));
	return s;
}

/// trim from right
inline std::string str_trim_right(std::string s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(),
		[](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
	return s;
}

/// trim from both sides
inline std::string str_trim(const std::string& s)
{
	return str_trim_left(str_trim_right(s));
}


/// factorial of a number
inline size_t factorial(size_t num)
{
	size_t acc = 1;
	for (size_t i = 2; i <= num; ++i) { acc *= i; }
	return acc;
}


/// computes all partial permutations of 'num' elements from the elements in
/// 'vec', together with a list of unused elements of 'vec'
template <typename T>
std::vector<std::pair<std::vector<T>, std::vector<T>>>
	partial_permutations_ext(const std::vector<T>& vec, unsigned num)
{ // {{{
	assert(vec.size() >= num);

	DEBUG_PRINT_LN("partial permutation for " + std::to_string(vec) + " of " + std::to_string(num) + " elements");

	// vector of pairs where the first element of the pair contains the so-far
	// constructed partial permutation (of length <= num) and the second element
	// contains a vector of elements that can still be chosen in the the partial
	// permutation
	using TmpStateType = std::vector<std::pair<std::vector<T>, std::vector<T>>>;

	TmpStateType tmp_set{{{}, vec}};
	for (unsigned i = 0; i < num; ++i) {
		TmpStateType next_tmp_set;
		for (const auto pair_of_vec : tmp_set) {
			const std::vector<T>& perm = pair_of_vec.first;
			const std::vector<T>& rest = pair_of_vec.second;
			assert(!rest.empty());
			DEBUG_PRINT_LN("perm = " + std::to_string(perm) + "; rest = " + std::to_string(rest));

			for (size_t k = 0; k < rest.size(); ++k) {
				std::vector<T> perm_cpy = perm;
				std::vector<T> rest_cpy = rest;
				perm_cpy.push_back(rest[k]);
				rest_cpy.erase(rest.begin() + k);
				std::pair<std::vector<T>, std::vector<T>> new_pair{perm_cpy, rest_cpy};
				next_tmp_set.emplace_back(std::move(new_pair));
			}
		}
		tmp_set = std::move(next_tmp_set);
	}

	assert(tmp_set.size() * kofola::factorial(num) == kofola::factorial(vec.size()));
	return tmp_set;
} // partial_permutations_ext() }}}


/// computes all partial permutations of 'num' elements from the elements in 'vec'
template <typename T>
std::vector<std::vector<T>> partial_permutations(const std::vector<T>& vec, unsigned num)
{ // {{{
	auto perms = partial_permutations_ext(vec, num);
	std::vector<std::vector<T>> res;
	for (auto& pair_of_vec : perms) {
		res.emplace_back(std::move(pair_of_vec.first));
	}

	return res;
} // partial_permutations() }}}



} // kofola }}}

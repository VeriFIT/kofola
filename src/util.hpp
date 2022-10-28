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

} // kofola }}}

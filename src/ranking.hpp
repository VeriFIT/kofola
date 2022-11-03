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

// a ranking for rank-based complementation

#pragma once

// kofola
#include "kofola.hpp"

namespace kofola
{ // {{{

/// a ranking for rank-based complementation
class ranking : public std::map<unsigned, unsigned>
{
private: // DATA MEMBERS

	unsigned max_rank = 0;

private:  // METHODS

	void set_max_rank(unsigned max_rank) { this->max_rank = max_rank; }

public:  // METHODS

	ranking() : std::map<unsigned, unsigned>() { }
	unsigned get_max_rank() const { return max_rank; }
	void check_tight(std::vector<ranking> rankings);
	bool is_bigger(ranking other);
	std::string to_string() const;

	/// returns all elements with the given rank (i.e., the pre-image of 'rank')
	std::set<unsigned> with_rank(unsigned rank) const;

	bool operator==(const ranking& rhs) const;

	bool operator<(const ranking& rhs) const;

	friend std::ostream& operator<<(std::ostream& os, const ranking& r)
	{
		os << r.to_string();
		return os;
	}
};

} // kofola }}}

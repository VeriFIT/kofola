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

// kofola
#include "ranking.hpp"

// standard library
#include <algorithm>

using namespace kofola;


std::string ranking::to_string() const
{
	return std::to_string(*static_cast<const std::map<unsigned, unsigned>*>(this));
}


std::set<unsigned> ranking::with_rank(unsigned rank) const
{
	std::set<unsigned> res;
	for (const auto& elem_rank_pair : *this) {
		if (elem_rank_pair.second == rank) {
			res.emplace(elem_rank_pair.first);
		}
	}
	return res;
}


bool ranking::operator==(const ranking& rhs) const
{ // {{{
	if (this->get_max_rank() != rhs.get_max_rank() ||
	    this->size() != rhs.size()) {
		return false;
	}

	return std::equal(this->begin(), this->end(),
		rhs.begin(), rhs.end());
} // operator == }}}


bool ranking::operator<(const ranking& rhs) const
{
	assert(false); // the implementation below is probably wrong -

	if (this->size() != rhs.size()) {
		return this->size() < rhs.size();
	}

	if (this->get_max_rank() != rhs.get_max_rank()) {
		return this->get_max_rank() < rhs.get_max_rank();
	}

	return std::lexicographical_compare(this->begin(), this->end(),
		rhs.begin(), rhs.end());
}


// bool compare_ranks(std::tuple<int, int, bool> first, std::tuple<int, int, bool> second)
// {
// 	return (std::get<1>(first) < std::get<1>(second));
// }

// std::vector<ranking> cart_product(std::vector<ranking> rankings, std::tuple<int, int, bool> state)
// {
//     int state_name = std::get<0>(state);
//     int max_rank = std::get<1>(state);
//     bool accepting = std::get<2>(state);
//
//     std::vector<ranking> result;
//
//     for (int i=0; i<max_rank; (accepting ? i+=2 : i++))
//     {
//         std::vector<ranking> new_rankings(rankings.begin(), rankings.end());
//         for (auto r : new_rankings)
//         {
//             r[state_name] = i;
//             if (r.get_max_rank() < i)
//                 r.set_max_rank(i);
//             result.push_back(r);
//         }
//     }
//
//     return result;
// }

// void check_tight(std::vector<ranking>& rankings)
// {
//     std::vector<ranking> result;
//     for (auto r : rankings)
//     {
//         bool skip = false;
//         if (r.get_max_rank() % 2 == 1)
//         {
//             // check box
//             int max_rank = -1;
//             if (r.find(-1) != r.end())
//             {
//                 if (r.get_max_rank() != r[-1])
//                     skip = true;
//                 max_rank = r[-1];
//             }
//
//             if (not skip)
//             {
//                 std::set<int> ranks;
//                 for (auto pr : r)
//                 {
//                     if (pr.first != -1 and pr.second == max_rank)
//                     {
//                         skip = true;
//                         break;
//                     }
//                     ranks.insert(pr.second);
//                 }
//                 if (not skip)
//                 {
//                     for (int i=1; i<r.get_max_rank(); i+=2)
//                     {
//                         if (ranks.find(i) == ranks.end())
//                         {
//                             skip = true;
//                             break;
//                         }
//                     }
//                     if (not skip)
//                         result.push_back(r);
//                 }
//             }
//         }
//     }
//
//     rankings = result;
// }
//
// std::vector<ranking> get_tight_rankings(std::vector<std::tuple<int, int, bool>> mp)
// {
//     std::vector<ranking> rankings;
//
//     // get max rank
//     if (mp.size() > 0)
//     {
//         auto max = std::max_element(mp.begin(), mp.end(), compare_ranks);
//         int max_rank = std::get<1>(*max);
//         if (max_rank > 2*mp.size())
//             max_rank = 2*mp.size();
//
//         // odd rank
//         if (max_rank == 0)
//             return rankings;
//         if (max_rank % 2 == 0)
//             max_rank--;
//
//         for (auto state : mp)
//         {
//             if (rankings.size() == 0)
//             {
//                 for (int i=0; i<std::get<1>(state); (std::get<2>(state) ? i+=2 : i++))
//                 {
//                     ranking r;
//                     r[std::get<0>(state)] = i;
//                     r.set_max_rank(i);
//                     rankings.push_back(r);
//                 }
//             }
//             else
//             {
//                 rankings = cart_product(rankings, state);
//             }
//         }
//
//         check_tight(rankings);
//     }
//
//     return rankings;
// }

bool ranking::is_bigger(ranking other)
{
    for (auto it=this->begin(); it!=this->end(); it++)
    {
        if (it->second < other[it->first])
            return false;
    }
    return true;
}

#pragma once

#include <set>
#include <map>
#include <vector>
#include <string>

#include "kofola.hpp"

namespace kofola
{
class ranking : public std::map<int, int>
{
private:
  unsigned max_rank = 0;

public:
  ranking() : std::map<int, int>(){ }
  std::string to_string() const;
  unsigned get_max_rank() const { return max_rank; }
  void set_max_rank(unsigned max_rank){this->max_rank = max_rank; }
  void check_tight(std::vector<ranking> rankings);
  bool is_bigger(ranking other);

  bool operator==(ranking &other) const
  {
      if (this->max_rank != other.max_rank)
          return false;
      for (auto it=this->begin(); it!=this->end(); it++)
      {
          if (it->second != other[it->first])
              return false;
      }
      return true;
  }

  bool operator<(ranking &other) const
  {
    if (this->max_rank == other.max_rank)
    {
        for (auto it=this->begin(); it!=this->end(); it++)
        {
            if (it->second != other[it->first])
                return it->second < other[it->first];
        }
        return false;
    }
    else
    {
        return this->max_rank < other.max_rank;
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const ranking& r)
  {
    os << r.to_string();
    return os;
  }
};


// bool compare_ranks(std::tuple<int, int, bool> first, std::tuple<int, int, bool> second);
// std::vector<ranking> get_tight_rankings(std::vector<std::tuple<int, int, bool>> mp);
// std::vector<ranking> cart_product(std::vector<ranking> rankings, std::tuple<int, int, bool> state);
// std::vector<ranking> get_succ_rankings(std::vector<std::tuple<int, int, bool>> restr, std::set<unsigned> reachable, bdd letter);
}

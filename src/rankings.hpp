#pragma once

#include <set>
#include <map>
#include <vector>
#include <string>

#include "kofola.hpp"

namespace cola
{
    class ranking : public std::map<int, int>
    {
        private:
            unsigned max_rank = 0;

        public:
            ranking() : std::map<int, int>(){};
            std::string get_name();
            unsigned get_max_rank(){return max_rank;};
            void set_max_rank(unsigned max_rank){this->max_rank = max_rank;};
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
    };

    struct rank_state
    {
        std::set<int> reachable;
        ranking f;
        std::set<int> O;
        int i=-1;
        bool track = true;

        bool operator==(const rank_state &other) const
        {
            if (this->reachable == other.reachable)
            {
                if (this->f == other.f)
                {
                    if (this->O == other.O)
                    {
                        if (this->i == other.i)
                        {
                            if (this->track == other.track)
                            {
                                return true;
                            }
                            else
                            {
                                return false;
                            }
                        }
                        else
                        {
                            return false;
                        }
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        bool operator<(const rank_state &other) const
        {
            if (this->reachable == other.reachable)
            {
                if (this->f == other.f)
                {
                    if (this->O == other.O)
                    {
                        if (this->i == other.i)
                        {
                            return this->track < other.track;
                        }
                        else
                        {
                            return this->i < other.i;
                        }
                    }
                    else
                    {
                        return this->O < other.O;
                    }
                }
                else
                {
                    return this->f < other.f;
                }
            }
            else
            {
                return this->reachable < other.reachable;
            }
        }
    };

    static bool compare_ranks(std::tuple<int, int, bool> first, std::tuple<int, int, bool> second);
    std::vector<ranking> get_tight_rankings(std::vector<std::tuple<int, int, bool>> mp);
    std::vector<ranking> cart_product(std::vector<ranking> rankings, std::tuple<int, int, bool> state);
    std::vector<ranking> get_succ_rankings(std::vector<std::tuple<int, int, bool>> restr, std::set<unsigned> reachable, bdd letter);
}

#pragma once

#include "kofola.hpp"
#include "complement_class.hpp"
#include "rankings.hpp"

namespace cola
{
    class rank_comp : public complement_class
    {
    private:
        std::vector<std::set<unsigned>> predecessors_;

    public:
        rank_comp(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, complement_mstate &mstate, compl_decomp_options decomp_options, bdd symbol, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector, std::vector<bool> is_accepting) : complement_class(aut, scc_index, scc_info, mstate, decomp_options, symbol, true_index, dir_sim, reachable_vector, is_accepting)
        {
            std::vector<std::set<unsigned>> predecessors(reachable_vector_.size());
            for (unsigned i = 0; i < reachable_vector.size(); i++)
            {
                unsigned scc_index = scc_info.scc_of(i);
                for (auto state : reachable_vector[i])
                {
                    unsigned succ_index = scc_info.scc_of(state);
                    if (scc_index != succ_index)
                        predecessors[succ_index].insert(scc_index);
                }
            }

            predecessors_ = predecessors;
        }

        complement_mstate getInit() { return complement_mstate(scc_info_, 0); };

        std::vector<std::pair<complement_mstate, bool>> get_succ_active();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active();

        std::set<int> get_successors_with_box(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index);

        void get_max(std::vector<ranking> &rankings);
        std::vector<ranking> get_maxrank(std::set<unsigned> reachable, rank_state state, bdd letter, unsigned scc_index);
        std::vector<ranking> get_succ_rankings(ranking r, std::vector<std::tuple<int, int, bool>> restr, std::set<unsigned> reachable, bdd letter, unsigned scc_index);
    };
}

#pragma once

#include "kofola.hpp"
#include "complement_class.hpp"

namespace cola
{
    // (C, S, B) for complementing DACs
    const int NCSB_C = 2;
    const int NCSB_S = 4;
    const int NCSB_B = 3;

    class ncsb_compl : public complement_class
    {
    public:
        ncsb_compl(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, complement_mstate &mstate, compl_decomp_options decomp_options, bdd symbol, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector, std::vector<bool> is_accepting) : complement_class(aut, scc_index, scc_info, mstate, decomp_options, symbol, true_index, dir_sim, reachable_vector, is_accepting) {}

        complement_mstate getInit() { return complement_mstate(scc_info_, 0); };

        std::vector<std::pair<complement_mstate, bool>> get_succ_active();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active();

        // void csb_successors(const std::vector<state_rank> &curr_det_states, int scc_index, std::vector<int> &next_scc_indices, std::vector<std::map<unsigned, int>> &succ_maps, std::vector<bool> &acc_succs, std::set<unsigned> &next_detstates, std::unordered_map<unsigned, std::vector<std::pair<bool, unsigned>>> &det_cache, unsigned active_index);
    };
}

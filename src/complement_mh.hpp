#pragma once

#include "kofola.hpp"
#include "complement_class.hpp"

namespace cola
{
    class mh_compl : public complement_class
    {
        public:
        mh_compl(const spot::const_twa_graph_ptr aut, std::vector<unsigned> scc_index, spot::scc_info &scc_info, complement_mstate &mstate, compl_decomp_options decomp_options, bdd symbol, unsigned true_index, std::vector<std::pair<unsigned, unsigned>> dir_sim, std::vector<std::set<int>> reachable_vector, std::vector<bool> is_accepting) : complement_class(aut, scc_index, scc_info, mstate, decomp_options, symbol, true_index, dir_sim, reachable_vector, is_accepting) {}

        complement_mstate getInit(){return complement_mstate(scc_info_, 0);};

        std::vector<std::pair<complement_mstate, bool>> get_succ_active();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track();
        std::vector<std::pair<complement_mstate, bool>> get_succ_track_to_active();
    };
}

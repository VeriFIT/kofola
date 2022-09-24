#include "complement_class.hpp"

namespace cola
{
    std::set<unsigned>
    complement_class::get_all_successors(std::vector<unsigned> current_states, bdd symbol)
    {
        std::set<unsigned> successors;

        for (unsigned s : current_states)
        {
            for (const auto &t : aut_->out(s))
            {
                if (!bdd_implies(symbol, t.cond))
                    continue;

                successors.insert(t.dst);
            }
        }

        return successors;
    }

    std::set<unsigned>
    complement_class::get_all_successors_in_scc_same(std::vector<unsigned> current_states, bdd symbol)
    {
        std::set<unsigned> successors;

        for (unsigned s : current_states)
        {
            for (const auto &t : aut_->out(s))
            {
                if (!bdd_implies(symbol, t.cond))
                    continue;

                if (scc_info_.scc_of(t.src) == scc_info_.scc_of(t.dst) and std::find(scc_index_.begin(), scc_index_.end(), scc_info_.scc_of(t.dst)) != scc_index_.end())
                    successors.insert(t.dst);
            }
        }

        return successors;
    }

    std::set<unsigned>
    complement_class::get_succ_acc_trans_scc(std::vector<unsigned> current_states, bdd symbol)
    {
        std::set<unsigned> successors;
        spot::acc_cond::mark_t acc = {0};

        for (unsigned s : current_states)
        {
            for (const auto &t : aut_->out(s))
            {
                if (!bdd_implies(symbol, t.cond))
                    continue;

                // only care about acc transitions inside the same scc
                if (t.acc == acc and scc_info_.scc_of(t.src) == scc_info_.scc_of(t.dst) and std::find(scc_index_.begin(), scc_index_.end(), scc_info_.scc_of(t.dst)) != scc_index_.end())
                    successors.insert((int)t.dst);
            }
        }

        return successors;
    }

    std::set<int>
    complement_class::get_all_successors_acc(std::set<unsigned> current_states, bdd symbol, unsigned scc_index)
    {
        std::set<int> successors;
        spot::acc_cond::mark_t acc = {0};

        for (unsigned s : current_states)
        {
            for (const auto &t : aut_->out(s))
            {
                if (!bdd_implies(symbol, t.cond))
                    continue;

                if (t.acc == acc and scc_info_.scc_of(t.dst) == scc_index)
                    successors.insert((int)t.dst);
            }
        }

        return successors;
    }

    std::set<unsigned>
    complement_class::get_all_successors(std::set<unsigned> current_states, bdd symbol)
    {
        std::set<unsigned> successors;

        for (unsigned s : current_states)
        {
            for (const auto &t : aut_->out(s))
            {
                if (!bdd_implies(symbol, t.cond))
                    continue;

                successors.insert(t.dst);
            }
        }

        return successors;
    }
}

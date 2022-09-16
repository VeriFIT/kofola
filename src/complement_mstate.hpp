#pragma once

#include "kofola.hpp"
#include "types.hpp"
#include "rankings.hpp"

namespace cola
{
    enum slice_mark
    {
        None = -1,
        Inf = 0,
        New = 1,
        Fin = 2
    };

    class complement_mstate
    {
    public:
        // the number of states num, default values, and number of NACs
        complement_mstate(spot::scc_info &si, unsigned num_det_sccs)
            : si_(si)
        {
            for (unsigned i = 0; i < num_det_sccs; i++)
            {
                detscc_ranks_.emplace_back(std::vector<state_rank>());
            }
        }

        complement_mstate(spot::scc_info &si) : si_(si){}

        complement_mstate(const complement_mstate &other)
            : si_(other.si_)
        {
            this->break_set_.clear();
            this->break_set_.insert(other.break_set_.begin(), other.break_set_.end());
            this->weak_set_.clear();
            this->weak_set_.insert(other.weak_set_.begin(), other.weak_set_.end());

            this->detscc_ranks_.clear();
            for (unsigned i = 0; i < other.detscc_ranks_.size(); i++)
            {
                std::vector<state_rank> copy = other.detscc_ranks_[i];
                this->detscc_ranks_.emplace_back(copy);
            }

            this->nondetscc_ranks_.clear();
            this->nondetscc_marks_.clear();
            for (unsigned i = 0; i < other.nondetscc_ranks_.size(); i++)
            {
                std::set<unsigned> copy = other.nondetscc_ranks_[i];
                this->nondetscc_ranks_.emplace_back(copy);
                this->nondetscc_marks_.push_back(other.nondetscc_marks_[i]);
            }

            this->set_iw_sccs(other.iw_sccs_);
            this->set_iw_break_set(other.iw_break_set_);
            this->det_break_set_ = other.det_break_set_;
            this->set_active_index(other.active_index_);
            this->set_acc_detsccs(other.acc_detsccs_);
            this->curr_reachable_ = other.curr_reachable_;

            this->na_sccs_ = other.na_sccs_;
        }

        std::set<unsigned>
        get_reach_set() const;

        std::set<unsigned>
        get_weak_set() const;

        bool is_empty() const;

        int get_max_rank() const;

        bool operator<(const complement_mstate &other) const;
        bool operator==(const complement_mstate &other) const;

        complement_mstate &
        operator=(const complement_mstate &other)
        {
            this->si_ = other.si_;
            this->break_set_.clear();
            this->break_set_.insert(other.break_set_.begin(), other.break_set_.end());
            this->weak_set_.clear();
            this->weak_set_.insert(other.weak_set_.begin(), other.weak_set_.end());

            this->detscc_ranks_.clear();
            for (unsigned i = 0; i < other.detscc_ranks_.size(); i++)
            {
                std::vector<state_rank> copy = other.detscc_ranks_[i];
                this->detscc_ranks_.emplace_back(copy);
            }

            this->nondetscc_ranks_.clear();
            this->nondetscc_marks_.clear();
            for (unsigned i = 0; i < other.nondetscc_ranks_.size(); i++)
            {
                std::set<unsigned> copy = other.nondetscc_ranks_[i];
                this->nondetscc_ranks_.emplace_back(copy);
                this->nondetscc_marks_.push_back(other.nondetscc_marks_[i]);
            }

            this->set_iw_sccs(other.iw_sccs_);
            this->set_iw_break_set(other.iw_break_set_);
            this->det_break_set_ = other.det_break_set_;
            this->set_active_index(other.active_index_);
            this->set_acc_detsccs(other.acc_detsccs_);
            this->curr_reachable_ = other.curr_reachable_;

            this->na_sccs_ = other.na_sccs_;

            return *this;
        }

        size_t hash() const;

        // SCC information
        spot::scc_info &si_;
        // 1. NAC by slice-based complementation
        std::vector<std::set<unsigned>> nondetscc_ranks_;
        std::vector<slice_mark> nondetscc_marks_; // marks for in O or not?

        // 2. DAC by determinization or NCSB
        std::vector<std::vector<state_rank>> detscc_ranks_;

        // 3. IWC states point to RANK_WEAK
        // breakpoint construction for weak accepting SCCs
        std::set<unsigned> weak_set_;
        std::set<unsigned> break_set_;

        std::vector<std::vector<unsigned>> iw_sccs_;
        std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs_;
        std::vector<unsigned> iw_break_set_;
        std::vector<unsigned> det_break_set_;
        int active_index_ = 0;
        std::vector<unsigned> curr_reachable_;

        std::vector<rank_state> na_sccs_;

        void
        set_iw_sccs(std::vector<std::vector<unsigned>> iw_sccs)
        {
            this->iw_sccs_ = iw_sccs;
        }

        void
        set_acc_detsccs(std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>> acc_detsccs)
        {
            this->acc_detsccs_ = acc_detsccs;
        }

        void
        set_iw_break_set(std::vector<unsigned> iw_break_set)
        {
            this->iw_break_set_ = iw_break_set;
        }

        void
        set_active_index(int index)
        {
            this->active_index_ = index;
        }
    };

    struct complement_mstate_hash
    {
        size_t
        operator()(const complement_mstate &s) const noexcept
        {
            return s.hash();
        }
    };


}

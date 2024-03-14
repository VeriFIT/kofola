#pragma once

// kofola
#include "kofola.hpp"
#include "abstract_complement_alg.hpp"
#include "complement_sync.hpp"
#include "abstract_successor.hpp"

// spot
#include <spot/twa/twa.hh>

namespace kofola
{
    class emptiness_check {
    public:
        emptiness_check(abstract_successor *as);

        bool empty();

        /// performs emptiness check: whether aut_A⊆aut_B using tarjan's algo
        void tarjan_is_empty(const std::unique_ptr<abstract_successor::mstate> &src_mstate);

    private:
        abstract_successor *abstr_succ_;

        const int UNDEFINED = -1;

        /// tarjan variables
        std::map<abstract_successor::mstate, signed> dfs_num_;
        std::map<abstract_successor::mstate, bool> on_stack_;
        signed index_ = 0;
        std::stack<abstract_successor::mstate> tarjan_stack_;
        std::stack<abstract_successor::mstate> SCCs_;
        /// end of tarjan variables

        /// to stop searching when counter-example
        bool decided_ = false;
        bool empty_ = true;
    };
}


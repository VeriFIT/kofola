#pragma once

// kofola
#include "kofola.hpp"
#include "abstract_complement_alg.hpp"
#include "complement_sync.hpp"
#include "abstract_successor.hpp"

// spot
#include <spot/twa/twa.hh>

#define INCLUSION 1
#define HYPERLTL_MC_EMPTINESS 2

namespace kofola
{
    class emptiness_check {
    public:
        emptiness_check(abstract_successor *as, int type, bool use_early_subsums);

        bool empty();

        /// performs emptiness check: whether aut_A⊆aut_B using tarjan's algo
        void tarjan_is_empty(const std::shared_ptr<abstract_successor::mstate> &src_mstate, spot::acc_cond::mark_t path_cond);

    private:
        abstract_successor *abstr_succ_;
        int type_;
        unsigned cnt_ = 0;

        const int UNDEFINED = -1;

        struct shared_ptr_comparator {
            template<typename T>
            bool operator()(const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const {
                return *lhs < *rhs;
            }
        };

        /// tarjan variables
        std::map<std::shared_ptr<abstract_successor::mstate>, signed, shared_ptr_comparator> dfs_num_;
        std::map<std::shared_ptr<abstract_successor::mstate>, bool, shared_ptr_comparator> on_stack_;
        signed index_ = 0;
        std::vector<std::shared_ptr<abstract_successor::mstate>> tarjan_stack_;
        std::stack<std::shared_ptr<abstract_successor::mstate>> SCCs_;
        /// end of tarjan variables

        std::vector<std::shared_ptr<abstract_successor::mstate>> empty_lang_states_;

        /// to stop searching when counter-example
        bool decided_ = false;
        bool empty_ = true;

        bool use_early_subsums_ = true;
    };
}


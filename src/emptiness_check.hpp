/**
 * @file emptiness_check.hpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief Declarations for on the fly emptiness check procedure
 * @version 0.1
 * @date 2024-05-03
 * 
 * 
 */

#pragma once

// kofola
#include "kofola.hpp"
#include "abstract_complement_alg.hpp"
#include "complement_sync.hpp"
#include "inclusion_check.hpp"

// spot
#include <spot/twa/twa.hh>

#define INCLUSION 1
#define HYPERLTL_MC_EMPTINESS 2

#define MAX_SUBSUM_BUCKET 100

namespace kofola
{
    class inclusion_check;
    class inclusion_mstate;
    /// @brief  main class that implements on the fly emptiness check
    class emptiness_check {
    public:
        /// constructor that takes abstract_successor to obtain intial states and successors
        emptiness_check(inclusion_check *incl_checker);

        /// returns true if the aut given by abstract_successor is empty (should be called after constructor is called)
        bool empty();

        void update_structures(const std::shared_ptr<inclusion_mstate>& src_mstate);

        /// implements the edited Gaiser and Schwoon algorithm suggested in the thesis
        /// path_cond can be omitted
        /// TODO too deep of a recursion can cause mem. problems, rewrite to iteration
        bool gs_edited(std::shared_ptr<inclusion_mstate> src_mstate);

        /// implements Gaiser and Schwoon algorithm, with the possibility of subsumptions usage
        /// path_cond can be omitted
        /// TODO too deep of a recursion can cause mem. problems, rewrite to iteration
        bool gs(std::shared_ptr<inclusion_mstate> src_mstate);

        /// decides whether there is a state p on the searchpath such that src_mstate is simul. (early or +1) less than p,
        /// without seeing acc. trans. if yes => true
        bool simulation_prunning(const std::shared_ptr<inclusion_mstate> & src_mstate);

        /// returns true if there is simul. bigger state in empty_lang_states_ than dst_mstate, if yes => true
        bool empty_lang(const std::shared_ptr<inclusion_mstate> &dst_mstate);

        /// part of the Gaiser and Schwoon algorithm (adjusted for TGBA) that merges acc. marks within SCC while reaching dst_mstate
        bool merge_acc_marks(const std::shared_ptr<inclusion_mstate> &dst_mstate);

        /// remove states of the SCC with root src_mstate from the stack
        void remove_SCC(const std::shared_ptr<inclusion_mstate> & src_mstate);

        /// used by gs_edited to guess no further exploration (based on thesis)
        bool check_simul_less(const std::shared_ptr<inclusion_mstate> &dst_mstate);

    private:
        inclusion_check *incl_checker_; /// instance of abstract_successor that provides states and transitions
        // unsigned cnt_ = 0; /// number of states, for benchmarks

        const int UNDEFINED = -1;

        struct shared_ptr_comparator {
            template<typename T>
            bool operator()(const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const {
                return *lhs < *rhs;
            }
        };

        /// GS algorithm variables
        std::map<std::shared_ptr<inclusion_mstate>, signed, shared_ptr_comparator> dfs_num_;
        std::map<std::shared_ptr<inclusion_mstate>, bool, shared_ptr_comparator> on_stack_;
        signed index_ = 0;
        std::vector<std::shared_ptr<inclusion_mstate>> tarjan_stack_;
        std::stack<std::shared_ptr<inclusion_mstate>> SCCs_;
        /// end of GS algorithm variables

        std::vector<std::pair<std::shared_ptr<inclusion_mstate>, spot::acc_cond::mark_t>> dfs_acc_stack_; /// this stack could probably be omitted and use SCCs_ instead

        /// stores states from non-trivial SCCs with empty language. Indexed by the states of the first automaton
        std::map<unsigned, std::vector<std::shared_ptr<inclusion_mstate>>> empty_lang_states_;

        /// to stop searching when counter-example
        bool decided_ = false;
        bool empty_ = true;

        #ifdef ENABLE_COUNTER
            unsigned cnt_ = 0;
        #endif

        std::map<std::shared_ptr<inclusion_mstate>,std::set<std::shared_ptr<inclusion_mstate>>> state_jumps_to_cutoffs_; /// for the new approach to jump to cut off states
    };
}


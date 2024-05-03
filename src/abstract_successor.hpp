/**
 * @file abstract_successor.hpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief Declaration and definitions of abstract successor
 * @version 0.1
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma  once

#include <spot/twa/twa.hh>
#include <string>
#include <vector>

namespace kofola {
    /// abstract class for emptiness check procedure to obtain successors in different scenarios
    class abstract_successor
    {
    public:
        /// abstract class for the macrostate used by the abstract_successor
        class mstate { // {{{
        protected:
            spot::acc_cond::mark_t acc_; /// acc. marks on transition incoming to the state
            bdd trans_cond_; /// transition condition to come to this state
            bool encountered_ = false; /// to avoid incoming transition to the root
        public: // METHODS
            spot::acc_cond::mark_t get_acc() {return acc_; }
            void set_acc(spot::acc_cond::mark_t new_acc) {acc_ = new_acc; }

            void set_encountered(bool val) {encountered_ = val;}
            bool get_encountered() {return encountered_;}

            virtual bool eq(const mstate& rhs) const = 0;

            virtual bool lt(const mstate& rhs) const = 0;

            mstate& operator=(const mstate& other) = default;

            virtual ~mstate() { }
        }; // mstate }}}

    protected: // DATA MEMBERS

    public: // METHODS
        virtual std::vector<std::shared_ptr<mstate>> get_initial_states() = 0;

        virtual std::vector<std::shared_ptr<mstate>> get_succs(const std::shared_ptr<abstract_successor::mstate> &src) = 0;

        virtual bool is_accepting(spot::acc_cond::mark_t cond) = 0;

        virtual bool subsum_less_early(const std::shared_ptr<abstract_successor::mstate> a, const std::shared_ptr<abstract_successor::mstate> b) {return a->eq(*b);};

        virtual bool subsum_less_early_plus(const std::shared_ptr<abstract_successor::mstate> a, const std::shared_ptr<abstract_successor::mstate> b) {return a->eq(*b);};

        virtual void print_mstate(const std::shared_ptr<abstract_successor::mstate> a) = 0;
    }; // abstract_complement_alg }}}

    bool operator==(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs);

    bool operator!=(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs);

    bool operator<(const abstract_successor::mstate& lhs,
                   const abstract_successor::mstate& rhs);

} // kofola


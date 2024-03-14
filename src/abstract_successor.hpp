#pragma  once

#include <spot/twa/twa.hh>
#include <string>
#include <vector>

namespace kofola {
    class abstract_successor
    {
    public:
        class mstate { // {{{
        protected:
            spot::acc_cond::mark_t acc_;
            bdd trans_cond_;
        public: // METHODS
            spot::acc_cond::mark_t get_acc() {return acc_; }
            spot::acc_cond::mark_t set_acc(spot::acc_cond::mark_t new_acc) {acc_ = new_acc; }

            virtual bool eq(const mstate& rhs) const = 0;

            virtual bool neq(const mstate& rhs) const = 0;

            virtual bool lt(const mstate& rhs) const = 0;

            virtual ~mstate() { }
        }; // mstate }}}

    protected: // DATA MEMBERS

    public: // METHODS
        virtual std::vector<std::unique_ptr<mstate>> get_initial_states() = 0;

        virtual std::vector<std::unique_ptr<mstate>> get_succs(const std::unique_ptr<abstract_successor::mstate> &src) = 0;

        virtual bool is_accepting(spot::acc_cond::mark_t cond) = 0;
    }; // abstract_complement_alg }}}

    bool operator==(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs);

    bool operator!=(const abstract_successor::mstate& lhs,
                    const abstract_successor::mstate& rhs);

    bool operator<(const abstract_successor::mstate& lhs,
                   const abstract_successor::mstate& rhs);

} // kofola


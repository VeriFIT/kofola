// implementation of NCSB-based complementation algorithm for deterministic SCCs

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of NCSB-based complementation algorithm for deterministic SCCs
    class complement_subs_tuple : public abstract_complement_alg
    { // {{{
    public: // METHODS

        /// constructor
        complement_subs_tuple(const cmpl_info& info, unsigned part_index);

        virtual mstate_set get_init() override;

        virtual mstate_col_set get_succ_track(
                const std::set<unsigned>&  glob_reached,
                const mstate*              src,
                const bdd&                 symbol) override;

        virtual mstate_set lift_track_to_active(const mstate* src) override;

        virtual mstate_col_set get_succ_active(
                const std::set<unsigned>&  glob_reached,
                const mstate*              src,
                const bdd&                 symbol,
                bool resample = true) override;

        virtual bool use_round_robin() const override { return false; }

        virtual bool use_shared_breakpoint() const override { return this->info_.shared_breakpoint_; }

        virtual spot::acc_cond get_acc_cond() override
        { return spot::acc_cond(1, spot::acc_cond::inf({0})); }

        virtual unsigned get_min_colour() const override { return 0; }

        virtual ~complement_subs_tuple() override;

        mstate_col_set upper_succ(
                const std::set<unsigned>&  glob_reached,
                const mstate*              src,
                const bdd&                 symbol,
                bool track);

        mstate_col_set lower_succ(
                const std::set<unsigned>&  glob_reached,
                const mstate*              src,
                const bdd&                 symbol,
                bool track);

        //mstate get_upper_succ(const mstate* src);

        //mstate get_lower_succ(const mstate* src);
    }; // complement_subs_tuple }}}
} // namespace kofola }}}

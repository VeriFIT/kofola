// Copyright (C) 2022  The Kofola Authors
//
// Kofola is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Kofola is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// implementation of rank-based complementation

#pragma once

#include "abstract_complement_alg.hpp"

namespace kofola { // {{{

/// implementation of rank-based complementation algorithm nondeterministic
/// accepting SCCs
class complement_rank2 : public abstract_complement_alg
{ // {{{
public: // pimpl

	class impl;
	std::unique_ptr<impl> pimpl_;

public: // METHODS

	/// constructor
	complement_rank2(const cmpl_info& info, unsigned part_index);

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

	virtual bool use_round_robin() const override { return true; }

	virtual bool use_shared_breakpoint() const override { return false; }

	virtual unsigned get_min_colour() const override { return 0; }

	virtual spot::acc_cond get_acc_cond() override
	{ return spot::acc_cond(1, spot::acc_cond::inf({0})); }

	virtual ~complement_rank2() override;
}; // complement_rank2 }}}
} // namespace kofola }}}

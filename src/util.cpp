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

#include "util.hpp"

namespace
{ // {{{

/// returns true iff all output transitions of 'current_state' that are in SCC
/// 'scc' are accepting
bool all_out_trans_in_scc_acc(
	const spot::const_twa_graph_ptr&  aut,
	unsigned                          current_state,
	unsigned                          scc,
	const spot::scc_info&             si)
{ // {{{
	(void)scc;
	auto current_scc = si.scc_of(current_state);

	for (auto &t : aut->out(current_state)) {
		if (si.scc_of(t.dst) == current_scc && !t.acc) {
			return false;
		}
	}

	return true;
} // all_out_trans_in_scc_acc() }}}

} //  }}}


spot::twa_graph_ptr kofola::saturate(
	const spot::const_twa_graph_ptr&  aut,
	const spot::scc_info&             si)
{
	spot::twa_graph_ptr aut_new = spot::make_twa_graph(aut, spot::twa::prop_set::all());
	for (unsigned i = 0; i < si.scc_count(); i++) {
		bool changed = true;
		while (changed) {
			changed = false;
			for (auto state : si.states_of(i)) {
				if (all_out_trans_in_scc_acc(aut_new, state, i, si)) {
					for (auto s : si.states_of(i)) {
						for (auto &t : aut_new->out(s)) {
							if (t.dst == state && !t.acc) {
								t.acc = spot::acc_cond::mark_t{0};
								changed = true;
							}
						}
					}
				}
			}
		}
	}

	return aut_new;
}

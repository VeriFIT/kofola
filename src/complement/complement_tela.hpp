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

#pragma once

// kofola
#include "../util/helpers.hpp"

// spot
#include <spot/twa/twa.hh>
#include <spot/twaalgos/isdet.hh>


namespace kofola
{ // {{{

/// complements a transition-based Emerson-Lei automaton (TELA)
spot::twa_graph_ptr complement_tela(const spot::twa_graph_ptr& aut);

/// complements a deterministic automaton by making it complete and complementing acceptance
spot::twa_graph_ptr complement_deterministic(const spot::twa_graph_ptr& aut);

/// applies postprocessing to an automaton based on configuration parameters
spot::twa_graph_ptr apply_postprocessing(const spot::twa_graph_ptr& aut, const spot::twa_graph_ptr& original_aut);

/// complements a TELA using the synchronous algorithm (cf. paper)
spot::twa_graph_ptr complement_sync(const spot::twa_graph_ptr& aut);

/**
 * @brief Check if an automaton is suitable for reduction.
 *
 * The automaton is considered unsuitable if the number of atomic
 * propositions (APs) exceeds a configured threshold (currently 10) or
 * if the automaton is deterministic. These heuristics avoid expensive
 * or pointless reduction attempts on large or already-deterministic
 * automata.
 *
 * @param aut The input transition-based Emerson-Lei automaton (TELA) to check.
 * @return true if the automaton is suitable for reduction, false otherwise.
 */
bool is_reduction_suitable(const spot::twa_graph_ptr& aut);

/**
 * @brief Check whether an automaton is suitable for post-reduction processing.
 *
 * This performs heuristics used after an initial reduction to determine
 * if further reduction steps or specific post-processing should be applied.
 * The exact criteria mirror the reduction suitability checks but are
 * evaluated in the post-reduction context.
 *
 * @param aut The input transition-based Emerson-Lei automaton (TELA) to check.
 * @return true if the automaton is suitable for post-reduction processing, false otherwise.
 */
bool is_post_reduction_suitable(const spot::twa_graph_ptr& aut);

} // kofola }}}

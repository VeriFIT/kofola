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
#include "kofola.hpp"

// spot
#include <spot/twa/twa.hh>


namespace kofola
{ // {{{

/// complements a transition-based Emerson-Lei automaton (TELA)
spot::twa_graph_ptr complement_tela(
	const spot::twa_graph_ptr&  aut,
	const kofola::options&      options);

/// complements a TELA using the synchronous algorithm (cf. paper)
spot::twa_graph_ptr complement_sync(
	const spot::twa_graph_ptr&  aut,
	const kofola::options&      options);

} // kofola }}}

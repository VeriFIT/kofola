#pragma once

// spot
#include <spot/twa/twa.hh>

namespace kofola
{
    // Determinization of transition-based Emerson-Lei automata (TELA).
    // Template/interface only: no determinization logic yet.
    class tela_determinize
    {
    private:
        spot::twa_graph_ptr aut_;

    public:
        explicit tela_determinize(const spot::twa_graph_ptr& aut);

        /// Run determinization; currently a stub returning the input automaton.
        spot::twa_graph_ptr run_new();
    };
}

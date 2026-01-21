#pragma once

// spot
#include <spot/twa/twa.hh>
#include <spot/twaalgos/sccinfo.hh>

#include <memory>

namespace kofola {
    /// Determinize a transition-based Emerson-Lei automaton (TELA).
    ///
    /// This is currently a thin wrapper around `tela_determinize`.
    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut);

    // Determinization of transition-based Emerson-Lei automata (TELA).
    // Template/interface only: no determinization logic yet.
    class tela_determinize
    {
    private:
        spot::twa_graph_ptr aut_;
        std::shared_ptr<spot::scc_info> scc_;

    public:
        tela_determinize(const spot::twa_graph_ptr& aut,
                         std::shared_ptr<spot::scc_info> scc);

        /// Run determinization; currently a stub returning the input automaton.
        spot::twa_graph_ptr run_new();
    };
}

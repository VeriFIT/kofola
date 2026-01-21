#include "determinize_tela.hpp"

namespace kofola
{
    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut) {
        tela_determinize det(aut);
        return det.run_new();
    }

    tela_determinize::tela_determinize(const spot::twa_graph_ptr& aut)
        : aut_(aut) {
    }

    spot::twa_graph_ptr tela_determinize::run_new() {
        // Interface-only stub: determinization not implemented yet.
        return aut_;
    }
}

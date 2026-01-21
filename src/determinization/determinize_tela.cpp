#include "determinize_tela.hpp"

namespace kofola
{
    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut) {
        auto scc = std::make_shared<spot::scc_info>(aut);
        tela_determinize det(aut, std::move(scc));
        return det.run_new();
    }

    tela_determinize::tela_determinize(const spot::twa_graph_ptr& aut,
                                       std::shared_ptr<spot::scc_info> scc)
        : aut_(aut), scc_(std::move(scc)) {
    }

    spot::twa_graph_ptr tela_determinize::run_new() {
        // Interface-only stub: determinization not implemented yet.
        return aut_;
    }
}

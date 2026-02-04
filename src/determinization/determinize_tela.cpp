#include "determinize_tela.hpp"

#include "../complement/complement_sync.hpp"

namespace kofola
{
    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut) {
        auto scc = std::make_shared<spot::scc_info>(aut, spot::scc_info_options::ALL);

        // SCC partitioning (via the same helper used by complementation).
        // Return tuple items:
        //  (0) `num_partitions`     - number of created partitions
        //  (1) `part_to_type_map`   - partition index -> PartitionType
        //  (2) `st_to_part_map`     - state index -> partition index
        //  (3) `scc_to_part_map`    - SCC index -> partition index
        //  (4) `part_to_acc_map`    - partition index -> restricted acceptance condition
        auto partitions = helpers::tnba_complement::create_partitions(*scc, kofola::OPTIONS);

        tela_determinize det(aut, std::move(scc), std::move(partitions));
        return det.run_new();
    }

    tela_determinize::tela_determinize(const spot::twa_graph_ptr& aut,
                                       std::shared_ptr<spot::scc_info> scc,
                                       kofola::scc_partitions_t partitions)
        : aut_(aut), scc_(std::move(scc)), partitions_(std::move(partitions)) {
    }

    spot::twa_graph_ptr tela_determinize::run_new() {
        // Interface-only stub: determinization not implemented yet.
        return aut_;
    }
}

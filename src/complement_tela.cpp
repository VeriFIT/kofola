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

// kofola
#include "complement_tela.hpp"
#include "util.hpp"
#include "decomposer.hpp"

// Spot
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complete.hh>

// standard library
#include <queue>

spot::twa_graph_ptr kofola::complement_tela(const spot::twa_graph_ptr& aut)
{
	spot::twa_graph_ptr aut_reduced;
	std::vector<bdd> implications;
	spot::twa_graph_ptr aut_tmp = nullptr;
	if (aut_tmp)
		aut_reduced = aut_tmp;
	else
		aut_reduced = aut;

	spot::scc_info scc(aut_reduced, spot::scc_info_options::ALL);

	if (kofola::has_value("postponed", "yes", kofola::OPTIONS.params)) { // postponed procedure
		// saturation
		if (kofola::has_value("saturate", "yes", kofola::OPTIONS.params)) {
			aut_reduced = kofola::saturate(aut_reduced, scc);
			spot::scc_info scc_sat(aut_reduced, spot::scc_info_options::ALL);
			scc = scc_sat;
		}

		// decompose source automaton - TODO: this should be done properly
		cola::decomposer decomp(aut_reduced);
		auto decomposed = decomp.run(
			true,
			kofola::has_value("merge_iwa", "yes", kofola::OPTIONS.params),
			kofola::has_value("merge_det", "yes", kofola::OPTIONS.params));

		if (decomposed.size() > 0) {
			std::vector<spot::twa_graph_ptr> part_res;

			spot::postprocessor p_pre;
			p_pre.set_type(spot::postprocessor::Buchi);
			p_pre.set_level(spot::postprocessor::High);

			spot::postprocessor p_post;
			p_post.set_type(spot::postprocessor::Generic);
			if (kofola::has_value("low_red_interm", "yes", kofola::OPTIONS.params)) {
				p_post.set_level(spot::postprocessor::Low);
			} else {
				p_post.set_level(spot::postprocessor::High);
			}

			// comparison for priority queue - smallest automata should be at top
			auto aut_cmp = [](const auto& lhs, const auto& rhs){ return lhs->num_states() > rhs->num_states();};
			std::priority_queue<spot::twa_graph_ptr,
				std::vector<spot::twa_graph_ptr>, decltype(aut_cmp)> aut_queue(aut_cmp);
			for (auto aut : decomposed)
			{
				// if (decomp_options.scc_compl_high)
				//   p.set_level(spot::postprocessor::High);
				// else
				//   p.set_level(spot::postprocessor::Low);
				// complement each automaton
				auto aut_preprocessed = p_pre.run(aut);
				spot::scc_info part_scc(aut_preprocessed, spot::scc_info_options::ALL);

				auto res = kofola::complement_sync(aut_preprocessed);
				// postprocessing for each automaton
				// part_res.push_back(p_post.run(dec_aut));
				aut_queue.push(p_post.run(res));
			}

			assert(!aut_queue.empty());
			while (aut_queue.size() > 1) { // until single aut remains
				auto first_aut = aut_queue.top();
				aut_queue.pop();
				DEBUG_PRINT_LN("first_aut size = " + std::to_string(first_aut->num_states()));
				auto second_aut = aut_queue.top();
				aut_queue.pop();
				DEBUG_PRINT_LN("second_aut size = " + std::to_string(second_aut->num_states()));
				auto result = spot::product(first_aut, second_aut);
				result = p_post.run(result);
				aut_queue.push(result);
			}

			return aut_queue.top();
		} else {
			assert(false);
		}
	}

	// make sure the input is a BA
	spot::postprocessor p;
	p.set_type(spot::postprocessor::Buchi);
	p.set_level(spot::postprocessor::High);
	spot::twa_graph_ptr aut_to_compl;
	aut_to_compl = p.run(aut_reduced);

	auto res = kofola::complement_sync(aut_to_compl);
	DEBUG_PRINT_LN("finished call to run_new()");

	// postprocessing  TODO: should also consider other options
	if (!kofola::has_value("raw", "yes", kofola::OPTIONS.params)) {
		spot::postprocessor p_post;
		if ("buchi" == kofola::OPTIONS.output_type) {
			p_post.set_type(spot::postprocessor::Buchi);
		}
        else if("tgba" == kofola::OPTIONS.output_type) {
            p_post.set_type(spot::postprocessor::GeneralizedBuchi);
        } else {
			p_post.set_type(spot::postprocessor::Generic);
		}

		p_post.set_level(spot::postprocessor::Low);
		res = p_post.run(res);
	}

	return res;
}

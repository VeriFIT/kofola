#include "sets.hpp"

#include <deque>
#include <utility>
#include <vector>

namespace kofola {

std::vector<partition> nondet_split_set(const std::set<unsigned>& universum, unsigned num_sets) {
	std::vector<partition> partitions{};
	if (num_sets == 0) {
		if (universum.empty()) {
			partitions.push_back(partition{});
		}
		return partitions;
	}

	std::deque<std::pair<partition, unsigned>> queue;
	queue.push_back({partition(num_sets), 0});

	std::vector<unsigned> universum_vec(universum.begin(), universum.end());

	while (!queue.empty()) {
		auto [macrostate, state_idx] = queue.front();
		queue.pop_front();

		if (state_idx >= universum_vec.size()) {
			partitions.push_back(std::move(macrostate));
			continue;
		}

		for (unsigned i = 0; i < num_sets; i++) {
			partition new_macrostate{macrostate};
			new_macrostate[i].insert(universum_vec[state_idx]);
			queue.push_back({std::move(new_macrostate), state_idx + 1});
		}
	}

	return partitions;
}

} // namespace kofola

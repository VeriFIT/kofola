
#pragma once

#include <functional>
#include <limits>
#include <vector>
#include <set>

namespace kofola { // {{{

using partition = std::vector<std::set<unsigned>>;

/**
 * Split the input set `universum` into `num_sets` parts in a
 * nondeterministic manner. The function returns a vector of possible
 * partitions; each partition is a `partition` (i.e. a `std::vector<std::set<unsigned>>`)
 * containing exactly `num_sets` sets whose union equals `universum`.
 *
 * Parameters:
 * - universum: the set to split.
 * - num_sets: the number of parts in each partition.
 *
 * Returns:
 * - a `std::vector<partition>` with candidate partitions. The exact ordering
 *   and contents are not specified and may vary between calls.
 */
std::vector<partition> nondet_split_set(const std::set<unsigned>& universum, unsigned num_sets);

/**
 * Compute the Cartesian product of `xs` and `ys`, mapping each pair `(x,y)`
 * to an `R` via the provided `make` function. If either input vector is
 * empty, the result is an empty vector.
 *
 * Template parameters:
 * - T: element type of the input vectors.
 * - R: result element type.
 *
 * Parameters:
 * - xs: first input vector.
 * - ys: second input vector.
 * - make: function that produces an `R` from `(const T&, const T&)`.
 *
 * Returns:
 * - `std::vector<R>` containing `make(x, y)` for every `x` in `xs` and every `y` in `ys`.
 */
template <class T, class R>
std::vector<R> cartesian_product(const std::vector<T>& xs, const std::vector<T>& ys, const std::function<R(const T&, const T&)>& make) {
	std::vector<R> out;
	if (xs.empty() || ys.empty()) {
		return out;
	}

	const auto max = std::numeric_limits<std::size_t>::max();
	if (ys.size() != 0 && xs.size() <= max / ys.size()) {
		out.reserve(xs.size() * ys.size());
	}

	for (const auto& x : xs) {
		for (const auto& y : ys) {
			out.push_back(make(x, y));
		}
	}
	return out;
}

} // namespace kofola }}}
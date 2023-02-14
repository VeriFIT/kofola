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
#include "complement_alg_rank2.hpp"
#include "ranking.hpp"
#include "util.hpp"

using namespace kofola;
using mstate_set = abstract_complement_alg::mstate_set;
using mstate_col_set = abstract_complement_alg::mstate_col_set;

namespace { // {{{

/// representation of all other runs (outside the partition block)
const unsigned BOX = UINT_MAX;

/// partial macrostate for the given component
class mstate_rank : public abstract_complement_alg::mstate
{ // {{{
private: // DATA MEMBERS

	bool is_waiting_;                 // true - waiting, false - tight
	std::set<unsigned> states_;       // {} when !is_waiting_
	std::set<unsigned> breakpoint_;   // {} when !active_ || is_waiting_
	kofola::ranking f_;               // {} when is_waiting_
	int i_;                           // -1 when !active_ || is_waiting_
	bool active_;

private: // METHODS

	/// constructor - private; use create_* functions to create new mstate
	mstate_rank(
		bool                       is_waiting,
		const std::set<unsigned>&  states,
		const std::set<unsigned>&  breakpoint,
		const ranking&             f,
		int                        i,
		bool                       active
	) :
		is_waiting_(is_waiting),
		states_(states),
		breakpoint_(breakpoint),
		f_(f),
		i_(i),
		active_(active)
	{ }

public: // METHODS

	virtual std::string to_string() const override;
	virtual bool is_active() const override { return this->active_; }
	virtual bool eq(const mstate& rhs) const override;
	virtual bool lt(const mstate& rhs) const override;
	virtual ~mstate_rank() override { }

	virtual const std::set<unsigned>& get_breakpoint() const override { return this->breakpoint_; }
  	virtual void set_breakpoint(const std::set<unsigned>& breakpoint) override { this->breakpoint_ = breakpoint; }

	/// checks whether internal invariants hold
	bool invariants_hold() const;

	static mstate_rank create_waiting_ms(
		const std::set<unsigned>&  states,
		bool                       active_);

	static mstate_rank create_tight_active_ms(
		const std::set<unsigned>&  breakpoint,
		const ranking&             f,
		int                        i);

	static mstate_rank create_tight_passive_ms(const ranking& f);

	friend class kofola::complement_rank2::impl;
}; // mstate_rank }}}


bool mstate_rank::invariants_hold() const
{ // {{{
  if (!this->is_waiting_ && !this->states_.empty()) {
    return false;
  }
  if ((!this->active_ || this->is_waiting_ ) && !this->breakpoint_.empty()) {
    return false;
  }
  if (this->is_waiting_ && !this->f_.empty()) {
    return false;
  }
  if ((!this->active_ || this->is_waiting_ ) && this->i_ != -1) {
    return false;
  }

  return true;
} // invariants_hold() }}}

std::string mstate_rank::to_string() const
{ // {{{
  std::string res = std::string("[RANK(") + ((this->active_)? "A" : "T") + "): ";
  res += "S=";
  if (this->is_waiting_) { // WAITING
    res += std::to_string(this->states_);
  } else { // TIGHT
    std::set<int> dom_of_f;
    for (const auto& pair : this->f_) {
      dom_of_f.insert(pair.first);
    }
    res += std::to_string(dom_of_f);
    if (this->active_) {
      res += ", O=" + std::to_string(this->breakpoint_);
    }
    res += ", f=" + std::to_string(this->f_);
    if (this->active_) {
      res += ", i=" + std::to_string(this->i_);
    }
  }
  res += "]";
  return res;
} // to_string() }}}


bool mstate_rank::eq(const mstate& rhs) const
{ // {{{
  const mstate_rank* rhs_rank = dynamic_cast<const mstate_rank*>(&rhs);
  assert(rhs_rank);

  return this->active_ == rhs_rank->active_ &&
    this->is_waiting_ == rhs_rank->is_waiting_ &&
    this->states_ == rhs_rank->states_ &&
    this->breakpoint_ == rhs_rank->breakpoint_ &&
    this->f_ == rhs_rank->f_ &&
    this->i_ == rhs_rank->i_;
} // mstate_rank::eq() }}}

bool mstate_rank::lt(const mstate& rhs) const
{ // {{{
  const mstate_rank* rhs_rank = dynamic_cast<const mstate_rank*>(&rhs);
  assert(rhs_rank);

  if (this->active_ == rhs_rank->active_) {
    if (this->is_waiting_ == rhs_rank->is_waiting_) {
      if (this->i_ == rhs_rank->i_) {
        if (this->states_ == rhs_rank->states_) {
          if (this->breakpoint_ == rhs_rank->breakpoint_) {
            return this->f_ < rhs_rank->f_;
          } else {
            return this->breakpoint_ < rhs_rank->breakpoint_;
          }
        } else {
          return this->states_ < rhs_rank->states_;
        }
      } else {
        return this->i_ < rhs_rank->i_;
      }
    } else {
      return this->is_waiting_ < rhs_rank->is_waiting_;
    }
  } else {
    return this->active_ < rhs_rank->active_;
  }
} // mstate_rank::lt() }}}

mstate_rank mstate_rank::create_waiting_ms(
	const std::set<unsigned>&  states,
	bool                       active)
{
	return mstate_rank(
		true,              // is it Waiting?
		states,            // reachable states (S)
		{},                // breakpoint (O)
		{},                // ranking (f)
		-1,                // index of tracked rank (i)
		active);           // active
}

mstate_rank mstate_rank::create_tight_active_ms(
	const std::set<unsigned>&  breakpoint,
	const ranking&             f,
	int                        i)
{
	return mstate_rank(
		false,             // is it Waiting?
		{},                // reachable states (S)
		breakpoint,        // breakpoint (O)
		f,                 // ranking (f)
		i,                 // index of tracked rank (i)
		true);             // active
}

mstate_rank mstate_rank::create_tight_passive_ms(const ranking& f)
{
	return mstate_rank(
		false,          // is it Waiting?
		{},             // reachable states (S)
		{},             // breakpoint (O)
		f,              // ranking (f)
		-1,             // index of tracked rank (i)
		false);         // active
}

} // anonymous namespace }}}


class kofola::complement_rank2::impl : public abstract_complement_alg
{ // {{{
private:  // DATA MEMBERS

	// TODO

private: // METHODS

	/// gets a vector of order-maximal tight rankings for a given macrostate
	std::vector<ranking> get_max_tight_rankings(const std::set<unsigned>& states);

	/// gets a vector of order-maximal tight rankings with the given rank (i.e.,
	/// the maximum rank of a state) for a given macrostate
	std::vector<ranking> get_max_tight_rankings_with_rank(
		const std::set<unsigned>&  states,
		unsigned                   rank,
		const ranking&             bounds);

	/// gets the rank bounds for each state in a macrostate
	ranking get_rank_bounds(const std::set<unsigned>& states);

	/// restricts states in a set ot the current partition
	std::set<unsigned> restr_states_to_part(const std::set<unsigned>& states) const;

	/// returns the maxrank successor (see the paper) if it exists or nothing
	std::optional<ranking> maxrank(
		const std::set<unsigned>&  glob_reached,
		const ranking&             f,
		const bdd&                 symbol);

public:  // METHODS

	/// constructor
	impl(const cmpl_info& info, unsigned part_index);

	virtual mstate_set get_init() override;

	virtual mstate_col_set get_succ_track(
		const std::set<unsigned>&  glob_reached,
		const mstate*              src,
		const bdd&                 symbol) override;

	virtual mstate_set lift_track_to_active(const mstate* src) override;

	virtual mstate_col_set get_succ_active(
		const std::set<unsigned>&  glob_reached,
		const mstate*              src,
		const bdd&                 symbol,
		bool resample = true) override;

	virtual bool use_round_robin() const override { return false; }

	virtual bool use_shared_breakpoint() const override { return false; }

	virtual spot::acc_cond get_acc_cond() override
	{ return spot::acc_cond(1, spot::acc_cond::inf({0})); }

	virtual unsigned get_min_colour() const override { return 0; }

	virtual ~impl() override { }
}; // complement_rank2::impl }}}


complement_rank2::impl::impl(const cmpl_info& info, unsigned part_index)
	: abstract_complement_alg(info, part_index)
{ }


mstate_set complement_rank2::impl::get_init()
{ // {{{
	DEBUG_PRINT_LN("init RANK for partition " + std::to_string(this->part_index_));
	std::set<unsigned> init_state;

	unsigned orig_init = this->info_.aut_->get_init_state_number();
	if (this->info_.st_to_part_map_.at(orig_init) == this->part_index_) {
		init_state.insert(orig_init);
	} else { // TODO: we assert that the partition block is reachable from initial
		// FIXME: the box should not be necessary in cases such as a NAC reachable
		// from an initial state via an acyclic path
		init_state.insert(BOX);
	}

	std::shared_ptr<mstate> ms(new mstate_rank(
		mstate_rank::create_waiting_ms(init_state, false)));
  mstate_set result = {ms};
  return result;
} // get_init() }}}



std::optional<ranking> complement_rank2::impl::maxrank(
	const std::set<unsigned>&  glob_reached,
	const ranking&             f,
	const bdd&                 symbol)
{ // {{{

	// TODO: use bounds!

	assert(false);

	return std::nullopt;
} // maxrank() }}}


mstate_col_set complement_rank2::impl::get_succ_track(
	const std::set<unsigned>&  glob_reached,
	const mstate*              src,
	const bdd&                 symbol)
{ // {{{
	const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
	assert(src_rank);
	assert(!src_rank->active_);
	assert(src_rank->invariants_hold());

	mstate_col_set result;
	DEBUG_PRINT_LN("getting tracking successor for " + std::to_string(*src) +
		" over symbol " + std::to_string(symbol));

	if (src_rank->is_waiting_) {
		std::set<unsigned> succ_states = this->restr_states_to_part(glob_reached);
		assert(false);

	} else { // tight
		std::optional<ranking> g = this->maxrank(glob_reached, src_rank->f_, symbol);
		assert(false);
	}


} // get_succ_track() }}}


std::set<unsigned> complement_rank2::impl::restr_states_to_part(
	const std::set<unsigned>& states) const
{ // {{{
	std::set<unsigned> res;
	for (unsigned st : states) {
		if (this->info_.st_to_part_map_.at(st) == this->part_index_) {
			res.insert(st);
		}
	}

	return res;
} // restr_states_to_part() }}}


ranking complement_rank2::impl::get_rank_bounds(const std::set<unsigned>& states)
{ // {{{
	assert(states.size() > 0);
	assert(states.size() != 1 || *states.begin() != BOX);   // BOX is not the only element

	unsigned max_rank = 2* states.size() - 1;    // TODO: refine

	ranking bounds;

	DEBUG_PRINT_LN("FIXME: get tighter bounds!");

	if (kofola::is_in(BOX, states)) {
		auto it_bool_pair = bounds.emplace(BOX, max_rank);
		assert(it_bool_pair.second);
		max_rank -= 2;
	}

	for (unsigned st : states) {
		if (BOX != st) {
			auto it_bool_pair = bounds.emplace(st, max_rank);
			assert(it_bool_pair.second);
		}
	}

	return bounds;
} // get_rank_bounds() }}}


std::vector<ranking> complement_rank2::impl::get_max_tight_rankings_with_rank(
	const std::set<unsigned>&  states,
	unsigned                   rank,
	const ranking&             bounds)
{ // {{{
	assert(rank % 2 == 1);               // rank should be odd

	// here, CORE will be the states that keep the tight ranking of a run, cf.
	// our CONCUR'21 paper; RANKED CORE are states of a CORE with assigned ranks
	// (possibly implicitly, e.g., by order in a vector)
	unsigned core_size = (rank+1) / 2;
	assert(core_size <= states.size());  // core_size should be reasonable

	std::vector<unsigned> vec_states(states.begin(), states.end());

	std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>>
		ranked_cores = kofola::partial_permutations_ext(vec_states, core_size);

	DEBUG_PRINT_LN("ranked cores: " + std::to_string(ranked_cores));

	std::vector<ranking> vec_rankings;
	for (const auto& core_rest_pair : ranked_cores) {
		const std::vector<unsigned>& r_core = core_rest_pair.first;
		const std::vector<unsigned>& rest = core_rest_pair.second;
		ranking core_ranking;
		for (size_t i = 0; i < r_core.size(); ++i) {
			core_ranking.insert({r_core[i], 2 * i + 1});
		}

		DEBUG_PRINT_LN("core ranking: " + std::to_string(core_ranking));
		if (!rest.empty())
		{
			assert(false);
		}

		vec_rankings.push_back(core_ranking);
	}

	return vec_rankings;
} // get_max_tight_rankings_with_rank() }}}


std::vector<ranking> complement_rank2::impl::get_max_tight_rankings(
	const std::set<unsigned>& states)
{ // {{{
	// get the upper bound on the ranks of each state
	ranking bounds = this->get_rank_bounds(states);
	assert(!bounds.empty());
	DEBUG_PRINT_LN("rank bounds: " + std::to_string(bounds));

	// get maximum possible rank
	auto it_max = std::max_element(bounds.begin(), bounds.end(),
		[](const auto& lhs, const auto& rhs){ return lhs.second < rhs.second; });
	assert(it_max != bounds.end());
	int max_rank = static_cast<int>(it_max->second);
	if (max_rank % 2 == 0) { --max_rank; }    // make max_rank odd (danger: can become negative)
	DEBUG_PRINT_LN("max rank: " + std::to_string(max_rank));

	std::vector<ranking> vec_rankings;
	for (; max_rank > 0; --max_rank) {
		std::vector<ranking> rankings_for_rank =
			this->get_max_tight_rankings_with_rank(states, max_rank, bounds);
		DEBUG_PRINT_LN("max tight rankings for rank " + std::to_string(max_rank) +
			": " + std::to_string(rankings_for_rank));
		vec_rankings.insert(vec_rankings.end(),
			rankings_for_rank.begin(), rankings_for_rank.end());
	}

	DEBUG_PRINT_LN("max tight rankings: " + std::to_string(vec_rankings));
	return vec_rankings;
} // get_max_tight_rankings() }}}


mstate_set complement_rank2::impl::lift_track_to_active(const mstate* src)
{ // {{{
	const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
	assert(src_rank);
	assert(!src_rank->active_);
	assert(src_rank->invariants_hold());

	DEBUG_PRINT_LN("lifting " + std::to_string(*src_rank));

	mstate_set result;
	if (src_rank->is_waiting_) { // src is from WAITING
		DEBUG_PRINT_LN("lifting WAITING state");

		std::shared_ptr<mstate_rank> src_cpy(new mstate_rank(*src_rank));
		src_cpy->active_ = true;
		result.push_back(src_cpy);          // one option is to stay in WAITING

		if (src_rank->states_.empty() || // no interesting successors
				(src_rank->states_.size() == 1 && BOX == *src_rank->states_.begin())) {
			return result;
		}

		std::vector<ranking> max_rankings = this->get_max_tight_rankings(src_rank->states_);

		DEBUG_PRINT_LN("max tight rankings: " + std::to_string(max_rankings));

		for (const auto& f : max_rankings) {
			std::set<unsigned> breakpoint = f.with_rank(0);

			std::shared_ptr<mstate> ms(new mstate_rank(mstate_rank::create_tight_active_ms(
				breakpoint,     // breakpoint (O)
				f,              // ranking (f)
				0)));           // index of tracked rank (i)
			result.push_back(ms);
		}

		// // and let's compute the successors that move to TIGHT
		// std::vector<std::tuple<int, int, bool>> r;
		// assert(false);
    //
		// // FIXME is this correct? to consider only states from this partition
		// // auto bound = rank_restr_[std::set<unsigned>(mstate.curr_reachable_.begin(), mstate.curr_reachable_.end())];
    //
		// std::set<unsigned> states_no_box;
		// for (unsigned state : src_rank->states_) {
		// 	if (BOX != state) {
		// 		states_no_box.insert(state);
		// 	}
		// }
    //
		// // TODO: add RankRestr
		// DEBUG_PRINT_LN("powerset: " + std::to_string(states_no_box))
		// auto bound = get_rank_bound(rank_restr_, states_no_box, this->info_);
    //
		// DEBUG_PRINT_LN("rank bound: " + std::to_string(bound))
    //
		// for (auto s : src_rank->states_) {
		// 	bool accepting = ((s != BOX)? this->info_.state_accepting_[s] : false);
		// 	r.push_back(std::make_tuple(s, bound, accepting));
		// }
    //
		// std::vector<ranking> rankings = get_tight_rankings(r);
		// rankings = get_max(rankings);
		// for (auto rnking : rankings)
		// {
		// 	std::set<unsigned> breakpoint;
		// 	for (auto pr : rnking) { // construct breakpoint
		// 		if (pr.second == 0) { breakpoint.insert(pr.first); }
		// 	}
		// 	std::shared_ptr<mstate> ms(new mstate_rank(
		// 		{},                // reachable states (S)
		// 		false,             // is it Waiting?
		// 		breakpoint,        // breakpoint (O)
		// 		rnking,            // ranking (f)
		// 		0,                 // index of tracked rank (i)
		// 		true));            // active
		// 	result.push_back(ms);
		// }
	} else { // src is from TIGHT
		DEBUG_PRINT_LN("lifting TIGHT state");
		int new_i = 0;
		std::set<unsigned> breakpoint = src_rank->f_.with_rank(new_i);

		std::shared_ptr<mstate> ms(new mstate_rank(mstate_rank::create_tight_active_ms(
			breakpoint,        // breakpoint (O)
			src_rank->f_,      // ranking (f)
			new_i)));          // index of tracked rank (i)
		result.push_back(ms);
	}

	DEBUG_PRINT_LN("complement_rank::lift returning " + std::to_string(result));
	return result;
} // lift_track_to_active() }}}


mstate_col_set complement_rank2::impl::get_succ_active(
	const std::set<unsigned>&  glob_reached,
	const mstate*              src,
	const bdd&                 symbol,
	bool resample)
{ // {{{
	const mstate_rank* src_rank = dynamic_cast<const mstate_rank*>(src);
	assert(src_rank);
	assert(src_rank->active_);
	assert(src_rank->invariants_hold());

	DEBUG_PRINT_LN("getting active successor for " + std::to_string(*src) +
		" over symbol " + std::to_string(symbol));

	std::unique_ptr<mstate_rank> tmp_track;
	if (src_rank->is_waiting_) {
		tmp_track.reset(new mstate_rank(mstate_rank::create_waiting_ms(src_rank->states_, false)));
	} else {
		tmp_track.reset(new mstate_rank(mstate_rank::create_tight_passive_ms(src_rank->f_)));
	}
	mstate_col_set track_succ = this->get_succ_track(glob_reached, tmp_track.get(), symbol);

	assert(false);
} // get_succ_active() }}}


//  ********************* complement_rank2 **************************

complement_rank2::complement_rank2(const cmpl_info& info, unsigned part_index) :
	abstract_complement_alg(info, part_index),
	pimpl_(new impl(info, part_index))
{ }


mstate_set complement_rank2::get_init() { return this->pimpl_->get_init(); }


mstate_col_set complement_rank2::get_succ_track(
	const std::set<unsigned>&  glob_reached,
	const mstate*              src,
	const bdd&                 symbol)
{ return this->pimpl_->get_succ_track(glob_reached, src, symbol); }


mstate_set complement_rank2::lift_track_to_active(const mstate* src)
{ return this->pimpl_->lift_track_to_active(src); }


mstate_col_set complement_rank2::get_succ_active(
	const std::set<unsigned>&  glob_reached,
	const mstate*              src,
	const bdd&                 symbol,
	bool resample)
{ return this->pimpl_->get_succ_active(glob_reached, src, symbol, resample); }


complement_rank2::~complement_rank2() { }

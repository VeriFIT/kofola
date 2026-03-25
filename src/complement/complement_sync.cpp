// Copyright (C) 2017-2019 Laboratoire de Recherche et Développement
// de l'Epita.
// Copyright (C) 2022  The COLA Authors
//
// COLA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// COLA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// kofola
#include "../util/helpers.hpp"
#include "../types/types.hpp"
#include "complement_tela.hpp"
#include "decomposer.hpp"
#include "../util/util.hpp"
#include "complement_sync.hpp"
#include "elevatorization.hpp"

#include "../algorithms/abstract_complement_alg.hpp"
#include "complement_alg_mh.hpp"
#include "complement_alg_ncsb.hpp"
#include "complement_alg_ncsb_delay.hpp"
#include "complement_alg_safra.hpp"
// #include "complement_alg_rank.hpp"
#include "complement_alg_rank2.hpp"
#include "complement_alg_init_det.hpp"
#include "complement_alg_subs_tuple.hpp"
#include "complement_alg_sd_tela.hpp"
#include "complement_alg_sd_inductive.hpp"
#include "complement_alg_iadacs.hpp"

#include <deque>
#include <map>
#include <set>
#include <stack>
#include <queue>

#include <spot/misc/hashfunc.hh>
#include <spot/twaalgos/dot.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/sccinfo.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/degen.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/parity.hh>
#include <spot/twaalgos/cleanacc.hh>
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/misc/bddlt.hh>
#include <spot/parseaut/public.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/misc/version.hh>
#include <spot/twa/acc.hh>


// Complementation of Buchi automara based on SCC decomposition
// We classify three types of SCCs in the input NBA:
// 1. inherently weak SCCs (IWCs): every cycle in the SCC will not visit accepting transitions or every cycle visits an accepting transition
// 2. deterministic accepting SCCs (DACs): states in the SCC have at most one successor remain in the same SCC for a letter
// 3. nondeterministic accepting SCCs (NACs): has an accepting transition and nondeterministic

namespace helpers {

    tnba_complement::tnba_complement(const spot::twa_graph_ptr &aut, spot::scc_info &si)
            : aut_(aut),
              non_const_aut_(aut),
              dir_sim_(),
              reachable_vector_(),
              si_(si),
              info_(nullptr),
              alg_vec_(),
              num_colours_(0),
              final_code_(),
              used_infs_(),
              vec_acc_code_(),
              part_col_offset_(),
              sink_acc_code_(),
              partitions_(),
              nb_states_(aut->num_states()),
              res_(nullptr),
              support_(nb_states_),
              compat_(nb_states_),
              is_accepting_(aut->num_states(), false),
              scc_types_(),
              show_names_(),
              uberstate_to_num_map_(),
              num_to_uberstate_map_() {
        // Generate bdd supports and compatible options for each state.
        // Also check if all its transitions are accepting.
        for (unsigned i = 0; i < nb_states_; ++i) {
            bdd res_support = bddtrue;
            bdd res_compat = bddfalse;
            bool accepting = true;
            bool has_transitions = false;
            for (const auto &out: aut->out(i)) {
                has_transitions = true;
                res_support &= bdd_support(out.cond);
                res_compat |= out.cond;
                if (!out.acc)
                    accepting = false;
            }
            support_[i] = res_support;
            compat_[i] = res_compat;
            is_accepting_[i] = accepting && has_transitions;
        }
        // obtain the types of each SCC
        scc_types_ = get_scc_types(si_);
        // find out the DACs and NACs

        if (kofola::LOG_VERBOSITY > 0) {
            DEBUG_PRINT_LN("Complementing the following aut:");
            spot::print_hoa(std::cerr, this->aut_);
            std::cerr << "\n\n\n\n";
        }

        // we don't use pruning of macrostates for large automata
        if(aut->num_states() >= 8000) {
            kofola::OPTIONS.params["sim-ms-prune"] = "no";
        }
        
        if(kofola::has_value("sim-ms-prune", "yes", kofola::OPTIONS.params)) {
            // Perform reduction and compute simulation on the current automaton.
            // This may change the automaton (state count, edges, etc.).
            this->reduce_and_compute_simulation();

            // The automaton potentially changed; refresh dependent structures.
            // 1) Update the number of states to match the reduced automaton.
            this->nb_states_ = this->aut_->num_states();
            // 2) Recompute SCC information for the reduced automaton.
            this->si_ = spot::scc_info(this->aut_, spot::scc_info_options::ALL);
            // For TELA automata, Spot's is_accepting_scc may report "unknown" for
            // Fin-based acceptance conditions unless determine_unknown_acceptance() is
            // called.  Without this, create_partitions() sees all SCCs as non-accepting
            // and returns num_partitions == 0, causing select_algorithms() to do nothing.
            if (kofola::has_value("tela", "yes", kofola::OPTIONS.params)) {
                this->si_.determine_unknown_acceptance();
            }
            // 3) Recompute SCC types.
            this->scc_types_ = get_scc_types(this->si_);
            // 4) Resize per-state vectors to the new size; they will be filled below.
            this->support_.assign(this->nb_states_, bddtrue);
            this->compat_.assign(this->nb_states_, bddfalse);
            this->is_accepting_.assign(this->nb_states_, false);
        }

        // Show macrostate labels in state names when raw=yes or show-macrostate-labels=yes.
        // With raw=yes postprocessing is skipped so labels are preserved in HOA output.
        this->show_names_ = kofola::has_value("raw", "yes", kofola::OPTIONS.params) ||
                            kofola::has_value("show-macrostate-labels", "yes", kofola::OPTIONS.params);

        // compute vector of accepting states, supports, etc.
        for (unsigned i = 0; i < this->aut_->num_states(); ++i) {
            bdd res_support = bddtrue;
            bdd res_compat = bddfalse;
            bool accepting = true;
            bool has_transitions = false;
            for (const auto &out: this->aut_->out(i)) {
                has_transitions = true;
                res_support &= bdd_support(out.cond);
                res_compat |= out.cond;
                if (!out.acc) {
                    accepting = false;
                }
            }
            this->support_[i] = res_support;
            this->compat_[i] = res_compat;
            this->is_accepting_[i] = accepting && has_transitions;
        }

        // Build reachability vector on the (possibly reduced) automaton.
        // Do not assume any particular SCC numbering direction.
        if(kofola::has_value("sim-ms-prune", "yes", kofola::OPTIONS.params)) {
            std::vector<std::set<int>> aux_reach = this->get_reachable_vector();
            for (unsigned st = 0; st < this->aut_->num_states(); ++st) {
                const std::set<int> &old_set = aux_reach[st];
                std::set<unsigned> new_set;
                for (const auto &el: old_set) {
                    assert(el >= 0);
                    new_set.insert(static_cast<unsigned>(el));
                }

                this->reachable_vector_.push_back(new_set);
            }
        }
        partitions_ = create_partitions(this->si_, kofola::OPTIONS);
        info_ = std::make_unique<kofola::cmpl_info>(
                this->aut_,             // automaton
                std::get<0>(partitions_),         // number of partitions
                std::get<1>(partitions_),       // partition types
                std::get<2>(partitions_),            // state to partition map
                this->reachable_vector_,// vector of reachable states
                create_part_to_scc_map(std::get<3>(partitions_)),        // map of partitions to sets of SCCS they contain
                create_scc_to_pred_sccs_map(this->si_, this->reachable_vector_),   // maps SCCs to the sets of their predecessors
                std::get<4>(partitions_),            // partitions to acceptance condition map
                this->si_,              // SCC information
                this->dir_sim_,         // direct simulation
                this->is_accepting_,    // vector for acceptance of states
                kofola::has_value("sh-break", "yes", kofola::OPTIONS.params));
    }

    unsigned
    tnba_complement::get_num_states() {
        return this->nb_states_;
    }

    spot::scc_info &
    tnba_complement::get_scc_info() {
        return this->si_;
    }

    void tnba_complement::reduce_and_compute_simulation() {
        // compute simulation on current automaton (may change aut_)
        std::vector<bdd> implications;
        this->non_const_aut_ = spot::simulation(this->aut_, &implications, -1);
        this->aut_ = this->non_const_aut_;

        // Recompute SCC info on the reduced automaton to evaluate reachability correctly
        spot::scc_info tmp_si(this->aut_);

        // get vector of simulated states
        std::vector<std::vector<char>> implies(
                implications.size(),
                std::vector<char>(implications.size(), 0));
        for (unsigned i = 0; i != implications.size(); ++i) {
            if (!tmp_si.reachable_state(i))
                continue;
            for (unsigned j = 0; j != implications.size(); ++j) {
                // reachable states
                if (!tmp_si.reachable_state(j))
                    continue;
                // j simulates i
                bool i_implies_j = bdd_implies(implications[i], implications[j]);
                if (i_implies_j) {
                    dir_sim_.push_back({i, j});
                }
            }
        }
    }

    std::set<int> tnba_complement::reachable_vertices(std::vector<std::vector<int>> list, std::set<int> from) {
        std::set<int> all(from);
        std::stack<int> stack;
        int item;
        for (int it: all)
            stack.push(it);

        while (stack.size() > 0) {
            item = stack.top();
            stack.pop();
            for (int dst: list[item]) {
                if (all.find(dst) == all.end()) {
                    stack.push(dst);
                    all.insert(dst);
                }
            }
        }
        return all;
    }

    std::vector<std::set<int>> tnba_complement::get_reachable_vector() {
        std::vector<std::set<int>> reachable_vector;

        std::vector<std::set<int>> list_set(aut_->num_states());
        std::vector<std::vector<int>> list_vector(aut_->num_states());

        for (unsigned s = 0; s < aut_->num_states(); s++) {
            reachable_vector.push_back(std::set<int>());
            list_set[s] = std::set<int>();

            // iterate over all transitions from s
            for (const auto &t: aut_->out(s)) {
                list_set[s].insert(t.dst);
            }

            list_vector[s] = std::vector<int>(list_set[s].begin(), list_set[s].end());
        }

        for (unsigned int s = 0; s < aut_->num_states(); s++) {
            std::set<int> tmp({static_cast<int>(s)});
            reachable_vector[s] = reachable_vertices(list_vector, tmp);
        }

        return reachable_vector;
    }

    // ######################################################################
    // NEW INTERFACE
    // ######################################################################

    /// constructor
    tnba_complement::uberstate::uberstate(const std::set<unsigned> &reached_states,
                                          const vec_macrostates &part_macrostates,
                                          int active_scc,
                                          const std::set<unsigned> &shared_breakpoint) :
            reached_states_(reached_states),
            part_macrostates_(part_macrostates),
            active_scc_(active_scc),
            shared_breakpoint_(shared_breakpoint) {}

    /// deleted copy constructor and assignment operator
    tnba_complement::uberstate::uberstate(const uberstate &us) :
            reached_states_(us.reached_states_),
            part_macrostates_(us.part_macrostates_),
            active_scc_(us.active_scc_),
            shared_breakpoint_(us.shared_breakpoint_) {}

    /// converts to string
    std::string tnba_complement::uberstate::to_string() const { // {{{
        std::string result;
        result += "<" + std::to_string(this->reached_states_);
        result += " |" + std::to_string(this->active_scc_) + "| ";
        result += std::to_string(this->shared_breakpoint_) + "| ";

        for (size_t i = 0; i < this->part_macrostates_.size(); ++i) {
            result += "c" + std::to_string(i) + ": " +
                      this->part_macrostates_[i]->to_string();
            if (this->part_macrostates_.size() != i + 1) {
                result += ", ";
            }
        }

        result += ">";
        return result;
    } // to_string() }}}

    /// returns the set of all reached states
    const std::set<unsigned> &tnba_complement::uberstate::get_reach_set() const { return this->reached_states_; }

    /// returns the partial macrostates
    const tnba_complement::vec_macrostates &
    tnba_complement::uberstate::get_part_macrostates() const { return this->part_macrostates_; }

    /// returns the index of the active SCC (INACTIVE_SCC if no active)
    int tnba_complement::uberstate::get_active_scc() const { return this->active_scc_; }

    /// get shared breakpoint
    const std::set<unsigned> &
    tnba_complement::uberstate::get_shared_breakpoint() const { return this->shared_breakpoint_; }

    /// total ordering operator to allow use in std::set and std::map
    bool tnba_complement::uberstate::operator<(const uberstate &rhs) const { // {{{
        assert(this->part_macrostates_.size() == rhs.part_macrostates_.size());

        // let's start by comparing active SCC indices
        if (this->active_scc_ != rhs.active_scc_) {
            return this->active_scc_ < rhs.active_scc_;
        }

        // continue by comparing reached_states_
        // FIXME: inefficient - in C++20, we would use the flying saucer
        // operator <=>
        if (this->reached_states_ != rhs.reached_states_) {
            return this->reached_states_ < rhs.reached_states_;
        }

        // then, compare the partial macrostates lexicographically
        const size_t length = this->part_macrostates_.size();
        for (size_t i = 0; i < length; ++i) {
            if (*(this->part_macrostates_[i]) != *(rhs.part_macrostates_[i])) {
                return *(this->part_macrostates_[i]) < *(rhs.part_macrostates_[i]);
            }
        }

        // compare according to shared breakpoint
        if (this->shared_breakpoint_ != rhs.shared_breakpoint_) {
            return this->shared_breakpoint_ < rhs.shared_breakpoint_;
        }

        // they are equal
        return false;
    } // operator< }}}

    bool tnba_complement::uberstate::operator==(const uberstate &rhs) const { // {{{
        assert(this->part_macrostates_.size() == rhs.part_macrostates_.size());

        if (this->active_scc_ != rhs.active_scc_ ||
            this->reached_states_ != rhs.reached_states_) {
            return false;
        }

        // then, compare the partial macrostates
        const size_t length = this->part_macrostates_.size();
        for (size_t i = 0; i < length; ++i) {
            if (*(this->part_macrostates_[i]) != *(rhs.part_macrostates_[i])) {
                return false;
            }
        }

        // compare according to shared breakpoint
        if (this->shared_breakpoint_ != rhs.shared_breakpoint_) {
            return false;
        }

        return true;
    } // operator== }}}


    // Here we have a bidirectional map between uberstates and state
    // identifiers (unsigned).  The uberstates are physically stored only at
    // 'num_to_uberstate_map_', the reason being that they contain vectors of
    // unique_ptr (no copy is therefore allowed).

    /// maps uberstates to state numbers
    std::map<const helpers::tnba_complement::uberstate *, unsigned, helpers::tnba_complement::uberstate_ptr_less_ftor> uberstate_to_num_map_;
    /// maps state numbers to uberstates
    std::vector<std::shared_ptr<helpers::tnba_complement::uberstate>> num_to_uberstate_map_;
    /// counter of states (to be assigned to uberstates) - 0 is reserved for sink
    unsigned cnt_state_ = 0;

    // reserved colours
    // static const unsigned SINK_COLOUR = 0;
    enum {
        SINK_COLOUR = 0
    };

    /// accessor into the uberstate table
    unsigned helpers::tnba_complement::uberstate_to_num(const helpers::tnba_complement::uberstate &us) const { // {{{
        auto it = this->uberstate_to_num_map_.find(&us);
        if (this->uberstate_to_num_map_.end() != it) {
            return it->second;
        } else {
            assert(false);
            return 0;
        }
    } // uberstate_to_num() }}}

    /// translates state number to uberstate
    const helpers::tnba_complement::uberstate &helpers::tnba_complement::num_to_uberstate(unsigned num) const { // {{{
        assert(num < this->num_to_uberstate_map_.size());
        assert(this->num_to_uberstate_map_[num]);
        return *this->num_to_uberstate_map_[num];
    } // num_to_uberstate() }}}

    std::pair<bdd, bdd> helpers::tnba_complement::uberstate_support(unsigned uberstateId) const {
        if(get_is_sink_created() && get_sink_state() == uberstateId) {
            return {bddfalse, bddfalse};
        }
        const auto &us = num_to_uberstate(uberstateId);
        bdd support = bddtrue;
        bdd compat = bddfalse;
        for (unsigned st: us.get_reach_set()) {
            support &= support_[st];
            compat |= compat_[st];
        }
        return {support, compat};
    }

    /// inserts an uberstate (by moving) and returns its assigned number (if
    /// not present), or just returns the number of an equal uberstate (if
    /// present)
    unsigned helpers::tnba_complement::insert_uberstate(const helpers::tnba_complement::uberstate &us) { // {{{
        DEBUG_PRINT_LN("inserting uberstate " + us.to_string());
        auto it = this->uberstate_to_num_map_.find(&us);
        if (this->uberstate_to_num_map_.end() == it) { // not found
            DEBUG_PRINT_LN("not found!");
            std::shared_ptr<uberstate> ptr(new uberstate(us));
            this->num_to_uberstate_map_.push_back(ptr);
            assert(this->num_to_uberstate_map_.size() == this->cnt_state_ + 1);  // invariant
            DEBUG_PRINT_LN("inserting at position " + std::to_string(this->cnt_state_));
            const uberstate *us_new = this->num_to_uberstate_map_[this->cnt_state_].get();
            assert(*ptr == *us_new);
            auto jt_bool_pair = this->uberstate_to_num_map_.insert({us_new, this->cnt_state_});
            //assert(jt_bool_pair.second);    // insertion happened
            ++this->cnt_state_;
            DEBUG_PRINT_LN("inserted as " + std::to_string(jt_bool_pair.first->second));
            return jt_bool_pair.first->second;
        } else { // found
            DEBUG_PRINT_LN("found as " + std::to_string(it->second));
            return it->second;
        }
    } // insert_uberstate() }}}


    int helpers::tnba_complement::get_next_active_scc(const helpers::tnba_complement::vec_algorithms &alg_vec, int prev) { // {{{
        for (size_t i = prev + 1; i < alg_vec.size(); ++i) {
            if (alg_vec[i]->use_round_robin()) {
                return i;
            }
        }

        for (int i = 0; i < std::min(prev + 1, static_cast<int>(alg_vec.size())); ++i) {
            if (alg_vec[i]->use_round_robin()) {
                return i;
            }
        }

        return INACTIVE_SCC;
    } // get_next_active_scc() }}}


    kofola::PartitionToSCCMap helpers::tnba_complement::create_part_to_scc_map(
            const kofola::SCCToPartitionMap &scc_to_part_map) { // {{{
        kofola::PartitionToSCCMap part_to_scc_map;
        for (const auto &scc_part_pair: scc_to_part_map) {
            if (scc_part_pair.second != -1) {
                unsigned part = scc_part_pair.second;
                auto it_bool_pair = part_to_scc_map.insert({part, {scc_part_pair.first}});
                if (!it_bool_pair.second) { // no insertion
                    it_bool_pair.first->second.insert(scc_part_pair.first);
                }
            }
        }
        return part_to_scc_map;
    } // create_part_to_scc_map() }}}


    kofola::SCCToSCCSetMap helpers::tnba_complement::create_scc_to_pred_sccs_map(
            const spot::scc_info &si,
            const kofola::ReachableVector &reach_vec) { // {{{
        DEBUG_PRINT_LN("in create_scc_to_pred_sccs_map");
        DEBUG_PRINT_LN("reach_vec = " + std::to_string(reach_vec));
        for (unsigned state = 0; state < si.get_aut()->num_states(); ++state) {
            DEBUG_PRINT_LN("SCC of state " + std::to_string(state) + " = " +
                           std::to_string(si.scc_of(state)));
        }

        kofola::SCCToSCCSetMap scc_to_pred_sccs_map;

        for (size_t st = 0; st < reach_vec.size(); ++st) {
            unsigned st_scc_index = si.scc_of(st);
            if (static_cast<int>(st_scc_index) < 0) { // SCC unreachable
                DEBUG_PRINT_LN("wrong scc");
                continue;
            }
            for (unsigned dst: reach_vec[st_scc_index]) {
                unsigned dst_scc_index = si.scc_of(dst);
                auto it_bool_pair = scc_to_pred_sccs_map.insert(
                        {dst_scc_index, {st_scc_index}});
                if (!it_bool_pair.second) { // if no insertion happened
                    it_bool_pair.first->second.insert(st_scc_index);
                }
            }
        }

        DEBUG_PRINT_LN("scc_to_pred_sccs_map: " + std::to_string(scc_to_pred_sccs_map));
        return scc_to_pred_sccs_map;
    } // create_scc_to_pred_sccs_map() }}}

    bool helpers::tnba_complement::subsum_less_early(unsigned a, unsigned b) {
        // If any operand is the sink state, avoid dereferencing the placeholder
        // nullptr stored in num_to_uberstate_map_ for the sink. Only the sink
        // trivially subsumes itself; otherwise, treat it as incomparable so we
        // do not prune incorrectly.
        if (get_is_sink_created() && (a == sink_state_ || b == sink_state_)) {
            return a == b;
        }

        auto uber_a = num_to_uberstate(a);
        auto uber_b = num_to_uberstate(b);

        return uber_a.subsum_less_early(uber_b);
    }

    bool helpers::tnba_complement::subsum_less_early_plus(unsigned a, unsigned b) {
        if (get_is_sink_created() && (a == sink_state_ || b == sink_state_)) {
            return a == b;
        }

        auto uber_a = num_to_uberstate(a);
        auto uber_b = num_to_uberstate(b);

        return uber_a.subsum_less_early_plus(uber_b);
    }

    /// gets all successors of an uberstate wrt a vector of algorithms and a
    /// symbol
    helpers::tnba_complement::vec_state_taggedcol helpers::tnba_complement::get_succ_uberstates(
            const helpers::tnba_complement::uberstate &src,
            const bdd &symbol) { // {{{
        DEBUG_PRINT_LN("Processing uberstate " + std::to_string(src) +
                       " for symbol " + std::to_string(symbol));
        using mstate = kofola::abstract_complement_alg::mstate;
        using mstate_set = kofola::abstract_complement_alg::mstate_set;
        using mstate_taggedcol = std::pair<std::shared_ptr<mstate>, std::set<std::pair<unsigned, unsigned>>>;
        using mstate_col_set = kofola::abstract_complement_alg::mstate_col_set;
        using mstate_taggedcol_set = std::vector<mstate_taggedcol>;

        assert(alg_vec_.size() == src.get_part_macrostates().size());
        const int active_index = src.get_active_scc();
        std::set<unsigned> all_succ = kofola::get_all_successors(
                this->aut_, src.get_reach_set(), symbol);

        DEBUG_PRINT_LN("all succ over " + std::to_string(symbol) + "= " + std::to_string(all_succ));
        if (kofola::has_value("sim-ms-prune", "yes", kofola::OPTIONS.params)) { // if doing simulation pruning
            std::set<unsigned> pruned_succ = all_succ;

            std::set<unsigned> to_remove;
            for (const auto &pr: this->dir_sim_) {
                unsigned smaller = pr.first;
                unsigned bigger = pr.second;
                if (smaller == bigger ||  // identity
                    !kofola::is_in(smaller, all_succ) ||
                    !kofola::is_in(bigger, all_succ)
                        ) { // the pair is irrelevant
                    continue;
                }

                // if (kofola::is_in(bigger, this->reachable_vector_[smaller]) &&
                //     !kofola::is_in(smaller, this->reachable_vector_[bigger])) {
                if (this->si_.scc_of(smaller) > this->si_.scc_of(bigger)) {
                    // we assume that if SCC A is reachable from SCC B, then #(A) < #(B) (this is how spot seems to order SCCs)
                    DEBUG_PRINT_LN("removing " + std::to_string(smaller) +
                                   " because of " + std::to_string(bigger));
                    pruned_succ.erase(smaller);
                }
            }
            DEBUG_PRINT_LN("pruned succ = " + std::to_string(pruned_succ));
            all_succ = pruned_succ;
        }


        const vec_macrostates &prev_part_macro = src.get_part_macrostates();
        // this container collects all sets of pairs of macrostates and colours,
        // later, we will turn it into the Cartesian product
        std::vector<mstate_taggedcol_set> succ_part_macro_col;
        const std::set<unsigned> &sh_break = src.get_shared_breakpoint();
        DEBUG_PRINT_LN("Shared breakpoint " + std::to_string(sh_break));
        for (size_t i = 0; i < alg_vec_.size(); ++i) {
            const kofola::abstract_complement_alg::mstate *ms = prev_part_macro[i].get();
            mstate_col_set mcs;

            // for shared breakpoint we need to propagete the breakpoint to each partial macrostate
            if (alg_vec_[i]->use_shared_breakpoint()) {
                const_cast<mstate *>(ms)->set_breakpoint(sh_break);
            }

            if (alg_vec_[i]->use_shared_breakpoint()) {
                mcs = alg_vec_[i]->get_succ_active(all_succ, ms, symbol, sh_break.empty());
            } else if (active_index == static_cast<int>(i) || !alg_vec_[i]->use_round_robin()) {
                mcs = alg_vec_[i]->get_succ_active(all_succ, ms, symbol);
            } else {
                mcs = alg_vec_[i]->get_succ_track(all_succ, ms, symbol);
            }

            if (mcs.empty()) { // one empty set of successor macrostates
                return {};
            }

            if (alg_vec_[i]->use_shared_breakpoint()) {
                const_cast<mstate *>(ms)->clear_breakpoint();
            }

            // tagged colours
            mstate_taggedcol_set tagged_mcs;
            for (const auto &mstate_col: mcs) {
                std::set<std::pair<unsigned, unsigned>> new_taggedcols;
                for (const auto &col: mstate_col.second) {
                    new_taggedcols.insert({i, col});
                }
                tagged_mcs.push_back({mstate_col.first, new_taggedcols});
            }

            succ_part_macro_col.emplace_back(std::move(tagged_mcs));
        }

        DEBUG_PRINT_LN("generated partial macrostates + colours: " +
                       std::to_string(succ_part_macro_col));

        // compute the Cartesian product of the partial macrostates (+ colours)
        std::vector<std::vector<mstate_taggedcol>> cp =
                compute_cartesian_prod(succ_part_macro_col);

        // generate uberstates
        vec_state_taggedcol result;
        for (const auto &vec: cp) {
            vec_macrostates vm;
            std::set<unsigned> breakpoint;
            std::set<std::pair<unsigned, unsigned>> cols;

            assert(alg_vec_.size() == vec.size());
            for (size_t i = 0; i < vec.size(); i++) {
                if (alg_vec_[i]->use_shared_breakpoint()) {
                    const std::set<unsigned> &vb = vec[i].first->get_breakpoint();
                    breakpoint.insert(vb.begin(), vb.end());
                    vec[i].first->clear_breakpoint();
                    assert(vec[i].first->get_breakpoint().size() == 0);

                    // We generate colors only if we go from the macrostate with empty shared breakpoint
                    if (sh_break.empty()) {
                        cols.insert(vec[i].second.begin(), vec[i].second.end());
                    }
                } else {
                    cols.insert(vec[i].second.begin(), vec[i].second.end());
                }
                vm.push_back(vec[i].first);
            }

            // TODO: so far behaves fishy when both round robing and shared breakpoint is active for an algorithm
            unsigned us_num = UINT_MAX;   // canary value
            if (INACTIVE_SCC == active_index) { // no round robin
                DEBUG_PRINT_LN("inserting")
                us_num = insert_uberstate(uberstate(all_succ, vm, INACTIVE_SCC, breakpoint));
                result.emplace_back(us_num, std::move(cols));
                DEBUG_PRINT_LN("inserted")
            } else { // round robin
                if (vm[active_index]->is_active()) { // the same SCC active
                    us_num = insert_uberstate(uberstate(all_succ, vm, active_index, breakpoint));
                    result.emplace_back(us_num, cols);
                } else { // another SCC active
                    int next_active = get_next_active_scc(alg_vec_, active_index);
                    DEBUG_PRINT_LN("next active index: " + std::to_string(next_active));
                    assert(INACTIVE_SCC != next_active);

                    // now we need to lift the macrostate for next_active from track to active
                    mstate_set active_macros = alg_vec_[next_active]->lift_track_to_active(vm[next_active].get());
                    for (const auto &am: active_macros) {
                        vm[next_active] = am;
                        us_num = insert_uberstate(uberstate(all_succ, vm, next_active, breakpoint));
                        result.emplace_back(us_num, cols);
                    }
                }
            }
        }

        DEBUG_PRINT_LN("computed successors: " + std::to_string(result));

        remove_duplicit(result);

        return result;
    } // get_succ_uberstates() }}}

    /// gets all initial uberstates wrt a vector of algorithms
    std::vector<unsigned>
    helpers::tnba_complement::get_initial_uberstates() { // {{{
        std::set<unsigned> initial_states = {aut_->get_init_state_number()};

        int init_active = get_next_active_scc(alg_vec_, INACTIVE_SCC);
        DEBUG_PRINT_LN("initial active partition: " + std::to_string(init_active));

        using mstate_set = kofola::abstract_complement_alg::mstate_set;
        std::vector<mstate_set> vec_mstate_sets;
        for (size_t i = 0; i < alg_vec_.size(); ++i) { // get outputs of all procedures
            mstate_set init_mstates = alg_vec_[i]->get_init();

            DEBUG_PRINT_LN("obtained initial state for algs " + std::to_string(i) +
                           ": " + std::to_string(init_mstates));
            if (static_cast<int>(i) == init_active || !alg_vec_[i]->use_round_robin()) { // make the partial macrostate active
                mstate_set new_mstates;
                for (const auto &st: init_mstates) {
                    mstate_set lifted = alg_vec_[i]->lift_track_to_active(st.get());
                    new_mstates.insert(new_mstates.end(), lifted.begin(), lifted.end());
                }
                remove_duplicit(new_mstates);
                init_mstates = std::move(new_mstates);
            }
            if (init_mstates.empty()) { return {}; }   // one empty set terminates
            vec_mstate_sets.emplace_back(std::move(init_mstates));
        }

        DEBUG_PRINT_LN("obtained vector of sets of partial macrostates: " + std::to_string(vec_mstate_sets));

        // compute the cartesian product from the vector of sets of macrostates
        std::vector<vec_macrostates> cp = compute_cartesian_prod(vec_mstate_sets);

        DEBUG_PRINT_LN("initial macrostates: " + std::to_string(cp));

        std::vector<unsigned> result;
        for (const auto &vec: cp) {
            // collect shared breakpoint from partial algorithms supporting it
            std::set<unsigned> sh_break;
            for (size_t i = 0; i < vec.size(); i++) {
                if (alg_vec_[i]->use_shared_breakpoint()) {
                    const std::set<unsigned> &vb = vec[i]->get_breakpoint();
                    sh_break.insert(vb.begin(), vb.end());
                    vec[i]->clear_breakpoint();
                }
            }

            unsigned us_num = insert_uberstate(uberstate(initial_states, vec, init_active, sh_break));
            result.push_back(us_num);
        }

        // sort and remove duplicates
        remove_duplicit(result);

        return result;
    } // get_initial_uberstates() }}}


    /**
     * Partitions the strongly connected components (SCCs) of the input automaton
     * according to decomposition options. This function is used in the complementation
     * algorithm to divide SCCs into groups (partitions) based on their acceptance type
     * (inherently weak, deterministic, nondeterministic) and user-specified options.
     *
     * @param scc_inf  SCC information of the automaton (spot::scc_info)
     * @param options  Decomposition and partitioning options (kofola::options)
     * @return A tuple containing:
     *         - num_partitions: number of partitions created
     *         - partition_types: map from partition index to PartitionType
     *         - state_to_partition_map: map from state index to partition index
     *         - scc_to_partition_map: map from SCC index to partition index
     *         - partition_to_acc_map: map from partition index to acceptance condition
     */
    std::tuple<size_t,
            kofola::PartitionToTypeMap,
            kofola::StateToPartitionMap,
            kofola::SCCToPartitionMap,
            kofola::PartitionToAccMap
    >
    helpers::tnba_complement::create_partitions(
            const spot::scc_info &scc_inf,
            const kofola::options &options) { // {{{
        using kofola::PartitionType;

        std::string scc_types = helpers::get_scc_types(scc_inf);
        size_t part_index = 0;

        kofola::PartitionToTypeMap part_to_type_map;
        kofola::StateToPartitionMap st_to_part_map;
        kofola::SCCToPartitionMap scc_to_part_map;   // -1 is invalid partition
        kofola::PartitionToAccMap part_to_acc_map {};
        std::unordered_map<unsigned, spot::acc_cond::mark_t> part_to_used {};

        int iwa_index = -1;
        int dac_index = -1;
        int iadac_index = -1;

        bool merge_iwa = kofola::has_value("merge_iwa", "yes", options.params);
        bool merge_det = kofola::has_value("merge_det", "yes", options.params);

        // When merging IWAs, move the IWA partition to the front.  This makes it
        // be active first, and may avoid larger state space generation.
        if (merge_iwa) {
            DEBUG_PRINT_LN("Merge IWA");
            for (size_t i = 0; i < scc_inf.scc_count(); ++i) {
                if (helpers::is_accepting_weakscc(scc_types, i)) { // if there is some IWA
                    iwa_index = part_index;
                    ++part_index;
                    part_to_type_map[iwa_index] = PartitionType::INHERENTLY_WEAK;
                    break;
                }
            }
        }

        // similar thing as above for DACs
        if (merge_det) {
            DEBUG_PRINT_LN("Merge DET");
            for (size_t i = 0; i < scc_inf.scc_count(); ++i) {
                if (helpers::is_accepting_detscc(scc_types, i) && !helpers::is_accepting_initial_almost_detscc(scc_types, i)) { // if there is some DAC
                    dac_index = part_index;
                    ++part_index;
                    part_to_type_map[dac_index] = PartitionType::DETERMINISTIC;
                    break;
                }
            }
        }

        // decide the partitioning
        for (size_t i = 0; i < scc_inf.scc_count(); ++i) {
            DEBUG_PRINT_LN("Processing SCC " + std::to_string(i));
            DEBUG_PRINT_LN("scc_partition map: " + std::to_string(scc_to_part_map));
            if (!helpers::is_accepting_scc(scc_types, i)) {
                scc_to_part_map[i] = -1;
                continue; // we don't care about nonaccepting SCCs
            }

            scc_to_part_map[i] = part_index;
            if (helpers::is_accepting_weakscc(scc_types, i)) {
                DEBUG_PRINT_LN("SCC " + std::to_string(i) + " is IWA");
                if (merge_iwa) { // merging IWAs
                    if (-1 == iwa_index) {
                        iwa_index = part_index;
                        part_to_type_map[iwa_index] = PartitionType::INHERENTLY_WEAK;
                        ++part_index;
                    }
                    scc_to_part_map[i] = iwa_index;
                } else { // not merging IWAs
                    part_to_type_map[part_index] = PartitionType::INHERENTLY_WEAK;
                    scc_to_part_map[i] = part_index;
                    ++part_index;
                }
            } else if (kofola::OPTIONS.operation != "inclusion" && helpers::is_accepting_initial_almost_detscc(scc_types, i)) {
                // for inclusion we treat  SCCs as DACs --> we don't have 
                // emptiness checking for general TELA
                DEBUG_PRINT_LN("SCC " + std::to_string(i) + " is deterministic border nondeterministic");
                if(merge_det) {
                    if (-1 == iadac_index) {
                        iadac_index = part_index;
                        part_to_type_map[iadac_index] = PartitionType::INITIAL_ALMOST_DETERMINISTIC;
                        ++part_index;
                    }
                    scc_to_part_map[i] = iadac_index;
                } else {
                    part_to_type_map[part_index] = PartitionType::INITIAL_ALMOST_DETERMINISTIC;
                    ++part_index;
                }
            }else if (helpers::is_accepting_detscc(scc_types, i) || (kofola::OPTIONS.operation == "inclusion" && helpers::is_accepting_initial_almost_detscc(scc_types, i))) {
                DEBUG_PRINT_LN("SCC " + std::to_string(i) + " is DAC");
                if (merge_det) { // merging DACs
                    if (-1 == dac_index) {
                        dac_index = part_index;
                        part_to_type_map[dac_index] = PartitionType::DETERMINISTIC;
                        ++part_index;
                    }
                    scc_to_part_map[i] = dac_index;
                } else { // not merging DACs
                    part_to_type_map[part_index] = PartitionType::DETERMINISTIC;
                    scc_to_part_map[i] = part_index;
                    ++part_index;
                }
            } else if (helpers::is_accepting_nondetscc(scc_types, i)) {
                DEBUG_PRINT_LN("SCC " + std::to_string(i) + " is NAC");
                part_to_type_map[part_index] = PartitionType::NONDETERMINISTIC;
                ++part_index;
            } else {
                throw std::runtime_error("Invalid SCC on the input");
            }

            // partition index
            int partition_index = scc_to_part_map[i];
            // part_to_used gathers all acceptance conditions used in the union of SCCs (partition)
            if(part_to_used.find(partition_index) == part_to_used.end()) {
                part_to_used[partition_index] = scc_inf.acc_sets_of(i);
            } else {
                part_to_used[partition_index] |= scc_inf.acc_sets_of(i);
            }
        }

        // restrict the acceptance condition to colors occurring in the partition
        for(const auto& [part, acc] : part_to_used) {
            part_to_acc_map[part] = scc_inf.get_aut()->acc().restrict_to(acc);
        }

        DEBUG_PRINT_LN("scc_partition map: " + std::to_string(scc_to_part_map));
        // map states to correct partitions
        for (size_t i = 0; i < scc_inf.get_aut()->num_states(); ++i) {
            if (kofola::is_in(scc_inf.scc_of(i), scc_to_part_map)) {
                st_to_part_map[i] = scc_to_part_map.at((scc_inf.scc_of(i)));
            }
        }

        DEBUG_PRINT_LN("number of partitions: " + std::to_string(part_index));
        DEBUG_PRINT_LN("partition to type map: " + std::to_string(part_to_type_map));
        DEBUG_PRINT_LN("state_to_partition map: " + std::to_string(st_to_part_map));
        DEBUG_PRINT_LN("scc_partition map: " + std::to_string(scc_to_part_map));

        return {part_index, part_to_type_map, st_to_part_map, scc_to_part_map, part_to_acc_map};
    } // create_partitions() }}}


    /// selects the algorithms to run on the SCCs
    void helpers::tnba_complement::select_algorithms()  { // {{{
        using kofola::PartitionType;

        for (size_t i = 0; i < this->info_->num_partitions_; ++i) { // determine which algorithms to run on each of the SCCs
            helpers::tnba_complement::abs_cmpl_alg_p alg;
            
            const PartitionType partition_type = this->info_->part_to_type_map_.at(i);

            switch (partition_type) {
                case PartitionType::INHERENTLY_WEAK:
                    alg = create_inherently_weak_algorithm(i);
                    break;
                case PartitionType::DETERMINISTIC:
                    alg = create_deterministic_algorithm(i);
                    break;
                case PartitionType::STRONGLY_DETERMINISTIC:
                    alg = create_strongly_deterministic_algorithm(i);
                    break;
                case PartitionType::NONDETERMINISTIC:
                    alg = create_nondeterministic_algorithm(i);
                    break;
                case PartitionType::INITIAL_ALMOST_DETERMINISTIC:
                    alg = create_initial_almost_deterministic_algorithm(i);
                    break;
                default:
                    throw std::runtime_error("Strange SCC type found!");
            }
            
            alg_vec_.push_back(std::move(alg));
        }

    } // select_algorithms() }}}

    /**
     * Creates the complementation algorithm for an inherently weak partition.
     *
     * This function returns a unique pointer to a `complement_mh` algorithm instance
     * configured for the specified partition index. The inherently weak algorithm is
     * used for SCCs (strongly connected components) that are classified as inherently weak.
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for inherently weak SCCs.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_inherently_weak_algorithm(size_t partition_index) { // {{{
        return std::make_unique<kofola::complement_mh>(*(this->info_.get()), partition_index);
    } // create_inherently_weak_algorithm() }}}

    /**
     * Creates the complementation algorithm for a deterministic partition.
     *
     * This function selects and returns a unique pointer to the appropriate deterministic
     * complementation algorithm for the given partition index, based on the acceptance
     * condition (Büchi or generalized). It may use NCSB, NCSB-delay, or SD-TELA algorithms.
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for deterministic SCCs.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_deterministic_algorithm(size_t partition_index) { // {{{
        // take the acceptance condition simplified according to the partition
        bool is_buchi = this->info_->part_to_acc_map_.at(partition_index).is_buchi();
        
        if (is_buchi) {
            if (kofola::has_value("ncsb-delay", "yes", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_ncsb_delay>(*(this->info_.get()), partition_index);
            } else {
                return std::make_unique<kofola::complement_ncsb>(*(this->info_.get()), partition_index);
            }
        } else {
            if (kofola::has_value("tela_det_alg", "inductive", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_sd_inductive>(*(this->info_.get()), partition_index);
            } else if(kofola::has_value("tela_det_alg", "dnf-tela", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_sd_tela>(*(this->info_.get()), partition_index);
            }
            return std::make_unique<kofola::complement_sd_tela>(*(this->info_.get()), partition_index);
        }
    } // create_deterministic_algorithm() }}}

    /**
     * Creates the complementation algorithm for a strongly deterministic partition.
     *
     * This function returns a unique pointer to the appropriate strongly deterministic
     * complementation algorithm for the specified partition index, using either NCSB or
     * SD-TELA depending on the acceptance condition.
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for strongly deterministic SCCs.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_strongly_deterministic_algorithm(size_t partition_index) { // {{{
        // take the acceptance condition simplified according to the partition
        bool is_buchi = this->info_->part_to_acc_map_.at(partition_index).is_buchi();

        if (is_buchi) {
            return std::make_unique<kofola::complement_ncsb>(*(this->info_.get()), partition_index);
        } else {
            if (kofola::has_value("tela_det_alg", "inductive", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_sd_inductive>(*(this->info_.get()), partition_index);
            } else if(kofola::has_value("tela_det_alg", "dnf-tela", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_sd_tela>(*(this->info_.get()), partition_index);
            }
            return std::make_unique<kofola::complement_sd_tela>(*(this->info_.get()), partition_index);
        }
    } // create_strongly_deterministic_algorithm() }}}

    /**
     * Creates the complementation algorithm for a nondeterministic partition.
     *
     * This function selects and returns a unique pointer to the appropriate nondeterministic
     * complementation algorithm for the given partition index, based on user options and
     * acceptance condition. It may use subs_tup, rank, or Safra-based algorithms.
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for nondeterministic SCCs.
     * @throws std::runtime_error if the algorithm for general components is not implemented.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_nondeterministic_algorithm(size_t partition_index) { // {{{
        bool is_buchi = this->info_->part_to_acc_map_.at(partition_index).is_buchi();

        if(!is_buchi) {
            throw std::runtime_error("Algorithm for general components is not implemented");
        }

        if (kofola::has_value("nac-alg", "subs_tup", kofola::OPTIONS.params)) {
            // use subs_tup for NACs
            return std::make_unique<kofola::complement_subs_tuple>(*(this->info_.get()), partition_index);
        } else if (kofola::has_value("nac-alg", "rank", kofola::OPTIONS.params)) {
            // Rank algorithm temporarily disabled due to compilation issues
            // Fall back to Safra-based complementation
            return std::make_unique<kofola::complement_safra>(*(this->info_.get()), partition_index);
        } else {
            // use determinization-based
            return std::make_unique<kofola::complement_safra>(*(this->info_.get()), partition_index);
        }
    } // create_nondeterministic_algorithm() }}}

    /**
     * Creates the complementation algorithm for the initial deterministic partition.
     *
     * This function returns a unique pointer to the `complement_init_det` algorithm instance
     * configured for the specified partition index. It is used for the initial deterministic
     * component in the modular complementation procedure.
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for the initial deterministic SCC.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_initial_deterministic_algorithm(size_t partition_index) { // {{{
        // initial deterministic component
        bool is_buchi = this->info_->part_to_acc_map_.at(partition_index).is_buchi();
        if (is_buchi) {
            return std::make_unique<kofola::complement_init_det>(*(this->info_.get()), partition_index);
        } else {
            if (kofola::has_value("tela_det_alg", "inductive", kofola::OPTIONS.params)) {
                return std::make_unique<kofola::complement_sd_inductive>(*(this->info_.get()), partition_index);
            } else {
                return std::make_unique<kofola::complement_sd_tela>(*(this->info_.get()), partition_index);
            }
        }
    } // create_initial_deterministic_algorithm() }}}

    /**
     * Creates the complementation algorithm for a initial almost deterministic partition.
     *
     * This function returns a unique pointer to a `complement_init_almost_det` algorithm instance
     * configured for the specified partition index. The initial almost deterministic algorithm is
     * used for SCCs (strongly connected components) that are classified as initial almost deterministic (IADACs).
     *
     * @param partition_index Index of the partition for which the algorithm is created.
     * @return Unique pointer to the abstract complementation algorithm for IADACs.
     */
    helpers::tnba_complement::abs_cmpl_alg_p 
    helpers::tnba_complement::create_initial_almost_deterministic_algorithm(size_t partition_index) { // {{{
        return std::make_unique<kofola::complement_init_almost_det>(*(this->info_.get()), partition_index);
    } // create_initial_almost_deterministic_algorithm() }}}

    bdd helpers::tnba_complement::get_support_at(unsigned s) {
        return support_[s];
    }

    bdd helpers::tnba_complement::get_compat_at(unsigned s) {
        return compat_[s];
    }

    unsigned int helpers::tnba_complement::get_cnt_state_() {
        return cnt_state_;
    }

    void helpers::tnba_complement::inc_cnt_state_() {
        ++this->cnt_state_;
    }

    void helpers::tnba_complement::handle_sink_state()
    {
        if (!is_sink_created_) {
            is_sink_created_ = true;
            sink_state_ = this->cnt_state_;
            ++this->cnt_state_;
            this->num_to_uberstate_map_.push_back(nullptr);

            DEBUG_PRINT_LN("creating a sink state: " + std::to_string(sink_state_));
        }
    }

    bool helpers::tnba_complement::get_is_sink_created() const
    {
        return is_sink_created_;
    }

    unsigned helpers::tnba_complement::get_sink_state() const
    {
        return sink_state_;
    }

    /// new modular complementation procedure
    spot::twa_graph_ptr
    helpers::tnba_complement::run_new() { // {{{
        DEBUG_PRINT_LN("selecting algorithms");
        // creates a vector of algorithms, for every SCC of aut one
        select_algorithms();

        DEBUG_PRINT_LN("algorithms selected");

        // our structure for the automaton (TODO: hash table might be better)
        std::map<unsigned, std::vector<std::pair<bdd, vec_state_taggedcol>>> compl_states;

        std::stack<unsigned> todo;

        // get initial uberstates
        auto init_vec{this->get_initial_uberstates()};
        for (unsigned state: init_vec) {
            compl_states.insert({state, {}});
            todo.push(state);
        }

        DEBUG_PRINT_LN("initial todo: " + std::to_string(todo));

        while (!todo.empty()) { // the main loop
            // get next uberstate
            unsigned us_num = todo.top();
            todo.pop();
            const uberstate &us = num_to_uberstate(us_num);

            // get the post of 'us'
            auto it = compl_states.find(us_num);
            assert(compl_states.end() != it);
            std::vector<std::pair<bdd, vec_state_taggedcol>> &us_post = it->second;

            DEBUG_PRINT_LN("processing " + std::to_string(us_num) + ": " + us.to_string());

            // compute support of all available states
            // TODO: this should be cached for each reach_set
            bdd msupport = bddtrue;
            bdd n_s_compat = bddfalse;
            const std::set<unsigned> &reach_set = us.get_reach_set();

            // compute the occurred variables in the outgoing transitions of ms, stored in msupport
            for (unsigned s: reach_set) {
                msupport &= get_support_at(s);
                n_s_compat |= get_compat_at(s);
            }

            // direct non-support symbols to sink
            bdd all = n_s_compat;
            if (all != bddtrue) {
                if (!is_sink_created_) { // first time encountering sink state
                    handle_sink_state();
                    // create a sink state (its transitions)
                    compl_states.insert({sink_state_, {{bddtrue, {{sink_state_, {{UINT_MAX, SINK_COLOUR}}}}}}});
                }
                vec_state_taggedcol succs = {{sink_state_, {}}};
                us_post.emplace_back(std::make_pair(!all, std::move(succs)));
            }

            // iterate over all symbols
            while (all != bddfalse) {
                bdd letter = bdd_satoneset(all, msupport, bddfalse);
                all -= letter;

                DEBUG_PRINT_LN("symbol: " + std::to_string(letter));

                vec_state_taggedcol succs = this->get_succ_uberstates(us, letter);
                us_post.emplace_back(std::make_pair(letter, succs));

                for (const auto &state_cols: succs) {
                    const unsigned &succ_state = state_cols.first;

                    auto it_bool_pair = compl_states.insert({succ_state, {}});
                    if (it_bool_pair.second) { // the successor state is new
                        DEBUG_PRINT_LN("inserted " + std::to_string(succ_state) + " into compl_states_");
                        todo.push(succ_state);
                    }
                }
            }
        }

        bool new_init_created = false;
        if (init_vec.size() > 1) { // handle multiple initial states
            new_init_created = true;
            DEBUG_PRINT_LN("handling multiple initial states: " +
                           std::to_string(init_vec));
            unsigned new_init = this->cnt_state_;
            ++this->cnt_state_;
            DEBUG_PRINT_LN("new init state: " + std::to_string(new_init));
            auto it_bool_pair = compl_states.insert({new_init, {}});
            assert(it_bool_pair.second);
            auto &new_init_trans_vec = it_bool_pair.first->second;

            for (const auto &state: init_vec) {
                const auto &trans_vec = compl_states.at(state);
                DEBUG_PRINT_LN("state " + std::to_string(state) +
                               ": " + std::to_string(trans_vec));
                new_init_trans_vec.insert(new_init_trans_vec.end(),
                                          trans_vec.begin(), trans_vec.end());
            }
            init_vec = {new_init};
        }

        DEBUG_PRINT_LN(std::to_string(compl_states));

        set_acc_cond();
        // finish colors
        if (is_sink_created_) {                // accept also in sink, if created
            final_code_ |= sink_acc_code_;
        }

        spot::acc_cond result_cond(num_colours_, final_code_);

        // convert the result into a spot automaton
        // FIXME: we should be directly constructing spot aut
        spot::twa_graph_ptr result = spot::make_twa_graph(this->aut_->get_dict());
        result->copy_ap_of(this->aut_);
        result->prop_copy(this->aut_,
                          {
                                  false,        // state based
                                  false,        // inherently_weak
                                  false, false, // deterministic
                                  false,        // complete
                                  false         // stutter inv
                          });

        result->set_acceptance(result_cond);
        DEBUG_PRINT_LN("Acc = " + std::to_string(result->get_acceptance()));


        std::vector<std::string> *state_names = nullptr;
        if (show_names_) { // show names
            state_names = new std::vector<std::string>();
            result->set_named_prop("state-names", state_names);
        }
        for (const auto &st_trans_pair: compl_states) {
            const unsigned &src = st_trans_pair.first;
            unsigned spot_state = result->new_state();
            (void) spot_state;
            assert(spot_state == src);

            for (const auto &bdd_vec_tgt_pair: st_trans_pair.second) {
                const bdd &symbol = bdd_vec_tgt_pair.first;
                for (const auto &tgt_col_pair: bdd_vec_tgt_pair.second) {
                    const unsigned &tgt = tgt_col_pair.first;
                    const std::set<std::pair<unsigned, unsigned>> &cols = tgt_col_pair.second;
                    std::vector<unsigned> new_cols;
                    for (const std::pair<unsigned, unsigned> &part_col_pair: cols) {
                        const unsigned part_index = part_col_pair.first;
                        const unsigned colour = part_col_pair.second;

                        if (src == sink_state_) { // sink state
                            new_cols.push_back(SINK_COLOUR);
                            continue;
                        } else if (UINT_MAX == colour) { // FIXME: this is a special way
                            // of dealing with determinization-based (this should set the
                            // colour to the maximum colour) - try to make more uniform

                            unsigned max_col = vec_acc_code_[part_index].num_sets() - 1;
                            new_cols.push_back(part_col_offset_.at(part_index) + max_col);
                        } else { // standard transition
                            unsigned shift = alg_vec_[part_index]->get_min_colour(); // how much to decrement the colour
                            new_cols.push_back(part_col_offset_.at(part_index) + colour - shift);
                        }
                    }
                    spot::acc_cond::mark_t spot_cols(new_cols.begin(), new_cols.end());
                    result->new_edge(src, tgt, symbol, spot_cols);
                }
            }

            if (this->show_names_) { // handle output state names
                if (is_sink_created_ && sink_state_ == src) {
                    state_names->push_back("SINK");
                } else if (new_init_created && init_vec[0] == src) {
                    assert(init_vec.size() == 1);
                    state_names->push_back("INIT");
                } else {
                    state_names->push_back(num_to_uberstate(src).to_string());
                }
            }
        }

        assert(init_vec.size() == 1);
        result->set_init_state(init_vec[0]);

        if (kofola::LOG_VERBOSITY > 0) {
            spot::print_hoa(std::cerr, result);
            std::cerr << "\n\n\n\n";
        }

        return result;
    } // run_new() }}}

    std::vector<spot::acc_cond> helpers::tnba_complement::get_vec_acc_cond()
    {
        return vec_acc_code_;
    }

    spot::acc_cond::acc_code helpers::tnba_complement::get_final_acc_code()
    {
        return final_code_;
    }

    unsigned int helpers::tnba_complement::get_alg_vec_mincolour_at_i(unsigned i)
    {
        return alg_vec_[i]->get_min_colour();
    }

    std::map<unsigned int, unsigned int> helpers::tnba_complement::get_part_col_offset()
    {
        return part_col_offset_;
    }

    std::set<unsigned> helpers::tnba_complement::set_acc_cond() {
        num_colours_ = RESERVED_COLOURS;
        int rr_colour = -1;     // colour for round robin
        int sh_br_colour = -2;  // colour for shared breakpoint

        sink_acc_code_ = spot::acc_cond::acc_code::inf({SINK_COLOUR});
        spot::acc_cond::acc_code alg_acc_code = spot::acc_cond::acc_code::t();
        for (size_t i = 0; i < this->info_->num_partitions_; ++i) { // sum up acceptance conditions
            const auto &alg = alg_vec_[i];

            const spot::acc_cond &cond = alg->get_acc_cond();
            vec_acc_code_.push_back(cond);
            spot::acc_cond::acc_code cond_code = cond.get_acceptance();
            if (alg->use_shared_breakpoint()) {
                if (sh_br_colour < 0) {
                    sh_br_colour = num_colours_;
                    ++num_colours_;
                    used_infs_.insert(static_cast<unsigned>(sh_br_colour));
                    alg_acc_code &= spot::acc_cond::acc_code::inf({static_cast<unsigned>(sh_br_colour)});
                }

                part_col_offset_[i] = sh_br_colour;
            } else if (alg->use_round_robin()) {
                if (rr_colour < 0) {  // the first round robin
                    rr_colour = num_colours_;
                    ++num_colours_;
                    used_infs_.insert(static_cast<unsigned>(rr_colour));
                    alg_acc_code &= spot::acc_cond::acc_code::inf({static_cast<unsigned>(rr_colour)});
                }

                part_col_offset_[i] = rr_colour;
            } else {
                used_infs_.insert(i+1);
                cond_code <<= num_colours_;
                alg_acc_code &= cond_code;
                part_col_offset_[i] = num_colours_;
                num_colours_ += cond.num_sets();
            }
        }

        DEBUG_PRINT_LN("colour offsets: " + std::to_string(part_col_offset_));
        DEBUG_PRINT_LN("vec_acc_code: " + std::to_string(vec_acc_code_));

        final_code_ = alg_acc_code;
        DEBUG_PRINT_LN("final code: " + std::to_string(final_code_));
        // std::cout << std::to_string(final_code_) << "\n";
        return used_infs_;
    }
}

spot::twa_graph_ptr kofola::complement_sync(const spot::twa_graph_ptr& aut)
{
    spot::scc_info si(aut, spot::scc_info_options::ALL);

    // if we work with TELA, we need to properly determine SCC acceptance
    // Spot's is_acceptance might say unknown for Fin conditions
    if (kofola::has_value("tela", "yes", kofola::OPTIONS.params)) {
        si.determine_unknown_acceptance();
    }
    
    auto comp = helpers::tnba_complement(aut, si);
    auto res = comp.run_new();

    return res;
}
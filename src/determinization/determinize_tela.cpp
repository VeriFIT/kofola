#include "determinize_tela.hpp"

#include "../complement/complement_sync.hpp"
#include "../complement/complement_tela.hpp"
#include "../complement/elevatorization.hpp"
#include "determinize_alg_dac.hpp"
#include "determinize_alg_mh.hpp"

#include <stack>

#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/postproc.hh>

namespace kofola
{
    // ######################################################################
    // uberstate
    // ######################################################################

    tela_determinize::uberstate::uberstate(const std::set<unsigned>& reached_states,
                                           const vec_macrostates& part_macrostates) :
            reached_states_(reached_states),
            part_macrostates_(part_macrostates)
    { }

    std::string tela_determinize::uberstate::to_string() const { // {{{
        std::string result;
        result += "<" + std::to_string(this->reached_states_) + "| ";

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

    const std::set<unsigned>& tela_determinize::uberstate::get_reach_set() const
    { return this->reached_states_; }

    const tela_determinize::vec_macrostates&
    tela_determinize::uberstate::get_part_macrostates() const
    { return this->part_macrostates_; }

    /// Note that the partial macrostates are compared through the *virtual*
    /// mstate operators (i.e., mstate::lt / mstate::eq), not by pointer
    /// identity: two independently constructed partial macrostates that are
    /// semantically equal must collapse into a single uberstate, otherwise the
    /// construction would not terminate.  The vectors are aligned by partition
    /// index, so the i-th entries always belong to the same partial algorithm
    /// and the dynamic_casts inside lt/eq are guaranteed to succeed.
    bool tela_determinize::uberstate::operator<(const uberstate& rhs) const { // {{{
        assert(this->part_macrostates_.size() == rhs.part_macrostates_.size());

        if (this->reached_states_ != rhs.reached_states_) {
            return this->reached_states_ < rhs.reached_states_;
        }

        for (size_t i = 0; i < this->part_macrostates_.size(); ++i) {
            if (*(this->part_macrostates_[i]) != *(rhs.part_macrostates_[i])) {
                return *(this->part_macrostates_[i]) < *(rhs.part_macrostates_[i]);
            }
        }

        return false;   // they are equal
    } // operator< }}}

    bool tela_determinize::uberstate::operator==(const uberstate& rhs) const { // {{{
        assert(this->part_macrostates_.size() == rhs.part_macrostates_.size());

        if (this->reached_states_ != rhs.reached_states_) { return false; }

        for (size_t i = 0; i < this->part_macrostates_.size(); ++i) {
            if (*(this->part_macrostates_[i]) != *(rhs.part_macrostates_[i])) { return false; }
        }

        return true;
    } // operator== }}}

    // ######################################################################
    // tela_determinize
    // ######################################################################

    tela_determinize::tela_determinize(const spot::twa_graph_ptr& aut,
                                       std::shared_ptr<spot::scc_info> scc,
                                       kofola::scc_partitions_t partitions)
        : aut_(aut),
          scc_(std::move(scc)),
          partitions_(std::move(partitions)),
          info_(nullptr),
          alg_vec_(),
          part_to_scc_map_(),
          reachable_vector_(),
          scc_to_pred_sccs_map_(),
          dir_sim_(),
          support_(aut->num_states(), bddtrue),
          compat_(aut->num_states(), bddfalse),
          is_accepting_(aut->num_states(), false),
          num_colours_(0),
          final_code_(spot::acc_cond::acc_code::f()),
          vec_acc_code_(),
          part_col_offset_(),
          uberstate_to_num_map_(),
          num_to_uberstate_map_(),
          cnt_state_(0),
          show_names_(kofola::has_value("raw", "yes", kofola::OPTIONS.params) ||
                      kofola::has_value("show-macrostate-labels", "yes", kofola::OPTIONS.params))
    {
        // compute supports, compatible symbols and accepting states
        for (unsigned i = 0; i < this->aut_->num_states(); ++i) {
            bdd res_support = bddtrue;
            bdd res_compat = bddfalse;
            bool accepting = true;
            bool has_transitions = false;
            for (const auto& out : this->aut_->out(i)) {
                has_transitions = true;
                res_support &= bdd_support(out.cond);
                res_compat |= out.cond;
                if (!out.acc) { accepting = false; }
            }
            this->support_[i] = res_support;
            this->compat_[i] = res_compat;
            this->is_accepting_[i] = accepting && has_transitions;
        }

        // 'cmpl_info' carries everything the partial algorithms need; it is
        // shared with complementation (cf. helpers::tnba_complement and
        // kofola::Elevatorization, which builds it the same way).
        //
        // Two things to keep in mind:
        //  (1) cmpl_info stores most of its arguments *by reference*, so every
        //      structure passed below has to be a member of this class and
        //      outlive info_ (that is why part_to_scc_map_ is a member and not
        //      a local variable).
        //  (2) reachable_vector_, scc_to_pred_sccs_map_ and dir_sim_ are left
        //      empty on purpose - no determinization algorithm consumes them
        //      yet.  An algorithm that needs them (e.g. for simulation-based
        //      pruning of macrostates) has to fill them in here first.
        this->part_to_scc_map_ =
            helpers::tnba_complement::create_part_to_scc_map(std::get<3>(this->partitions_));
        this->info_ = std::make_unique<kofola::cmpl_info>(
                this->aut_,                        // automaton
                std::get<0>(this->partitions_),    // number of partitions
                std::get<1>(this->partitions_),    // partition types
                std::get<2>(this->partitions_),    // state to partition map
                this->reachable_vector_,           // vector of reachable states
                this->part_to_scc_map_,            // map of partitions to their SCCs
                this->scc_to_pred_sccs_map_,       // maps SCCs to their predecessors
                std::get<4>(this->partitions_),    // partitions to acceptance condition map
                *(this->scc_),                     // SCC information
                this->dir_sim_,                    // direct simulation
                this->is_accepting_,               // vector for acceptance of states
                false);                            // no shared breakpoint
    }

    // ----------------------------------------------------------------------
    // algorithm selection
    // ----------------------------------------------------------------------

    /// Picks one partial determinization algorithm per partition.  The partition
    /// types come from the very same analysis that drives complementation
    /// (helpers::tnba_complement::create_partitions), so both constructions
    /// agree on how the input automaton is decomposed.
    ///
    /// Together, the inherently weak and the deterministic algorithms cover
    /// every accepting SCC of an *elevator* automaton.  A nondeterministic
    /// accepting SCC is what makes an automaton non-elevator, and it is the only
    /// case left unimplemented - the factory fails loudly rather than silently
    /// producing a wrong automaton.
    void tela_determinize::select_algorithms() { // {{{
        using kofola::PartitionType;

        for (size_t i = 0; i < this->info_->num_partitions_; ++i) {
            alg_p alg;

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

    tela_determinize::alg_p
    tela_determinize::create_inherently_weak_algorithm(size_t partition_index) {
        return std::make_unique<kofola::determinize_mh>(*(this->info_.get()), partition_index);
    }

    tela_determinize::alg_p
    tela_determinize::create_deterministic_algorithm(size_t partition_index) {
        return std::make_unique<kofola::determinize_dac>(*(this->info_.get()), partition_index);
    }

    tela_determinize::alg_p
    tela_determinize::create_strongly_deterministic_algorithm(size_t partition_index) {
        return std::make_unique<kofola::determinize_dac>(*(this->info_.get()), partition_index);
    }

    tela_determinize::alg_p
    tela_determinize::create_nondeterministic_algorithm(size_t partition_index) {
        (void)partition_index;
        throw std::runtime_error("determinization of nondeterministic partitions is not implemented yet");
    }

    tela_determinize::alg_p
    tela_determinize::create_initial_almost_deterministic_algorithm(size_t partition_index) {
        // an initial almost deterministic SCC is in particular deterministic
        // inside (that is what SCC_DET_BORDER_NONDET_TYPE requires), so the DAC
        // algorithm applies to it as well
        return std::make_unique<kofola::determinize_dac>(*(this->info_.get()), partition_index);
    }

    // ----------------------------------------------------------------------
    // uberstate bookkeeping
    // ----------------------------------------------------------------------

    const tela_determinize::uberstate&
    tela_determinize::num_to_uberstate(unsigned num) const { // {{{
        assert(num < this->num_to_uberstate_map_.size());
        assert(this->num_to_uberstate_map_[num]);
        return *this->num_to_uberstate_map_[num];
    } // num_to_uberstate() }}}

    /// Inserts an uberstate and returns its number, or returns the number of an
    /// already present equal uberstate.  The uberstates are physically owned by
    /// num_to_uberstate_map_; uberstate_to_num_map_ is keyed by a *pointer* into
    /// that storage (ordered by uberstate_ptr_less_ftor), so the key inserted
    /// below must be the address of the stored copy and never of the caller's
    /// temporary.
    unsigned tela_determinize::insert_uberstate(const uberstate& us) { // {{{
        auto it = this->uberstate_to_num_map_.find(&us);
        if (this->uberstate_to_num_map_.end() == it) { // not found
            std::shared_ptr<uberstate> ptr(new uberstate(us));
            this->num_to_uberstate_map_.push_back(ptr);
            assert(this->num_to_uberstate_map_.size() == this->cnt_state_ + 1);  // invariant
            const uberstate* us_new = this->num_to_uberstate_map_[this->cnt_state_].get();
            auto jt_bool_pair = this->uberstate_to_num_map_.insert({us_new, this->cnt_state_});
            ++this->cnt_state_;
            return jt_bool_pair.first->second;
        } else { // found
            return it->second;
        }
    } // insert_uberstate() }}}

    // ----------------------------------------------------------------------
    // the construction
    // ----------------------------------------------------------------------

    /// Unlike in complementation, there is a single initial uberstate: each
    /// partial algorithm is deterministic and thus offers exactly one initial
    /// partial macrostate, so no Cartesian product (and no subsequent merging
    /// of several initial states into a fresh one) is needed.
    unsigned tela_determinize::get_initial_uberstate() { // {{{
        std::set<unsigned> initial_states = {aut_->get_init_state_number()};

        vec_macrostates vm;
        for (size_t i = 0; i < alg_vec_.size(); ++i) {
            vm.push_back(alg_vec_[i]->get_init());
        }

        return insert_uberstate(uberstate(initial_states, vm));
    } // get_initial_uberstate() }}}

    /// Computes the (unique) successor of an uberstate.  This is where the
    /// construction is glued together: 'all_succ' is a plain subset
    /// construction step over the *whole* input automaton, shared by all the
    /// partial algorithms.  Each of them extracts what it needs from it (the
    /// partitions cover accepting SCCs only, so states outside of every
    /// partition contribute to the reach set alone) - the reach set therefore
    /// also relieves the partial algorithms from tracking how runs enter their
    /// partition.
    ///
    /// Since every partial algorithm returns exactly one successor, so does
    /// this function; that is precisely what makes the result deterministic.
    tela_determinize::state_taggedcol tela_determinize::get_succ_uberstate(
            const uberstate& src,
            const bdd& symbol) { // {{{
        assert(alg_vec_.size() == src.get_part_macrostates().size());

        std::set<unsigned> all_succ = kofola::get_all_successors(
                this->aut_, src.get_reach_set(), symbol);

        DEBUG_PRINT_LN("all succ over " + std::to_string(symbol) + "= " + std::to_string(all_succ));

        const vec_macrostates& prev_part_macro = src.get_part_macrostates();
        vec_macrostates vm;
        std::set<std::pair<unsigned, unsigned>> cols;

        for (size_t i = 0; i < alg_vec_.size(); ++i) {
            abstract_determinize_alg::mstate_col mc =
                alg_vec_[i]->get_succ(all_succ, prev_part_macro[i].get(), symbol);

            // The colours reported by the partial algorithms are local: they all
            // start at 0 (or at get_min_colour()).  We tag them with the index of
            // the partition and only translate them into global colours in
            // run_new(), once set_acc_cond() has laid out the colour blocks.
            for (unsigned col : mc.second) {
                cols.insert({i, col});
            }
            vm.push_back(mc.first);
        }

        return {insert_uberstate(uberstate(all_succ, vm)), cols};
    } // get_succ_uberstate() }}}

    /// Lays out the colours of the result and assembles its acceptance
    /// condition.  Every accepting run of the input automaton is eventually
    /// confined to a single accepting SCC, hence L(A) is the *union* of the
    /// languages of the partitions and the global condition is the disjunction
    /// of the partial ones - the dual of the conjunction used by
    /// complementation (cf. helpers::tnba_complement::set_acc_cond()).
    ///
    /// Each partition gets a contiguous block of colours starting at
    /// part_col_offset_[i]; shifting its local condition by that offset is what
    /// keeps the blocks independent.  There is no reserved colour here: unlike
    /// the complement, the determinization has no accepting sink.
    void tela_determinize::set_acc_cond() { // {{{
        num_colours_ = 0;
        spot::acc_cond::acc_code alg_acc_code = spot::acc_cond::acc_code::f();  // neutral for |=

        for (size_t i = 0; i < this->info_->num_partitions_; ++i) {
            const spot::acc_cond& cond = alg_vec_[i]->get_acc_cond();
            vec_acc_code_.push_back(cond);
            spot::acc_cond::acc_code cond_code = cond.get_acceptance();

            cond_code <<= num_colours_;     // shift the partition's condition
            alg_acc_code |= cond_code;
            part_col_offset_[i] = num_colours_;
            num_colours_ += cond.num_sets();
        }

        if (num_colours_ > SPOT_MAX_ACCSETS) {
            throw std::runtime_error("the determinization needs " +
                std::to_string(num_colours_) + " colours, more than the " +
                std::to_string(SPOT_MAX_ACCSETS) + " supported by Spot");
        }

        final_code_ = alg_acc_code;
        DEBUG_PRINT_LN("colour offsets: " + std::to_string(part_col_offset_));
        DEBUG_PRINT_LN("final code: " + std::to_string(final_code_));
    } // set_acc_cond() }}}

    /// Used when the input has no accepting SCC at all: no run can then be
    /// accepting, so the language is empty and there is nothing to determinize.
    spot::twa_graph_ptr tela_determinize::make_empty_aut() const { // {{{
        spot::twa_graph_ptr result = spot::make_twa_graph(this->aut_->get_dict());
        result->copy_ap_of(this->aut_);
        result->set_acceptance(spot::acc_cond(0, spot::acc_cond::acc_code::f()));
        result->new_state();
        result->set_init_state(0);
        result->prop_universal(true);
        result->prop_complete(false);
        return result;
    } // make_empty_aut() }}}

    spot::twa_graph_ptr tela_determinize::run_new() { // {{{
        DEBUG_PRINT_LN("selecting algorithms");
        if (0 == std::get<0>(this->partitions_)) {
            // no accepting SCC at all - the language is empty
            return make_empty_aut();
        }

        // creates a vector of algorithms, one for every partition
        select_algorithms();
        DEBUG_PRINT_LN("algorithms selected");

        // The transitions of the result; also serves as the set of visited
        // states.  It has to be an *ordered* map: insert_uberstate() hands out
        // state numbers 0, 1, 2, ... in the order of first discovery, and the
        // conversion below re-creates them with new_state() while iterating over
        // this map, relying on the numbers coming out in ascending order.
        std::map<unsigned, std::vector<std::pair<bdd, state_taggedcol>>> det_states;

        std::stack<unsigned> todo;

        unsigned init_state = this->get_initial_uberstate();
        det_states.insert({init_state, {}});
        todo.push(init_state);

        while (!todo.empty()) { // the main loop
            unsigned us_num = todo.top();
            todo.pop();
            const uberstate& us = num_to_uberstate(us_num);

            auto it = det_states.find(us_num);
            assert(det_states.end() != it);
            std::vector<std::pair<bdd, state_taggedcol>>& us_post = it->second;

            DEBUG_PRINT_LN("processing " + std::to_string(us_num) + ": " + us.to_string());

            // compute support of all available states
            bdd msupport = bddtrue;
            bdd n_s_compat = bddfalse;
            const std::set<unsigned>& reach_set = us.get_reach_set();

            for (unsigned s : reach_set) {
                msupport &= support_[s];
                n_s_compat |= compat_[s];
            }

            // Symbols outside of 'n_s_compat' have no successor at all, so we
            // simply omit those edges and leave the result incomplete.
            // Complementation instead redirects them to an *accepting* sink;
            // here the dual would be a rejecting sink, which is never needed to
            // accept a word.  Note that an uberstate with an empty reach set has
            // n_s_compat == bddfalse and hence becomes a dead end, as it should.
            bdd all = n_s_compat;
            while (all != bddfalse) { // iterate over all symbols
                bdd letter = bdd_satoneset(all, msupport, bddfalse);
                all -= letter;

                DEBUG_PRINT_LN("symbol: " + std::to_string(letter));

                state_taggedcol succ = this->get_succ_uberstate(us, letter);
                us_post.emplace_back(std::make_pair(letter, succ));

                auto it_bool_pair = det_states.insert({succ.first, {}});
                if (it_bool_pair.second) { // the successor state is new
                    DEBUG_PRINT_LN("inserted " + std::to_string(succ.first) + " into det_states");
                    todo.push(succ.first);
                }
            }
        }

        set_acc_cond();
        spot::acc_cond result_cond(num_colours_, final_code_);

        // convert the result into a spot automaton
        spot::twa_graph_ptr result = spot::make_twa_graph(this->aut_->get_dict());
        result->copy_ap_of(this->aut_);
        // none of the properties of the input carries over; determinism is
        // claimed explicitly at the end of this function instead
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

        std::vector<std::string>* state_names = nullptr;
        if (show_names_) { // show names
            state_names = new std::vector<std::string>();
            result->set_named_prop("state-names", state_names);
        }

        for (const auto& st_trans_pair : det_states) {
            const unsigned& src = st_trans_pair.first;
            unsigned spot_state = result->new_state();
            (void) spot_state;
            assert(spot_state == src);

            for (const auto& bdd_tgt_pair : st_trans_pair.second) {
                const bdd& symbol = bdd_tgt_pair.first;
                const unsigned& tgt = bdd_tgt_pair.second.first;
                const std::set<std::pair<unsigned, unsigned>>& cols = bdd_tgt_pair.second.second;

                // translate the partition-local colours into the global ones;
                // get_min_colour() rebases algorithms whose colours do not start
                // at 0 (they may only discover their range during the run, so
                // this must not happen before the main loop has finished)
                std::vector<unsigned> new_cols;
                for (const std::pair<unsigned, unsigned>& part_col_pair : cols) {
                    const unsigned part_index = part_col_pair.first;
                    const unsigned colour = part_col_pair.second;

                    unsigned shift = alg_vec_[part_index]->get_min_colour();
                    if (colour < shift ||
                        colour - shift >= vec_acc_code_[part_index].num_sets()) {
                        // the colour is outside of the range the algorithm
                        // finally declared, hence it occurs in no disjunct of
                        // its acceptance condition - see the contract of
                        // abstract_determinize_alg::get_acc_cond()
                        continue;
                    }
                    new_cols.push_back(part_col_offset_.at(part_index) + colour - shift);
                }
                spot::acc_cond::mark_t spot_cols(new_cols.begin(), new_cols.end());
                result->new_edge(src, tgt, symbol, spot_cols);
            }

            if (this->show_names_) { // handle output state names
                state_names->push_back(num_to_uberstate(src).to_string());
            }
        }

        result->set_init_state(init_state);
        result->prop_universal(true);

        if (kofola::LOG_VERBOSITY > 0) {
            spot::print_hoa(std::cerr, result);
            std::cerr << "\n\n\n\n";
        }

        assert(spot::is_deterministic(result));

        return result;
    } // run_new() }}}

    // ######################################################################

    /// Does the automaton have an accepting SCC that is neither inherently weak
    /// nor deterministic inside?  Those are exactly the partitions for which no
    /// partial determinization algorithm exists.  Only accepting SCCs matter:
    /// the partitions cover no other ones.
    static bool has_nondet_acc_scc(const spot::twa_graph_ptr& aut)
    { // {{{
        spot::scc_info si(aut, spot::scc_info_options::ALL);
        si.determine_unknown_acceptance();
        std::string scc_types = helpers::get_scc_types(si);

        for (unsigned scc = 0; scc < si.scc_count(); ++scc) {
            if (helpers::is_accepting_nondetscc(scc_types, scc)) { return true; }
        }

        return false;
    } // has_nondet_acc_scc() }}}

    /// Spot-based reduction of the determinized automaton 'det'.  Unlike
    /// apply_postprocessing() used by complementation, the requested output type
    /// is respected (Generic, i.e. Emerson-Lei, by default) and the result is
    /// never forced into a TGBA: that conversion only succeeds for DBA-
    /// recognizable languages and otherwise wastes a whole run of Spot on a
    /// nondeterministic automaton that is thrown away.  'orig' is the input of
    /// the determinization, used by the heuristic deciding whether the
    /// reduction is affordable.
    static spot::twa_graph_ptr postprocess_det(const spot::twa_graph_ptr& det,
                                               const spot::twa_graph_ptr& orig)
    { // {{{
        if (kofola::has_value("raw", "yes", kofola::OPTIONS.params) ||
            !kofola::is_post_reduction_suitable(orig)) {
            return det;
        }

        spot::postprocessor p;
        p.set_pref(spot::postprocessor::Deterministic);
        p.set_level(spot::postprocessor::Low);
        if ("buchi" == kofola::OPTIONS.output_type) {
            p.set_type(spot::postprocessor::Buchi);
        } else if ("tgba" == kofola::OPTIONS.output_type) {
            p.set_type(spot::postprocessor::GeneralizedBuchi);
        } else {
            p.set_type(spot::postprocessor::Generic);
        }

        return p.run(det);
    } // postprocess_det() }}}

    spot::twa_graph_ptr determinize_tela(const spot::twa_graph_ptr& aut) {
        if (spot::is_deterministic(aut)) { // nothing to do
            return aut;
        }

        // There is no partial determinization algorithm for a nondeterministic
        // accepting SCC, so a non-elevator automaton is first limit-determinized
        // into an elevator one.  complement_tela() preprocesses the same way,
        // but with only_non_buchi = true: complementation has NCSB for a Buchi
        // NAC and only needs the other ones removed.  Here every NAC has to go,
        // hence the default only_non_buchi = false.
        //
        // Elevatorization rewrites its automaton in place, so it is handed a
        // copy - 'aut' is still needed as the reference for the postprocessing
        // heuristic at the end.
        spot::twa_graph_ptr input = aut;
        if (has_nondet_acc_scc(aut)) {
            spot::twa_graph_ptr copy =
                spot::make_twa_graph(aut, spot::twa::prop_set::all());
            kofola::Elevatorization elev(copy);
            input = elev.elevatorize();

            // limit-determinization may already have made it deterministic
            if (spot::is_deterministic(input)) { return input; }
        }

        auto scc = std::make_shared<spot::scc_info>(input, spot::scc_info_options::ALL);

        // Spot's is_accepting_scc might say "unknown" for Fin conditions; without
        // this, create_partitions() would see no accepting SCC at all
        scc->determine_unknown_acceptance();

        // SCC partitioning (via the same helper used by complementation).
        // Return tuple items:
        //  (0) `num_partitions`     - number of created partitions
        //  (1) `part_to_type_map`   - partition index -> PartitionType
        //  (2) `st_to_part_map`     - state index -> partition index
        //  (3) `scc_to_part_map`    - SCC index -> partition index
        //  (4) `part_to_acc_map`    - partition index -> restricted acceptance condition
        // Unlike complementation, determinization does not benefit from merging
        // all the DACs into a single partition: the DAC algorithm needs one
        // block of colours per potentially live run, i.e. per state of the
        // partition, so a merged partition burns the (very limited) colour
        // budget of Spot for nothing.  Keep every DAC separate unless the user
        // asked otherwise.
        //
        // When kofola is run from the command line, main() has already set this
        // default (it is operation-aware), so the branch below only fires for
        // callers that use this function as a library and never populate the
        // parameter map - the unit tests, most notably.
        kofola::options det_options = kofola::OPTIONS;
        if (0 == det_options.params.count("merge_det")) {
            det_options.params["merge_det"] = "no";
        }

        auto partitions = helpers::tnba_complement::create_partitions(*scc, det_options);

        tela_determinize det(input, scc, std::move(partitions));
        auto res = det.run_new();

        // Postprocessing is an optimization, never a semantic step, so anything
        // that goes wrong in it leaves the result of the construction untouched.
        // Two things do go wrong in practice:
        //  - Spot's Deterministic preference is only a preference, so the
        //    reduction may hand back a nondeterministic automaton;
        //  - the reduction may need more than the SPOT_MAX_ACCSETS colours Spot
        //    can represent, even though the automaton we produced fits.
        try {
            auto post = postprocess_det(res, aut);
            if (spot::is_deterministic(post)) { return post; }
        }
        catch (const std::runtime_error&) { }

        return res;
    }
}

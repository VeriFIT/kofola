#include "inclusion_test.hpp"
#include <tuple>

// kofola
#include "complement_tela.hpp"
#include "util.hpp"
#include "decomposer.hpp"
#include "complement_sync.hpp"

// Spot
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complete.hh>

// standard library
#include <queue>

namespace kofola {
    inclusionTest::inclusionTest() {
        ;
    }

    spot::twa_graph_ptr inclusionTest::preprocess(const spot::twa_graph_ptr &aut) {
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
                auto aut_cmp = [](const auto &lhs, const auto &rhs) { return lhs->num_states() > rhs->num_states(); };
                std::priority_queue<spot::twa_graph_ptr,
                        std::vector<spot::twa_graph_ptr>, decltype(aut_cmp)> aut_queue(aut_cmp);
                for (auto aut: decomposed) {
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

        //spot::print_hoa(std::cerr, aut_to_compl);
        spot::complete_here(aut_to_compl);
        //spot::print_hoa(std::cerr, aut_to_compl);

        return aut_to_compl;
    }

    std::pair<bdd, bdd> inclusionTest::symbols_from_A(const spot::twa_graph_ptr &aut_A) {
        bdd res = bddfalse;
        bdd msupport = bddtrue;

        for (unsigned s = 0; s < aut_A->num_states(); s++) {
            for (const auto &t: aut_A->out(s)) {
                res |= t.cond;
                msupport &= bdd_support(t.cond);
                DEBUG_PRINT_LN("automaton alphabet symbol: " + std::to_string(t.cond));
            }
        }

        return std::make_pair(res, msupport);
    }

    bool inclusionTest::is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol) {
        for (const auto &t: aut_A->out(src)) {
            if (dst == t.dst && bdd_implies(symbol, t.cond)) {
                if (t.acc) { return true; }
                else return false;
            }
        }

        return false;
    }

    std::set<inclusionTest::state_col>
    inclusionTest::get_cross_prod(const spot::twa_graph_ptr &aut_A, unsigned aut_A_src, cola::tnba_complement &aut_B,
                                   std::set<unsigned> &states_A, cola::tnba_complement::vec_state_taggedcol &states_B,
                                   std::map<intersect_mstate, vec_state_col> &intersect_states, const bdd &letter) {
        std::set<inclusionTest::state_col> cross_prod;

        auto vec_acc_cond = aut_B.get_vec_acc_cond();
        auto part_col_offset = aut_B.get_part_col_offset();

        for (const auto &state_A: states_A) {
            for (const auto &state_cols: states_B) {
                const unsigned &state_B = state_cols.first;
                auto cols = state_cols.second;
                std::set<unsigned> new_cols;

                if (is_transition_acc(aut_A, aut_A_src, state_A, letter)) {
                    new_cols.insert(UINT_MAX-1); // TODO what colour here????
                }

                for (const std::pair<unsigned, unsigned> &part_col_pair: cols) {
                    const unsigned part_index = part_col_pair.first;
                    const unsigned colour = part_col_pair.second;

                    //else if (UINT_MAX == colour) { // FIXME: this is a special way
                    //    // of dealing with determinization-based (this should set the
                    //    // colour to the maximum colour) - try to make more uniform

                    //    unsigned max_col = vec_acc_cond[part_index].num_sets() - 1;
                    //    new_cols.insert(part_col_offset.at(part_index) + max_col);
                    //}
                    //else { // standard transition
                        unsigned shift = aut_B.get_alg_vec_mincolour_at_i(
                                part_index); // how much to decrement the colour
                        new_cols.insert(part_col_offset.at(part_index) + colour - shift);
                        inf_acc_cols_conj_.insert(part_col_offset.at(part_index) + colour - shift);
                    //}
                }

                auto succ = std::make_pair(std::make_pair(state_A, state_B), new_cols);
                cross_prod.insert(succ);

                //auto it_bool_pair =
                intersect_states.insert({{state_A, state_B},{}});
                /*if (it_bool_pair.second) { // the successor state is new
                    todo_.push({state_A, succ_state_B});
                }*/
            }
        }

        return cross_prod;
    }

    std::set<inclusionTest::state_col>
    inclusionTest::get_all_succs(const spot::twa_graph_ptr &aut_A, cola::tnba_complement &aut_B,
                                  intersect_mstate &mstate,
                                  std::map<intersect_mstate, vec_state_col> &intersect_states) {
        std::set<inclusionTest::state_col> cross_prod;
        std::set<inclusionTest::state_col> result;

        // prepare to insert post of current intersection macrostate
        auto it = intersect_states.find(mstate);
        assert(intersect_states.end() != it);
        vec_state_col &intersect_post = it->second;

        // extract A state and compl.B state from intersection macrostate to compute successors
        unsigned state_of_A = mstate.first;

        auto tmp_bdds = symbols_from_A(aut_A);
        bdd msupport = tmp_bdds.second;
        bdd n_s_compat = tmp_bdds.first;
        bdd all = n_s_compat;

        // iterate over all symbols
        while (all != bddfalse) {
            bdd letter = bdd_satoneset(all, msupport, bddfalse);
            all -= letter;

            DEBUG_PRINT_LN("symbol: " + std::to_string(letter));

            std::set<unsigned> succs_A = get_all_successors(aut_A, std::set<unsigned>{state_of_A}, letter);
            cola::tnba_complement::vec_state_taggedcol succs_B;

            if(!aut_B.get_is_sink_created() || mstate.second != aut_B.get_sink_state())
            {

                const cola::tnba_complement::uberstate &us_B = aut_B.num_to_uberstate(mstate.second);
                succs_B = aut_B.get_succ_uberstates(us_B, letter);
                //for(auto state_orig_B: us_B.get_reach_set())
                //{
                //    auto part_succ = get_all_successors(preprocessed_orig_aut_B_, std::set<unsigned>{state_orig_B}, letter);
                //    if(part_succ.empty())
                //    {
                //        aut_B.handle_sink_state();
                //        auto sink_state = aut_B.get_sink_state();
                //        succs_B.push_back({sink_state, {}});
                //    }
//
                //}
            }
            else
            {
                succs_B.push_back({aut_B.get_sink_state(), {}});
            }

            if(!succs_A.empty())
            {
                cross_prod = get_cross_prod(aut_A, state_of_A, aut_B, succs_A, succs_B, intersect_states,
                                            letter);
                result.insert(cross_prod.begin(), cross_prod.end());
                intersect_post.insert(intersect_post.end(), cross_prod.begin(), cross_prod.end());
            }
        }

        return result;
    }

    void inclusionTest::tarjan_is_empty(inclusionTest::intersect_mstate &src_mstate, const spot::twa_graph_ptr &aut_A,
                                         cola::tnba_complement &aut_B, unsigned ith_on_path, std::vector<std::set<unsigned>> cols_path) {
        /// STRONGCONNECT
        cols_visited_.insert({src_mstate, std::set<unsigned>()}); // init

        // mark position in current dfs
        if(dfs_state_pos_.count(src_mstate) == 0)
            dfs_state_pos_.insert({src_mstate, ith_on_path});
        else
            dfs_state_pos_.at(src_mstate) = ith_on_path;

        indices_.at(src_mstate) = index_;
        lowlinks_.at(src_mstate) = index_;
        index_++;
        tarjan_stack_.push(src_mstate);
        on_stack_.at(src_mstate) = true;

        auto succs = get_all_succs(aut_A, aut_B, src_mstate, intersect_states_);
        if(succs.empty())
            return;

        for (const auto &state_cols_pair: succs) {
            auto dst_mstate = state_cols_pair.first;
            auto cols = state_cols_pair.second;
            cols_path.emplace_back(cols);

            if(indices_.count(dst_mstate) == 0)
            {
                indices_.insert({dst_mstate, UNDEFINED});
                lowlinks_.insert({dst_mstate, UNDEFINED});
                on_stack_.insert({dst_mstate, false});
            }

            if (indices_.at(dst_mstate) == UNDEFINED)
            {
                tarjan_is_empty(dst_mstate, aut_A, aut_B, ++ith_on_path, cols_path);
                // propagate visited cols when backtracking
                if(lowlinks_.at(src_mstate) < lowlinks_.at(dst_mstate))
                {
                    cols_visited_.at(src_mstate).insert(cols_visited_.at(dst_mstate).begin(), cols_visited_.at(dst_mstate).end());
                }

                // acc condition check
                bool is_subset = std::includes(cols_visited_.at(src_mstate).begin(), cols_visited_.at(src_mstate).end(), inf_acc_cols_conj_.begin(), inf_acc_cols_conj_.end());;
                if(is_subset)
                {
                    decided_ = true;
                    empty_ = false;

                    return;
                }

                cols_path.at(cols_path.size() - 1).insert(cols_visited_.at(src_mstate).begin(), cols_visited_.at(src_mstate).end()); // add possible lasso colors to the current node on DFS traversal

                lowlinks_.at(src_mstate) = std::min(lowlinks_.at(src_mstate), lowlinks_.at(dst_mstate));

                if(decided_)
                    return;

            } else if (on_stack_.at(dst_mstate)) { // visiting some lasso
                lowlinks_.at(src_mstate) = std::min(lowlinks_.at(src_mstate), indices_.at(dst_mstate));

                if(cols_visited_.count(dst_mstate) == 0) // if not initialized
                    cols_visited_.insert({dst_mstate, std::set<unsigned>()});

                // collecting all colors seen in the lasso
                for(unsigned i = dfs_state_pos_.at(dst_mstate); i < cols_path.size(); i++)
                {
                    // union
                    cols_visited_.at(dst_mstate).insert(cols_path.at(i).begin(), cols_path.at(i).end());
                }

                // propagate visited cols
                cols_visited_.at(src_mstate).insert(cols_visited_.at(dst_mstate).begin(), cols_visited_.at(dst_mstate).end());
                // acc cond check
                bool is_subset = std::includes(cols_visited_.at(dst_mstate).begin(), cols_visited_.at(dst_mstate).end(), inf_acc_cols_conj_.begin(), inf_acc_cols_conj_.end());;
                if(is_subset)
                {
                    decided_ = true;
                    empty_ = false;

                    return;
                }
            }
        }

        if (lowlinks_.at(src_mstate) == indices_.at(src_mstate)) {
            intersect_mstate tmp_mstate; // for poping out of the stack
            std::set<unsigned> tmp_set; // to track colors in scc

            do {
                tmp_mstate = tarjan_stack_.top();
                tarjan_stack_.pop();

                // union of all colours seen in scc
                if(cols_visited_.count(tmp_mstate) != 0)
                    tmp_set.insert(cols_visited_.at(tmp_mstate).begin(), cols_visited_.at(tmp_mstate).end());

                on_stack_.at(tmp_mstate) = false;
            } while (src_mstate != tmp_mstate);

            //scc_cols_.emplace_back(tmp_set);

            bool is_subset = std::includes(tmp_set.begin(), tmp_set.end(), inf_acc_cols_conj_.begin(), inf_acc_cols_conj_.end());
            if(is_subset)
            {
                decided_ = true;
                empty_ = false;

                return;
            }
        }
    }

    bool inclusionTest::test(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B) {
        auto si_A = spot::scc_info(aut_A, spot::scc_info_options::ALL);
        auto preprocessed_aut_A = kofola::saturate(aut_A, si_A);
        preprocessed_aut_A->prop_state_acc(false);

        unsigned init_A = preprocessed_aut_A->get_init_state_number();

        preprocessed_orig_aut_B_ = preprocess(aut_B);

        spot::scc_info si_B(preprocessed_orig_aut_B_, spot::scc_info_options::ALL);
        auto comp_aut_B = cola::tnba_complement(preprocessed_orig_aut_B_, si_B);

        DEBUG_PRINT_LN("selecting algorithms");
        // creates a vector of algorithms, for every SCC of aut one
        comp_aut_B.select_algorithms();
        DEBUG_PRINT_LN("algorithms selected");

        // get initial uberstates
        auto init_vec{comp_aut_B.get_initial_uberstates()};
        std::stack<intersect_mstate> init_states;
        for (unsigned state: init_vec) {
            intersect_states_.insert({{init_A, state},
                                      {}});
            indices_.insert({{init_A, state}, UNDEFINED});
            lowlinks_.insert({{init_A, state}, UNDEFINED});
            on_stack_.insert({{init_A, state}, false});
            init_states.push({init_A, state});
        }

        std::vector<std::set<unsigned>> cols_path; // cache for colors seen during the dfs traversal
        inf_acc_cols_conj_ = comp_aut_B.set_acc_cond(); // TODO take acc cond from it
        inf_acc_cols_conj_.insert(UINT_MAX-1); // for automaton A
        //acceptance_cond_ = comp_aut_B.get_final_acc_code();

        /*******************************************************
         *                                                     *
         *                        Tarjan                       *
         *                                                     *
         *******************************************************/
        while (!init_states.empty()) { // search from all initial states
            intersect_mstate src_mstate = init_states.top();
            init_states.pop();
            if (indices_.at(src_mstate) == UNDEFINED) {
                // init
                if(dfs_state_pos_.count(src_mstate) == 0)
                    dfs_state_pos_.insert({src_mstate, 0});
                else
                    dfs_state_pos_.at(src_mstate) = 0;
                tarjan_is_empty(src_mstate, preprocessed_aut_A, comp_aut_B, 0, cols_path);
                if(decided_)
                    return empty_;
            }
        }

        return true;
    }
}// namespace KOFOLA
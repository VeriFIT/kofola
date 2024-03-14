// kofola
#include "inclusion_check.hpp"
#include "kofola.hpp"
#include "complement_tela.hpp"
#include "util.hpp"
#include "decomposer.hpp"
#include "complement_sync.hpp"

// spot
#include <spot/twaalgos/postproc.hh>
#include <spot/twaalgos/product.hh>
#include <spot/twaalgos/complete.hh>

namespace kofola {
    inclusion_check::inclusion_check(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B)
    : aut_B_compl_(init_compl_aut_b(aut_B)){
        auto si_A = spot::scc_info(aut_A, spot::scc_info_options::ALL);
        auto preprocessed_aut_A = kofola::saturate(aut_A, si_A);
        preprocessed_aut_A->prop_state_acc(false);
        aut_A_ = preprocessed_aut_A;

        unsigned init_A = preprocessed_aut_A->get_init_state_number();

        DEBUG_PRINT_LN("selecting algorithms");
        // creates a vector of algorithms, for every SCC of aut one
        aut_B_compl_.select_algorithms();
        DEBUG_PRINT_LN("algorithms selected");

        // get initial uberstates
        auto init_vec{aut_B_compl_.get_initial_uberstates()};
        for (unsigned state: init_vec) {
            auto ptr = std::make_shared<inclusion_mstate>();
            ptr->state_ = {init_A, state};

            init_states_.emplace_back(std::move(ptr));
            intersect_states_.insert({{init_A, state},
                                      {}});
        }

        // setting accepting cond to acc of complement & Inf(color of aut_A=UINTMAX-1)
        first_col_to_use_ = static_cast<unsigned int>(aut_B_compl_.set_acc_cond());
        acc_cond_ = aut_B_compl_.get_final_acc_code();
        acc_cond_ &= spot::acc_cond::acc_code::inf({first_col_to_use_});
    }

    cola::tnba_complement inclusion_check::init_compl_aut_b(const spot::twa_graph_ptr &aut_B) {
        preprocessed_orig_aut_B_ = preprocess(aut_B);

        spot::scc_info si_B(preprocessed_orig_aut_B_, spot::scc_info_options::ALL);
        cola::tnba_complement comp(preprocessed_orig_aut_B_, si_B);
        return comp;
    }

    bool inclusion_check::inclusion() {
        emptiness_check emptiness_checker(this, INCLUSION);
        auto res = emptiness_checker.empty();
        return res;
    }

    bool inclusion_check::is_accepting(spot::acc_cond::mark_t inf_cond) {
        return acc_cond_.accepting(inf_cond);
    }

    std::vector<std::shared_ptr<abstract_successor::mstate>> inclusion_check::get_initial_states() {
        std::vector<std::shared_ptr<abstract_successor::mstate>> base_mstates;

        for (auto& init : init_states_) {
            base_mstates.emplace_back(std::move(init));
        }

        return base_mstates;
    }

    spot::twa_graph_ptr inclusion_check::preprocess(const spot::twa_graph_ptr &aut) {
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

    std::pair<bdd, bdd> inclusion_check::symbols_from_A(const spot::twa_graph_ptr &aut_A) {
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

    std::vector<std::shared_ptr<abstract_successor::mstate>> inclusion_check::get_succs(const std::shared_ptr<abstract_successor::mstate> &src) {
        std::vector<std::shared_ptr<inclusion_mstate>> cartesian_prod;
        std::vector<std::shared_ptr<abstract_successor::mstate>> result;

        auto casted_src = dynamic_cast<inclusion_mstate*>(src.get());

        // extract A state and compl.B state from intersection macrostate to compute successors
        unsigned state_of_A = casted_src->state_.first;

        auto tmp_bdds = symbols_from_A(aut_A_);
        bdd msupport = tmp_bdds.second;
        bdd n_s_compat = tmp_bdds.first;
        bdd all = n_s_compat;

        // iterate over all symbols
        while (all != bddfalse) {
            bdd letter = bdd_satoneset(all, msupport, bddfalse);
            all -= letter;

            DEBUG_PRINT_LN("symbol: " + std::to_string(letter));

            std::set<unsigned> succs_A = get_all_successors(aut_A_, std::set<unsigned>{state_of_A}, letter);
            cola::tnba_complement::vec_state_taggedcol succs_B;

            if(!aut_B_compl_.get_is_sink_created() || casted_src->state_.second != aut_B_compl_.get_sink_state())
            {

                const cola::tnba_complement::uberstate &us_B = aut_B_compl_.num_to_uberstate(casted_src->state_.second);
                succs_B = aut_B_compl_.get_succ_uberstates(us_B, letter);
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
                succs_B.push_back({aut_B_compl_.get_sink_state(), {}});
            }

            if(!succs_A.empty())
            {
                cartesian_prod = get_cartesian_prod(state_of_A, succs_A, succs_B,
                                                    letter);
                for (auto& succ : cartesian_prod) {
                    // Downcast the shared pointer to hyperltl_mc_mstate
                    result.emplace_back(std::move(succ));
                }
            }
        }

        return result;
    }

    bool inclusion_check::is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol) {
        for (const auto &t: aut_A->out(src)) {
            if (dst == t.dst && bdd_implies(symbol, t.cond)) {
                if (t.acc) { return true; } // TODO is it valid???
                else return false;
            }
        }

        return false;
    }

    std::vector<std::shared_ptr<inclusion_mstate>>
    inclusion_check::get_cartesian_prod(unsigned aut_A_src, std::set<unsigned> &states_A,
                                        cola::tnba_complement::vec_state_taggedcol &states_B, const bdd &letter) {
        std::vector<std::shared_ptr<inclusion_mstate>> cartesian_prod;

        auto vec_acc_cond = aut_B_compl_.get_vec_acc_cond();
        auto part_col_offset = aut_B_compl_.get_part_col_offset();

        for (const auto &state_A: states_A) {
            for (const auto &state_cols: states_B) {
                const unsigned &state_B = state_cols.first;
                auto cols = state_cols.second;
                std::set<unsigned> new_cols;

                if (is_transition_acc(aut_A_, aut_A_src, state_A, letter)) {
                    new_cols.insert(first_col_to_use_); // TODO what colour here????
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
                    unsigned shift = aut_B_compl_.get_alg_vec_mincolour_at_i(
                            part_index); // how much to decrement the colour
                    new_cols.insert(part_col_offset.at(part_index) + colour - shift);
                    //}
                }

                //auto succ = std::make_pair(std::make_pair(state_A, state_B), new_cols);
                auto succ = std::make_shared<inclusion_mstate>();
                succ->state_ = std::make_pair(state_A, state_B);
                spot::acc_cond::mark_t spot_cols(new_cols.begin(), new_cols.end());
                succ->acc_ = spot_cols;
                cartesian_prod.emplace_back(std::move(succ));

                //auto it_bool_pair =
                intersect_states_.insert({{state_A, state_B},{}});
                /*if (it_bool_pair.second) { // the successor state is new
                    todo_.push({state_A, succ_state_B});
                }*/
            }
        }

        return cartesian_prod;
    }
}
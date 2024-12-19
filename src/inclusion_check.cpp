/**
 * @file inclusion_check.cpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief implementation of on the fly inclusion checking
 * @version 0.1
 * @date 2024-05-03
 * 
 * 
 */

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
    bool operator<(const inclusion_mstate& lhs,
                   const inclusion_mstate& rhs) {
            return lhs.lt(rhs);
    }

    inclusion_check::inclusion_check(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B)
    : aut_A_(init_aut_A(aut_A)), aut_B_compl_(init_compl_aut_b(aut_B)){
        auto tmp_bdds = symbols_from_A(aut_A_);
        msupport_ = tmp_bdds.second;
        n_s_compat_ = tmp_bdds.first;

        unsigned init_A = aut_A_->get_init_state_number();

        DEBUG_PRINT_LN("selecting algorithms");
        // creates a vector of algorithms, for every SCC of aut one
        aut_B_compl_.select_algorithms();
        DEBUG_PRINT_LN("algorithms selected");

        if(kofola::OPTIONS.params.count("dir_sim") != 0 && kofola::OPTIONS.params["dir_sim"] == "yes")
            compute_simulation(aut_A_, aut_B);

        // store initial uberstates
        auto init_vec{aut_B_compl_.get_initial_uberstates()};
        for (unsigned state: init_vec) {
            auto uberstate = aut_B_compl_.num_to_uberstate(state);
            auto B_reach_set = uberstate.get_reach_set();
            std::set<unsigned> A_B_intersect;

            if(dir_simul_.count(init_A)) {
                set_intersection(dir_simul_[init_A].begin(), dir_simul_[init_A].end(), B_reach_set.begin(), B_reach_set.end(),
                            std::inserter(A_B_intersect, A_B_intersect.begin()));
            }
            if(A_B_intersect.empty()) {
                auto ptr = std::make_shared<inclusion_mstate>();
                ptr->state_ = {init_A, state};

                init_states_.emplace_back(std::move(ptr));
                intersect_states_.insert({{init_A, state},
                                        {}});
            }
        }

        // setting accepting cond to acc of complement & Inf(the first unused color by aut_B_compl)
        infs_from_compl_ = (aut_B_compl_.set_acc_cond());
        first_col_to_use_ = infs_from_compl_.size() + 1;
        acc_cond_ = aut_B_compl_.get_final_acc_code();
        acc_cond_ &= spot::acc_cond::acc_code::inf({first_col_to_use_});
    }

    spot::twa_graph_ptr inclusion_check::init_aut_A(const spot::twa_graph_ptr &aut_A) {
        spot::twa_graph_ptr res = aut_A;
        if(kofola::OPTIONS.params.count("preproc_incl_A") != 0 && kofola::OPTIONS.params["preproc_incl_A"] == "high") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::High);
            res = p.run(aut_A);
        }
        else if(kofola::OPTIONS.params.count("preproc_incl_A") != 0 && kofola::OPTIONS.params["preproc_incl_A"] == "medium") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::Medium);
            res = p.run(aut_A);
        }
        else if(kofola::OPTIONS.params.count("preproc_incl_A") != 0 && kofola::OPTIONS.params["preproc_incl_A"] == "low") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::Low);
            res = p.run(aut_A);
        }

        return res;
    }

    cola::tnba_complement inclusion_check::init_compl_aut_b(const spot::twa_graph_ptr &aut_B) {
        spot::twa_graph_ptr aut_to_compl = aut_B;
        if(kofola::OPTIONS.params.count("preproc_incl_B") != 0 && kofola::OPTIONS.params["preproc_incl_B"] == "high") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::High);
            aut_to_compl = p.run(aut_B);
        }
        else if(kofola::OPTIONS.params.count("preproc_incl_B") != 0 && kofola::OPTIONS.params["preproc_incl_B"] == "medium") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::Medium);
            aut_to_compl = p.run(aut_B);
        }
        else if(kofola::OPTIONS.params.count("preproc_incl_B") != 0 && kofola::OPTIONS.params["preproc_incl_B"] == "low") {
            spot::postprocessor p;
            p.set_type(spot::postprocessor::Buchi);
            p.set_level(spot::postprocessor::Low);
            aut_to_compl = p.run(aut_B);
        }


        kofola::OPTIONS.output_type = "tgba";
        spot::scc_info si_B(aut_B, spot::scc_info_options::ALL);
        cola::tnba_complement comp(aut_to_compl, si_B);
        return comp;
    }

    spot::twa_graph_ptr inclusion_check::aut_union(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B) {
        auto res = spot::make_twa_graph(aut_A->get_dict());
        res->copy_ap_of(aut_A);
        res->set_acceptance(aut_A->acc());
        offset_ = aut_A->num_states();

        unsigned init_a;
        unsigned init_b;
        unsigned new_st;

        res->copy_state_names_from(aut_A);
        for(unsigned i = 0; i < aut_A->num_states(); i++) {
            new_st = res->new_state();
            if(i == aut_A->get_init_state_number()) {
                init_a = new_st;
            }
            for(const auto& t: aut_A->out(i)) {
                res->new_edge(new_st, t.dst, t.cond, t.acc);
            }
        }

        for(unsigned i = 0; i < aut_B->num_states(); i++) {
            new_st = res->new_state();
            if(i == aut_B->get_init_state_number()) {
                init_b = new_st;
            }
            for(const auto& t: aut_B->out(i)) {
                res->new_edge(new_st, offset_ + t.dst, t.cond, t.acc);
            }
        }

        // creating initial state such that it has transitions to both initial states of A and B respectively
        new_st = res->new_state();
        res->new_edge(new_st, init_a, bdd_true());
        res->new_edge(new_st, init_b, bdd_true());
        res->set_init_state(new_st);
        
        return res;
    }

    void inclusion_check::compute_simulation(const spot::twa_graph_ptr &aut_A, const spot::twa_graph_ptr &aut_B) {
        auto uni = aut_union(aut_B, aut_A);
        //spot::print_hoa(std::cout, uni);

        auto reduced = spot::simulation(uni);
        auto x = reduced->get_named_prop<std::vector<unsigned>>("simulated-states"); // to know which states were merged where
        auto orig_to_new = *x;

        // [...,last_of_B,...,last_of_A]
        for(unsigned i = offset_; i < orig_to_new.size() - 1; i++) {
            for(unsigned j = 0; j < offset_; j++) {
                //if state of A merged to state of B => simulated, not considering transition prunning (-1 in condition)
                if(orig_to_new[i] == orig_to_new[j] && orig_to_new[i] < offset_ && orig_to_new[i] != -1) {
                    dir_simul_[i - offset_].emplace_back(j);
                }
            }
        }
    }

    bool inclusion_check::inclusion() {
        emptiness_check emptiness_checker(this);
        auto res = emptiness_checker.empty();
        return res;
    }

    bool inclusion_check::is_accepting(spot::acc_cond::mark_t inf_cond) {
        return acc_cond_.accepting(inf_cond);
    }

    std::vector<std::shared_ptr<inclusion_mstate>> inclusion_check::get_initial_states() {
        std::vector<std::shared_ptr<inclusion_mstate>> base_mstates;

        for (auto& init : init_states_) {
            base_mstates.emplace_back(std::move(init));
        }

        return base_mstates;
    }

    spot::twa_graph_ptr inclusion_check::preprocess(const spot::twa_graph_ptr &aut) {
        spot::twa_graph_ptr aut_reduced = aut;
        std::vector<bdd> implications;

        // make sure the input is a BA
        spot::postprocessor p;
        p.set_type(spot::postprocessor::Buchi);
        p.set_level(spot::postprocessor::High);
        spot::twa_graph_ptr aut_to_compl;
        aut_to_compl = p.run(aut_reduced);

        //spot::print_hoa(std::cerr, aut_to_compl);

        /// complete to avoid sink state in complement
        // spot::complete_here(aut_to_compl);
        spot::complete_here(aut_reduced);

        //spot::print_hoa(std::cerr, aut_to_compl);

        // return aut_to_compl;
        return aut_reduced;
    }

    std::pair<bdd, bdd> inclusion_check::symbols_from_A(const spot::twa_graph_ptr &aut_A) {
        // taken from complement_sync
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

    void inclusion_check::print_mstate(const std::shared_ptr<inclusion_mstate> a) {
        auto casted_a = dynamic_cast<inclusion_mstate*>(a.get());

        std::cout << casted_a->acc_ << "\n";
        std::cout << casted_a->state_.first << ", " << std::to_string(aut_B_compl_.num_to_uberstate(casted_a->state_.second));
        std::cout << std::endl << "==================" << std::endl; 
    }

    bool inclusion_check::subsum_less_early(const std::shared_ptr<inclusion_mstate> a, const std::shared_ptr<inclusion_mstate> b) {
        auto casted_a = a;
        auto casted_b = b;

        // between states of A only identity
        return (casted_a->state_.first == casted_b->state_.first && aut_B_compl_.subsum_less_early(casted_a->state_.second, casted_b->state_.second));
    }

    bool inclusion_check::subsum_less_early_plus(const std::shared_ptr<inclusion_mstate> a, const std::shared_ptr<inclusion_mstate> b) {
        auto casted_a = a;
        auto casted_b = b;

        // between states of A only identity
        return (casted_a->state_.first == casted_b->state_.first && aut_B_compl_.subsum_less_early_plus(casted_a->state_.second, casted_b->state_.second));
    }

    cola::tnba_complement::vec_state_taggedcol inclusion_check::get_successors_compl(unsigned compl_state, const bdd& letter) {
        cola::tnba_complement::vec_state_taggedcol succs_B;
        
        // no succs yet
        if(compl_state_storage_.count(compl_state) == 0){
            compl_state_storage_.insert({compl_state,{}});
        }

        // find if there are defined
        for(auto &successors: compl_state_storage_[compl_state]) {
            bdd symb = successors.second;
            if(bdd_implies(letter, symb)) {
                succs_B = successors.first;
            }
        }
        // if not computed yet, compute
        if(succs_B.empty()) {
            // std::cout << std::to_string(compl_state) << " , " << std::to_string(letter) << "\n";
            const cola::tnba_complement::uberstate &us_B = aut_B_compl_.num_to_uberstate(compl_state);
            succs_B = aut_B_compl_.get_succ_uberstates(us_B, letter);
            compl_state_storage_[compl_state].emplace_back(std::pair(succs_B, letter));
        }

        return succs_B;
    }

    std::vector<std::shared_ptr<inclusion_mstate>> inclusion_check::get_succs(const std::shared_ptr<inclusion_mstate> &src) {
        std::vector<std::shared_ptr<inclusion_mstate>> cartesian_prod;
        std::vector<std::shared_ptr<inclusion_mstate>> result;

        auto casted_src = dynamic_cast<inclusion_mstate*>(src.get());

        // extract A state and compl.B state from intersection macrostate to compute successors
        unsigned state_of_A = casted_src->state_.first;
        unsigned compl_state = casted_src->state_.second;
        
        bdd all = n_s_compat_;
        struct BDDComparator {
            bool operator()(const bdd& a, const bdd& b) const {
                return a.id() < b.id();
            }
        };

        auto spot_s_A = aut_A_->state_from_number(state_of_A);
        std::map<bdd, std::set<unsigned>,BDDComparator> cond_to_states;
        for (auto i: aut_A_->succ(spot_s_A)) {
            if (cond_to_states.count(i->cond()) == 0)
                cond_to_states[i->cond()] = {};
        
            cond_to_states[i->cond()].emplace(aut_A_->state_number(i->dst()));
        }

        for (auto const& i: cond_to_states) {
            get_successors_compl(compl_state, i.first);
            std::set<unsigned> succs_A = i.second;
            if(succs_A.empty())
                continue;

            cola::tnba_complement::vec_state_taggedcol succs_B;
            if(!aut_B_compl_.get_is_sink_created() || compl_state != aut_B_compl_.get_sink_state())
            {
                succs_B = get_successors_compl(compl_state, i.first);
            }
            else
            {
                succs_B.push_back({aut_B_compl_.get_sink_state(), {}});
            }
            
            if(!succs_A.empty() && !succs_B.empty())
            {
                cartesian_prod = get_cartesian_prod(state_of_A, succs_A, succs_B,
                                                    i.first);
                for (auto& succ : cartesian_prod) {
                    result.emplace_back(std::move(succ));
                }
            }
        }
        return result;
    }

    bool inclusion_check::is_transition_acc(const spot::twa_graph_ptr &aut_A, unsigned src, unsigned dst, const bdd &symbol) {
        for (const auto &t: aut_A->out(src)) {
            if (dst == t.dst && bdd_implies(symbol, t.cond)) {
                if (t.acc) { return true; } 
                else return false;
            }
        }

        return false;
    }

    std::vector<std::shared_ptr<inclusion_mstate>>
    inclusion_check::get_cartesian_prod(unsigned aut_A_src, std::set<unsigned> &states_A,
                                        cola::tnba_complement::vec_state_taggedcol &states_B, const bdd &letter) {
        std::vector<std::shared_ptr<inclusion_mstate>> cartesian_prod;
        // COMPUTATION OF COLORS IS TAKEN FROM complement_sync.cpp might be better to create method in the complement_sync class

        // TODO the following two lines should not be called everytime this function is called !!!!
        auto vec_acc_cond = aut_B_compl_.get_vec_acc_cond();
        auto part_col_offset = aut_B_compl_.get_part_col_offset();

        for (const auto &state_A: states_A) {
            for (const auto &state_cols: states_B) {
                const unsigned &state_B = state_cols.first;
                auto cols = state_cols.second;
                std::set<unsigned> new_cols;

                if(aut_B_compl_.get_is_sink_created() && state_B == aut_B_compl_.get_sink_state()){
                    new_cols = infs_from_compl_; // TODO go to next iteration after the next if for A
                }
                if (is_transition_acc(aut_A_, aut_A_src, state_A, letter)) {
                    new_cols.insert(first_col_to_use_); // accepting mark of A
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


                auto uberstate = aut_B_compl_.num_to_uberstate(state_B);
                auto B_reach_set = uberstate.get_reach_set(); // (H,(C,S,B),...) -- B_reach_set is H
                std::set<unsigned> A_B_intersect;
                if(dir_simul_.count(state_A)) {
                    set_intersection(dir_simul_[state_A].begin(), dir_simul_[state_A].end(), B_reach_set.begin(), B_reach_set.end(),
                            std::inserter(A_B_intersect, A_B_intersect.begin()));
                }

                if(A_B_intersect.empty()) {
                    auto succ = std::make_shared<inclusion_mstate>();
                    succ->state_ = std::make_pair(state_A, state_B);
                    spot::acc_cond::mark_t spot_cols(new_cols.begin(), new_cols.end());
                    succ->acc_ = spot_cols;
                    succ->trans_cond_ = letter;
                    cartesian_prod.emplace_back(std::move(succ));

                    //auto it_bool_pair =
                    intersect_states_.insert({{state_A, state_B},
                                              {}});
                }
            }
        }

        return cartesian_prod;
    }
}
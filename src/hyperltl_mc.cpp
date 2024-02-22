#include "hyperltl_mc.h"
#include "inclusion_test.hpp"

namespace kofola {
    hyperltl_mc::hyperltl_mc(parsed_hyperltl_form_ptr parsed_hyperltl_f, spot::kripke_graph_ptr system):
    parsed_hyperltl_f_(parsed_hyperltl_f),
    system_(system),
    built_aut_(parsed_hyperltl_f->aut)
    {
        Quantification q;
        while(!parsed_hyperltl_f->q_list.empty()) {
            q = parsed_hyperltl_f->q_list.back();
            parsed_hyperltl_f->q_list.pop_back();

            if(q.type == static_cast<unsigned int>(QuantificationType::Exists)) {
                existential_projection(q.trace_vars);
            }
            else if(parsed_hyperltl_f->q_list.empty()) {
                // n-fold self composition and inclusion,
                // although even when there is only seq. of existential, we can
                // perform emptiness check on the fly
            }
            else {
                built_aut_ = kofola::complement_tela(built_aut_);
                existential_projection(q.trace_vars);
                built_aut_ = kofola::complement_tela(built_aut_);
            }
        }
    }

    spot::twa_graph_ptr hyperltl_mc::existential_projection(std::list<std::string> exist_trac_vars) {
        auto init_formula = parsed_hyperltl_f_->aut->get_init_state_number();
        auto init_system = system_->get_init_state_number();

        mc_macrostate initial_macrostate = std::make_pair(init_formula, init_system);
        std::queue<mc_macrostate> macrostates;
        macrostates.push(initial_macrostate);

        while(!macrostates.empty()) {

        }
    }


    AP_trace hyperltl_mc::parse_formula_AP(std::string input_ap) {
        auto pos = 0;
        unsigned curly_braces = 0;

        for(unsigned i = 0; i < input_ap.length(); i++) {
            if(input_ap[i] == '{')
                curly_braces++;
            else if(input_ap[i] == '}')
                curly_braces--;
            else if(input_ap[i] == '_' && curly_braces == 0) {
                pos = i;
                break;
            }
        }

        AP_trace res;
        res.atomic_prop = input_ap.substr(1, pos - 1); // exclude starting '{' and trailing '}_'
        res.trace_var = input_ap.substr(pos + 1); // exclude '_'

        return res;
    }
} // kofola
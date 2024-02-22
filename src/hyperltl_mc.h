#pragma once

#include "kofola.hpp"
#include "hyperltl_formula_processor.hpp"

namespace kofola {
    typedef struct {
        std::string atomic_prop;
        std::string trace_var;
    } AP_trace;

    class hyperltl_mc {
        using mc_macrostate = std::pair<unsigned, unsigned>; // first is system, second is formula

        parsed_hyperltl_form_ptr parsed_hyperltl_f_;
        spot::twa_graph_ptr built_aut_;
        spot::kripke_graph_ptr system_;
        std::map<std::string, AP_trace> ap_map_;
    public:
        hyperltl_mc(parsed_hyperltl_form_ptr parsed_hyperltl_f, spot::kripke_graph_ptr system);

        spot::twa_graph_ptr existential_projection(std::list<std::string> exist_trac_vars);

        AP_trace parse_formula_AP(std::string input_ap);
    };

} // kofola

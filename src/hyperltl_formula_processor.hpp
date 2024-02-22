#pragma once

// kofola
#include "kofola.hpp"
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/postproc.hh>
#include <spot/tl/formula.hh>
#include <spot/tl/parse.hh>
#include <spot/kripke/kripkegraph.hh>

namespace kofola {
    enum class QuantificationType : unsigned int {
        Exists = 0,
        Forall = 1
    };

    typedef struct {
        unsigned int type : 1;
        std::list<std::string> trace_vars; // for sequences of the same type
    } Quantification;

    struct parsed_hyperltl_form
    {
        spot::formula formula;
        spot::twa_graph_ptr aut;
        std::list<Quantification> q_list;
    };
    typedef std::shared_ptr<parsed_hyperltl_form> parsed_hyperltl_form_ptr;

    class hyperltl_formula_processor
    {
        std::string filename_;
        public:
            hyperltl_formula_processor(const std::string& file);

            void print_quantifications(const std::list<Quantification>& q_list);

            parsed_hyperltl_form_ptr parse_hyperltl_formula();
    };
}

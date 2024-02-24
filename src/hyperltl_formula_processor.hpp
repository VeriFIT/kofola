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
        unsigned int type : 1; /// 1-bit information of quantification type
        std::list<std::string> trace_vars; /// the order corresponds with the textual form on the input
    } Quantification;

    typedef struct {
        std::string atomic_prop; /// atomic proposition of the system
        std::string trace_var; /// of trace 'trace_var'
    } AP_trace;

    struct parsed_hyperltl_form
    {
        spot::formula formula; /// spot representation of the formula (if needed)
        spot::twa_graph_ptr aut; /// automaton for the formula
        std::list<Quantification> q_list; /// the order corresponds with the textual form on the input
        std::map<std::string, AP_trace> aps_map; /// to map formula automaton APs "{ap}_{trace}" -> (ap,trace) with ap being the system AP
    };
    typedef std::shared_ptr<parsed_hyperltl_form> parsed_hyperltl_form_ptr;

    class hyperltl_formula_processor
    {
        std::string filename_;
        public:
            hyperltl_formula_processor(const std::string& file);

            void print_quantifications(const std::list<Quantification>& q_list);

            /**
             * @brief For the string "{ap}_{trace}" creates struct 'AP_trace'
             *
             * */
            AP_trace parse_formula_AP(std::string input_ap);

            parsed_hyperltl_form_ptr parse_hyperltl_formula();
    };
}

/**
 * @file hyperltl_formula_processor.hpp
 * @author Ondrej Alexaj (xalexa09@stud.fit.vutbr.cz)
 * @brief Declarations for hyperltl formula processing
 * @version 0.1
 * @date 2024-05-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

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
        unsigned int type : 1; /// 1-bit information indicating quantification type
        std::vector<std::string> trace_vars; /// the order corresponds with the textual form on the input
    } Quantification;

    typedef struct {
        std::string atomic_prop; /// atomic proposition of the system
        std::string trace_var; /// of trace 'trace_var'
    } AP_trace;

    struct parsed_hyperltl_form
    {
        spot::formula formula; /// spot representation of the formula (if needed)
        spot::twa_graph_ptr aut; /// automaton for the formula
        unsigned qantifiers;
        std::list<Quantification> q_list; /// the order corresponds with the textual form on the input
        std::map<std::string, AP_trace> aps_map; /// to map formula automaton APs "{ap}_{trace}" -> (ap,trace) with ap being the system AP
        bool negate; /// before returning SAT or UNSAT, negate the answer (if negate = true)
    };
    typedef std::shared_ptr<parsed_hyperltl_form> parsed_hyperltl_form_ptr;

    class hyperltl_formula_processor
    {
        std::string filename_;
        public:
            /// constructor, takes file where is the hyperltl formula
            hyperltl_formula_processor(const std::string& file);

            /// to know if the negation should be on the top level, based on the quantifiers types and alternations
            void preprocess(parsed_hyperltl_form_ptr &formula_to_preproc);

            /// For the string "{ap}_{trace}" creates struct 'AP_trace'
            AP_trace parse_formula_AP(std::string input_ap);

            /// creates and fills structure with info about hyperltl formula
            parsed_hyperltl_form_ptr parse_hyperltl_formula();
    };
}

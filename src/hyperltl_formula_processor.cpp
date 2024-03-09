#include "hyperltl_formula_processor.hpp"

// Spot
#include <spot/twaalgos/postproc.hh>
#include <regex>

namespace kofola {
    hyperltl_formula_processor::hyperltl_formula_processor(const std::string& file) {
        filename_ = file;
    }

    void hyperltl_formula_processor::preprocess(parsed_hyperltl_form_ptr &formula_to_preproc) {
        if(formula_to_preproc->q_list.front().type == static_cast<unsigned int>(QuantificationType::Exists) &&
                formula_to_preproc->q_list.size() >= 2) {
            for (auto &q: formula_to_preproc->q_list) {
                if (q.type == static_cast<unsigned int>(QuantificationType::Exists))
                    q.type = static_cast<unsigned int>(QuantificationType::Forall);
                else
                    q.type = static_cast<unsigned int>(QuantificationType::Exists);
            }

            formula_to_preproc->negate = true;
        }
    }

    AP_trace hyperltl_formula_processor::parse_formula_AP(std::string input_ap) {
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
        res.atomic_prop = input_ap.substr(1, pos - 2); // exclude starting '{' and trailing '}_'
        res.trace_var = input_ap.substr(pos + 2, input_ap.size() - pos - 3); // exclude '_'

        return res;
    }

    parsed_hyperltl_form_ptr hyperltl_formula_processor::parse_hyperltl_formula() {
        std::string ltl_body;
        std::ifstream file(filename_);
        std::string quantifier;
        auto res = std::make_shared<parsed_hyperltl_form>();
        res->negate = false; // implicitly

        // cycle to collect each quantification with its trace variable
        while(file >> quantifier) {
            std::string trace_var;
            Quantification q;

            if(quantifier == "exists") {
                q.type = static_cast<unsigned int>(QuantificationType::Exists);
            }
            else if(quantifier == "forall") {
                q.type = static_cast<unsigned int>(QuantificationType::Forall);
            }
            else {
                // no more quantification, return ltl body
                ltl_body += quantifier; // quantifier holds some other string than quantifiers
                std::string tmp;

                preprocess(res);

                while(file >> tmp) {
                    ltl_body += tmp;
                }

                // search for APs
                std::regex pattern("\"([^\"]*)\"");
                std::smatch match;
                tmp = ltl_body;
                // Search for matches
                while (std::regex_search(tmp, match, pattern)) {
                    res->aps_map[match[1]] = parse_formula_AP(match[1]);
                    // Update input to search for the next occurrence
                    tmp = match.suffix();
                }

                if(res->negate)
                    ltl_body = "!(" + ltl_body + ")";
                break;
            }

            file >> trace_var;
            trace_var.pop_back(); // remove '.' following the trace variable
            // when alternation, keep order of quantifiers
            if(res->q_list.empty() || q.type != res->q_list.front().type) {
                q.trace_vars.emplace_back(trace_var);
                res->q_list.emplace_back(q);
            }
            else {
                res->q_list.front().trace_vars.emplace_back(trace_var);
            }
            res->qantifiers++;
        }

        res->formula = spot::parse_formula(ltl_body);
        spot::translator trans;
        trans.set_type(spot::postprocessor::Buchi);
        res->aut = trans.run(&(res->formula));

        return res;
    }
}
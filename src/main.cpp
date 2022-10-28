// Copyright (C) 2022  The Kofola Authors
//
// This file is a part of kofola, a tool for complementation of omega automata.
//
// cola is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// cola is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// kofola
#include "kofola.hpp"
#include "complement_tela.hpp"
#include "util.hpp"

// standard library headers
#include <unistd.h>
#include <ctime>
#include <string>

// spot
#include <spot/parseaut/public.hh>
#include <spot/misc/version.hh>

// Args.hxx
#include "../3rdparty/args.hxx"


//void output_input_type(spot::twa_graph_ptr aut)
//{
//  bool type = false;
//  if (spot::is_deterministic(aut))
//  {
//    type = true;
//    std::cout << "deterministic" << std::endl;
//  }
//  if (spot::is_semi_deterministic(aut))
//  {
//    type = true;
//    std::cout << "limit-deterministic" << std::endl;
//  }
//  if (cola::is_elevator_automaton(aut))
//  {
//    std::cout << "elevator" << std::endl;
//  }
//  if (cola::is_weak_automaton(aut))
//  {
//    std::cout << "inherently weak" << std::endl;
//  }
//  if (spot::is_unambiguous(aut))
//  {
//    std::cout << "unambiguous" << std::endl;
//  }
//  if (!type)
//  {
//    std::cout << "nondeterministic" << std::endl;
//  }
//}

void output_scc_info(spot::twa_graph_ptr aut)
{
  // strengther
  spot::scc_info si(aut, spot::scc_info_options::ALL);
  unsigned num_iwcs = 0;
  unsigned num_acc_iwcs = 0;
  unsigned num_iwcs_states = 0;
  unsigned num_max_iwcs_states = 0;
  unsigned num_acciwcs_states = 0;
  unsigned num_max_acciwcs_states = 0;
  unsigned num_dacs = 0;
  unsigned num_dacs_states = 0;
  unsigned num_max_dacs_states = 0;
  unsigned num_nacs = 0;
  unsigned num_nacs_states = 0;
  unsigned num_max_nacs_states = 0;

  std::string types = cola::get_scc_types(si);
  for (unsigned sc = 0; sc < si.scc_count(); sc++)
  {
    unsigned num = si.states_of(sc).size();
    if (cola::is_weakscc(types, sc))
    {
      num_iwcs_states += num;
      num_iwcs++;
      num_max_iwcs_states = std::max(num_max_iwcs_states, num);
    }
    if (cola::is_accepting_weakscc(types, sc))
    {
      num_acciwcs_states += num;
      num_acc_iwcs++;
      num_max_acciwcs_states = std::max(num_max_acciwcs_states, num);
    }

    if (cola::is_accepting_detscc(types, sc))
    {
      num_dacs_states += num;
      num_dacs++;
      num_max_dacs_states = std::max(num_max_dacs_states, num);
    }
    if (cola::is_accepting_nondetscc(types, sc))
    {
      num_nacs_states += num;
      num_nacs++;
      num_max_nacs_states = std::max(num_max_nacs_states, num);
    }
  }
  std::cout << "Number of IWCs: " << num_iwcs << " with " << num_iwcs_states << " states, in which max IWC with " << num_max_iwcs_states << " states\n";
  std::cout << "Number of ACC_IWCs: " << num_acc_iwcs << " with " << num_acciwcs_states << " states, in which max IWC with " << num_max_acciwcs_states << " states\n";
  std::cout << "Number of DACs: " << num_dacs << " with " << num_dacs_states << " states, in which max DAC with " << num_max_dacs_states << " states\n";
  std::cout << "Number of NACs: " << num_nacs << " with " << num_nacs_states << " states, in which max NAC with " << num_max_nacs_states << " states\n";
}

namespace
{ // {{{

void print_version()
{
	std::cout << "VERSION\n";
	exit(EXIT_SUCCESS);
}

void print_version_long()
{
	std::cout << "VERSION\n";
	assert(false);
	exit(EXIT_SUCCESS);
}


/// parse string with parameters of the form "key=value;flag;..."
kofola::string_to_string_dict parse_params(const std::string& str)
{ // {{{
	DEBUG_PRINT_LN("parameter string: " + str);

	kofola::string_to_string_dict res;
	std::string tmp_str = kofola::str_trim(str);
	while (!tmp_str.empty()) {
		auto delim_pos = tmp_str.find(";");
		std::string one_param = kofola::str_trim(tmp_str.substr(0, delim_pos));
		if (std::string::npos == delim_pos) {
			tmp_str = "";
		} else {
			tmp_str = kofola::str_trim(tmp_str.substr(delim_pos + 1));
		}

		DEBUG_PRINT_LN("parameter: " + one_param);

		std::string key;
		std::string value;
		auto eq_pos = one_param.find("=");
		if (std::string::npos == eq_pos) {
			key = one_param;
			value = "";
		} else {
			key = kofola::str_trim(one_param.substr(0, eq_pos));
			value = kofola::str_trim(one_param.substr(eq_pos + 1));
		}

		DEBUG_PRINT_LN("parameter key = " + key + ", value = " + value);
		auto it_bool_pair = res.insert({key, value});
		if (!it_bool_pair.second) {
			throw std::runtime_error("an attempt to redefine parameter " + key + " with " + value);
		}
	}

	return res;
} // parse_params() }}}


/**
 * process command line arguments
 * @param[in]  argc  Number of arguments
 * @param[in]  argv  Array of arguments
 * @param[out] params  Output structure with arguments
 * @returns  EXIT_SUCCESS  iff successful  */
int process_args(int argc, char *argv[], kofola::options* params)
{ // {{{
	assert(nullptr !=params);

	args::ArgumentParser parser("kofola: modular complementation of omega-automata.");
	args::CompletionFlag completion(parser, {"complete"});

	// inputs
	args::Group input_group(parser, "Input:");
	args::PositionalList<std::string> filenames(input_group, "FILENAME",
	                                            "files with input omega automata "
	                                            "(in the HOA format or other formats that Spot can process)");

	// output form
	args::Group output_type_group(parser, "Output automaton type:", args::Group::Validators::AtMostOne);
	args::Flag buchi_flag(output_type_group, "buchi", "transition-based Buchi automaton", {'b', "buchi", "tba"});
	args::Flag tgba_flag(output_type_group, "tgba", "transition-based generalized Buchi automaton", {"tgba"});
	args::Flag tela_flag(output_type_group, "tela", "transition-based Emerson-Lei automaton", {"tela"});

	// command
	args::Group operation_group(parser, "Operation:", args::Group::Validators::AtMostOne);
	args::Flag complement_flag(operation_group, "complement", "complement the inputs (default)", {"complement"});
	args::Flag type_flag(operation_group, "type", "print out types of the inputs", {"type"});
	args::Flag scc_types_flag(operation_group, "scc-types", "print out types of SCCs in the inputs", {"scc-types"});
	args::ActionFlag version_flag(operation_group, "version", "print program version", {"version"}, print_version);
	args::ActionFlag version_long_flag(operation_group, "version-long", "print program version (long)", {"version-long"}, print_version_long);
	args::HelpFlag help_flag(operation_group, "help", "display this help menu", {'h', "help"});

	// miscellaneous flags
	args::Group misc_group(parser, "Miscellaneous options:");
	args::Flag debug_flag(misc_group, "debug", "output debugging information", {'d', "debug"});
	args::ValueFlag<std::string> params_flag(misc_group, "params",
	                                         "string with ';'-separated parameters of the form 'key=value' "
	                                         "(or just 'key' for yes/no flags), e.g., "
	                                         "'merge-iwa=yes;preproc-reduction=high;raw-output'",
	                                         {"params"});

	try {
		parser.ParseCLI(argc, argv);
	}
	catch (const args::Completion& e) {
		std::cout << e.what();
		return EXIT_SUCCESS;
	}
	catch (const args::Help&) {
		std::cout << parser;
		return EXIT_SUCCESS;
	}
	catch (const args::ValidationError& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return EXIT_FAILURE;
	}
	catch (const args::ParseError& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return EXIT_FAILURE;
	}

	if (!filenames) { // no filename
		if (isatty(STDIN_FILENO)) { // stdin is terminal
			std::cerr << "kofola: No automaton to process? Run 'kofola --help' for more help.\n";
			std::cerr << "Use 'kofola -' to force reading automata from the standard input.\n";
			return EXIT_FAILURE;
		} else { // stdin is e.g. piped
			params->filenames.emplace_back("-");  // "-" denotes stdin
		}
	} else {
		params->filenames = filenames.Get();
	}

	if (type_flag) {
		params->operation = "type";
	} else if (scc_types_flag) {
		params->operation = "scc-types";
	} else if (help_flag) {
		params->operation = "help";
	} else { // default
		params->operation = "complement";
	}

	if (buchi_flag) {
		params->output_type = "buchi";
	} else if (tgba_flag) {
		params->output_type = "tgba";
	} else { // default
		params->output_type = "tela";
	}

	if (debug_flag) {
		kofola::LOG_VERBOSITY = 42;
	}

	if (params_flag) {
		params->params = parse_params(params_flag.Get());
	}

	return EXIT_SUCCESS;
} // process_args() }}}

} // anonymous namespace }}}


/// entry point
int main(int argc, char *argv[])
{ // {{{
	kofola::options options;
	int rv = process_args(argc, argv, &options);
	if (EXIT_SUCCESS != rv) { return EXIT_FAILURE; }
	kofola::OPTIONS = options;   // set the global variable

	DEBUG_PRINT_LN("filenames: " + std::to_string(options.filenames));
	DEBUG_PRINT_LN("operation: " + std::to_string(options.operation));
	DEBUG_PRINT_LN("params: " + std::to_string(options.params));

	auto dict = spot::make_bdd_dict();

	for (const std::string& input_filename : options.filenames) {
		spot::parsed_aut_ptr parsed_aut = nullptr;
		try {
			spot::automaton_stream_parser parser(input_filename);

			while (true) { // keep reading from the file as long as possible
				parsed_aut = parser.parse(dict);
				if (parsed_aut->format_errors(std::cerr)) { return EXIT_FAILURE; }

				// input automaton
				spot::twa_graph_ptr aut = parsed_aut->aut;
				if (!aut) { break; }

				if (options.operation == "complement") {
					//clock_t c_start = clock();
					spot::twa_graph_ptr result = kofola::complement_tela(aut);
					//clock_t c_end = clock();
					//auto duration = (c_end - c_start) / CLOCKS_PER_SEC;

					spot::print_hoa(std::cout, result);
					std::cout << "\n";
				} else if (options.operation == "type") {
					assert(false);
				} else if (options.operation == "scc-types") {
					assert(false);
				} else {
					throw std::runtime_error("invalid operation: " + options.operation);
				}
			}
		}
		catch (const std::exception& ex) {
			std::cerr << "Error: " << ex.what() << "\n";
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;

  // options
//  bool aut_type = false;
//  bool print_scc = false;
//
//  postprocess_level preprocess = Low;
//  postprocess_level post_process = Low;
//
//  output_aut_type output_type = Generic;

//  for (int i = 1; i < argc; i++)
//  {
//    std::string arg = argv[i];
//    if (arg.find("--preprocess=") != std::string::npos)
//    {
//      unsigned level = 0;//parse_int(arg);
//      if (level == 0)
//      {
//        preprocess = None;
//      }
//      else if (level == 1)
//      {
//        preprocess = Low;
//      }
//      else if (level == 2)
//      {
//        preprocess = Medium;
//      }
//      else if (level == 3)
//      {
//        preprocess = High;
//      }
//    }
//    else if (arg == "--print-scc")
//    {
//      print_scc = true;
//    }
//    else if (arg == "--debug") {
//      kofola::LOG_VERBOSITY = 42;
//    }
//    else if (arg == "--type") {
//        aut_type = true;
//    }
//    else if (arg == "--low-red-interm") {
//      decomp_options.low_red_interm = true;
//    }
//    else if (arg == "--merge-iwa")
//    {
//      decomp_options.merge_iwa = true;
//    }
//    else if (arg == "--merge-det")
//    {
//      decomp_options.merge_det = true;
//    }
//    else if (arg == "--tgba")
//    {
//      decomp_options.tgba = true;
//    }
//    else if (arg == "--tba")
//    {
//      decomp_options.tba = true;
//    }
//    else if (arg == "--raw")
//    {
//      decomp_options.raw = true;
//    }
//    else if (arg == "--rank")
//    {
//      decomp_options.rank_for_nacs = true;
//    }
//    else if (arg == "--iw-sim")
//    {
//      decomp_options.iw_sim = true;
//    }
//    else if (arg == "--det-sim")
//    {
//      decomp_options.det_sim = true;
//    }
//    else if (arg == "--scc-compl")
//    {
//      decomp_options.scc_compl = true;
//    }
//    else if (arg == "--scc-high")
//    {
//      decomp_options.scc_compl_high = true;
//    }
//    else if (arg == "--no-sat")
//    {
//      decomp_options.sat = false;
//    }
//    else if (arg == "--dataflow")
//    {
//      decomp_options.dataflow = true;
//    }
//    else if (arg == "--version")
//    {
//      std::cout << "kofola x.y"
//                   " (using Spot "
//                << spot::version() << ")\n\n"
//                                      "Copyright (C) 2020  The cola Authors.\n"
//                                      "License GPLv3+: GNU GPL version 3 or later"
//                                      " <http://gnu.org/licenses/gpl.html>.\n"
//                                      "This is free software: you are free to change "
//                                      "and redistribute it.\n"
//                                      "There is NO WARRANTY, to the extent permitted by law.\n"
//                << std::flush;
//      return 0;
//    }
//  }


//  for (std::string &path_to_file : path_to_files)
//  {
//    spot::automaton_stream_parser parser(path_to_file);
//
//    for (;;)
//    {
//      spot::parsed_aut_ptr parsed_aut = parser.parse(dict);
//
//      if (parsed_aut->format_errors(std::cerr))
//        return 1;
//
//      // input automata
//      spot::twa_graph_ptr aut = parsed_aut->aut;
//
//      if (!aut)
//        break;
//
//      // Check if input is TGBA
////      if (aut->acc().is_generalized_buchi())
////      {
////        aut = spot::degeneralize_tba(aut);
////      }
//
//      if (!aut->acc().is_buchi())
//      {
//        std::cerr << "cola requires Buchi condition on input.\n";
//        return 1;
//      }
//
//      if (aut_type)
//      {
////        output_input_type(aut);
//        break;
//      }
//
//      if (print_scc)
//      {
//        output_scc_info(aut);
//        break;
//      }
//    }
//  }
} // main() }}}

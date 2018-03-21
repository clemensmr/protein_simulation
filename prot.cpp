 
/*
 * main.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include <string>
#include <vector>
#include <time.h>

#include "protdata.hpp"
#include "args.hxx"

//using namespace Ipopt;


int main(int argc, char** argv) {
    bool dim2 = false;


    int MPC_horizon = 100;


    string matrix_A = "../A.mtx";


    args::ArgumentParser parser("convection diffusion equation 1d.", "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu",{'h', "help"});
    args::CompletionFlag completion(parser,{"complete"});

    args::Flag dim2_(parser, "2 dimensional", "turn dimension 2 on",{'d', "dim2"});

    args::ValueFlag<int> MPC_horizon_(parser, "mpc horizon", "MPC horizon, N >= 1",{'N', "mpc"});


    try {
	parser.ParseCLI(argc, argv);
    }
    catch (args::Completion e) {
	std::cout << e.what();
	return 0;
    }
    catch (args::Help) {
	std::cout << parser;
	return 0;
    }
    catch (args::ParseError e) {
	std::cerr << e.what() << std::endl;
	std::cerr << parser;
	return 1;
    }

    //int
    if (MPC_horizon_) {
	MPC_horizon = args::get(MPC_horizon_);
    }
    
    std::cout << "this works" << std::endl;

}




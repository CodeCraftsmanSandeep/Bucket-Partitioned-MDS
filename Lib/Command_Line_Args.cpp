#include "Command_Line_Args.h"
#include "Utils.h"

#include <string> 
#include <memory>
#include <fstream> 
#include <iostream>
#include <filesystem> 

std::string Command_Line_Args::get_usage()
{
    /*
    * get_usage: Helper function to get how to run the executable 
    * @return: Returns string which contains usage information of the executable
    */

    std::ostringstream usage;

    usage << "Allowed arguments:\n"
          << "  --alpha   : Angle to partition the 2-D plane\n"
          << "  --rho     : Number of different DFS orderings for exploitation\n"
          << "  --input   : Input file path\n"
          << " --output   : Output file path\n\n"

          << "Mandatory arguments:\n"
          << "  alpha and rho\n\n"

          << "Optional arguments:\n"
          << "  input and output\n\n"

          << "Acceptable usage examples:\n"
          << "  ./executable --alpha=<alpha> --rho=<rho>\n"
          << "  ./executable --alpha=<alpha> --rho=<rho> --output=<output>\n"
          << "  ./executable --alpha=<alpha> --rho=<rho> --input=<input>\n"
          << "  ./executable --alpha=<alpha> --rho=<rho> --input=<input> --output=<output>\n";

    return usage.str();
}


void Command_Line_Args::set_alpha(int argc, char** argv)
{
    /* 
    * set_alpha: Helper function to get alpha from command line arguments and store it in data attribute
    * @param argc: Number of command-line arguments 
    * @param argv: Array of command-line arguments
    * @return: Returns nothing
    */

    bool alpha_passed = false;

    for(int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg.rfind("--alpha=", 0) == 0) 
        {
            if (alpha_passed)
            {
                HANDLE_ERROR("Alpha must be passed only once in command line arguments\n" + this->get_usage(), true);
            }
            alpha_passed = true;
            this->alpha = std::stod(arg.substr(8));
            if (this->alpha <= 0 || this->alpha > 360)
            {
                HANDLE_ERROR("Alpha must be in the range (0, 360]", true);
            }
        } 
    }
    
    if(!alpha_passed)
    {
        HANDLE_ERROR("Alpha is not passed\n" + this->get_usage(), true);
    }

    return;
}

void Command_Line_Args::set_rho(int argc, char** argv)
{
    /* 
    * set_rho: Helper function to get rho from command line arguments and store it in data attribute
    * @param argc: Number of command-line arguments 
    * @param argv: Array of command-line arguments
    * @return: Returns nothing
    */

    bool rho_passed = false;

    for(int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg.rfind("--rho=", 0) == 0) 
        {
            if (rho_passed)
            {
                HANDLE_ERROR("Rho must be passed exactly once in command line arguments\n" + this->get_usage(), true);
            }
            rho_passed = true;
            this->rho = std::stoi(arg.substr(6));
            if (this->rho <= 0)
            {
                HANDLE_ERROR("Rho must be positive", true);
            }
        } 
    }

    if(!rho_passed)
    {
        HANDLE_ERROR("Rho is not passed\n" + this->get_usage(), true);
    }

    return;
}

void Command_Line_Args::set_input_stream(int argc, char** argv)
{
    /* 
    * set_input_stream: Helper function to set up the input stream from command-line arguments, 
                        defaults to std::cin if no input file is provided
    * @param argc: Number of command-line arguments 
    * @param argv: Array of command-line arguments
    * @return: Returns nothing
    */

    input_stream = &std::cin; // Default value
    bool input_file_passed = false;

    for (int i = 1; i < argc; i++) 
    {
        std::string arg = argv[i];
        if (arg.rfind("--input=", 0) == 0) 
        {
            if(input_file_passed)
            {
                HANDLE_ERROR("Input file path must be passed at most once\n" + this->get_usage(), true);
            }
            input_file_passed = true;
            std::string filename = arg.substr(8);
            input_file_owner = std::make_unique<std::ifstream>(filename);
            if (!input_file_owner->is_open())
            {
                HANDLE_ERROR("Cannot open input file: " + filename, true);
            }
            input_stream = input_file_owner.get();
        } 
    }
}

void Command_Line_Args::set_output_stream(int argc, char** argv)
{
    /* 
    * set_output_stream: Helper function to set up the output stream from command-line arguments, 
                        defaults to std::cout if no output file is provided
    * @param argc: Number of command-line arguments 
    * @param argv: Array of command-line arguments
    * @return: Returns nothing
    */

    output_stream = &std::cout; // Default value
    bool output_file_passed = false;
    for (int i = 1; i < argc; i++) 
    {
        std::string arg = argv[i];
        if (arg.rfind("--output=", 0) == 0) 
        {
            if(output_file_passed)
            {
                HANDLE_ERROR("Output file path must be passed at most once\n" + this->get_usage(), true);
            }
            output_file_passed = true;
            std::string filename = arg.substr(9);
            std::filesystem::path out_path(filename);

            // create directories if needed
            if (out_path.has_parent_path()) 
            {
                std::filesystem::create_directories(out_path.parent_path());
            }
            
            output_file_owner = std::make_unique<std::ofstream>(out_path, std::ios::trunc);
            if (!output_file_owner->is_open()) 
            {
                HANDLE_ERROR("Cannot open output file: " + filename, true);
            }
            output_stream = output_file_owner.get();
        } 
    }
}

Command_Line_Args::Command_Line_Args(int argc, char** argv) 
{
    /* 
    * Command_Line_Args: Constructor parses the command-line arguments and initializes the class members 
    * @param argc: Number of command-line arguments
    * @param argv: Array of command-line arguments 
    */

    // Parsing mandatory function parameters 
    this->set_alpha(argc, argv);
    this->set_rho(argc, argv);

    // Parsing optional function parameters
    this->set_input_stream(argc, argv);
    this->set_output_stream(argc, argv);
}

std::istream& Command_Line_Args::input() const 
{ 
    /*
    * input: Function to access input stream
    */

    return *input_stream; 
}

std::ostream& Command_Line_Args::output() const
{ 
    /*
    * output: Function to access output stream
    */

    return *output_stream; 
}

double Command_Line_Args::get_alpha() const 
{ 
    /*
    * get_alpha: Function to get alpha value
    */

    return alpha; 
}

int Command_Line_Args::get_rho() const 
{ 
    /*
    * get_rho: Function to get rho value
    */

    return rho; 
}
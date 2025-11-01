#ifndef COMMAND_LINE_ARGS_H
#define COMMAND_LINE_ARGS_H

#include <string>
#include <istream>
#include <ostream>
#include <memory>
#include <fstream>  
#include <iostream>

class Command_Line_Args 
{
    /*
    * Command_Line_Args: Class for managing:
        i) parsing command line arguments
        ii) managing I/O
    */
private:
    std::unique_ptr<std::ifstream> input_file_owner;
    std::unique_ptr<std::ofstream> output_file_owner;
    std::istream* input_stream; // points to either file or std::cin
    std::ostream* output_stream; // points to either file or std::cout
    double alpha; // in degrees
    int lambda;
    int rho;
 
    void set_input_stream(int, char**);
    void set_output_stream(int, char**);
    void set_alpha(int, char**);
    void set_lambda(int, char**);
    void set_rho(int, char**);
    std::string get_usage();
public:
    Command_Line_Args(int, char**);
    // Accessors
    std::istream& input() const;
    std::ostream& output() const;
    // Getters
    double get_alpha() const;
    int get_lambda() const; 
    int get_rho() const;
};

#endif
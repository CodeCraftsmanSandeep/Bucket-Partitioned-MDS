#include "Utils.h"
#include <iostream>

std::ostream& ERROR_FILE = std::cout;
void handle_error(
    const char*       file,
    int               line,
    std::ostream&     out,
    const std::string msg,
    const bool        exit_flag) 
{
    out << "❌ Error: " << msg << "\n"
        << "   In file: " << file << "\n"
        << "   At line: " << line << std::endl;

    if(exit_flag)
    {
        std::exit(EXIT_FAILURE);
    }
}
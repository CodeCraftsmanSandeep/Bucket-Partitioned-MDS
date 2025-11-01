#include "Command_Line_Args.h"
#include <iomanip>
#include <omp.h>

void initialize(const Command_Line_Args& args)
{
    args.output() << std::fixed << std::setprecision(4);
    omp_set_nested(3); // Enable nested parallelsim
}
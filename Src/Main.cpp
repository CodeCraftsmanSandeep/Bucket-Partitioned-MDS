#include "Command_Line_Args.h"
#include "Initializer.h" 
#include "Bucket_Partitioned_MDS.h"

int main(int argc, char* argv[])
{ 
    // Parse command line arguments 
    const Command_Line_Args args(argc, argv);

    // Initialize the process execution 
    initialize(args);

    // Create CVRP object
    const Bucket_Partitioned_MDS::CVRP cvrp(args.input());

    // Create Bucket-Partitioned-MDS solver
    const Bucket_Partitioned_MDS::Solver solver(args.get_alpha(), args.get_rho());

    // Solve CVRP using Bucket-Partitioned-MDS solver
    Bucket_Partitioned_MDS::Solution cvrp_solution = solver.solve(cvrp);

    // Verify the solution
    cvrp_solution.verify(cvrp);
 
    // Print the solution
    cvrp_solution.print(args.output());

    // Return success
    return 0;
}
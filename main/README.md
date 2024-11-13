FORWARD

FORWARD algorithm is composed by the succession of 4 sub-process:
 1. PreProcessor
 2. Islander
 3. NetConcad
 4. Sampler
    
Which together generate a feasible polytree as a solution.

The algorithm can be run just using Forward.jl function with the indicated inputs or breaking it into the four sub-processes as it is shown in demo/RandomSmallGraphExample.ipynb.

Auxiliary.jl file containts auxiliary functions for these four main processes.

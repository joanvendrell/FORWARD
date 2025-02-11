"""
######################## Conditioned Markov Chain #########################
Select integer variables as a random walk in the solution space. There is
no need to deactivate some constraints in order to compute a solution.
"""
module ConditionedMarkovChain

include("../forward/Forward.jl")
include("../forward/Auxiliary.jl")
using .Forward
using .Auxiliary
using Graphs
using MetaGraphs
using Combinatorics

function random_walk(G::MetaDiGraph{Int64, Float64},
                     Omega::Vector{Int},
                     kappa::Int,
                     optimizer,
                     strategy::String,
                     limit_iters = 100,
                     outlev = 0)
    # define the maximum number of iterations
    limit_iters = min(limit_iters,((kappa*length(Omega))/(length(Omega)-4*(kappa-1)))*(2+log(kappa)))
    # prepare the search space
    V_g = []; S = nothing; f_S = Inf; termination_status = 0
    solution_space = collect(combinations(Omega, kappa))
    # find an initial point
    while termination_status==0
        V_g = rand(solution_space)       
        G_o = deepcopy(G)
        for non_candidate in setdiff(Omega,V_g)
            set_prop!(G_o,non_candidate,:p,0)
        end
        S, termination_status, f_S = Forward.forward(G_o,optimizer,strategy,outlev)   
        filter!(x -> x != V_g, solution_space)
    end
    # start iterations
    termination_condition = 0; iters = 0
    while (termination_condition==0)&&(iters<=limit_iters)
        V_g_i = rand(solution_space)
        G_i = deepcopy(G)
        for non_candidate in setdiff(Omega,V_g_i)
            set_prop!(G_i,non_candidate,:p,0)
        end
        S_i, termination_status_i, f_S_i = Forward.forward(G_i,optimizer,strategy,outlev)   
        filter!(x -> x != V_g_i, solution_space)
        # update solution
        if f_S_i <= f_S && termination_status_i == 1
            S = S_i
            f_S = f_S_i
            termination_status = termination_status_i
            V_g = V_g_i
        end
        # check loop solution
        iters = iters + 1
        if isempty(solution_space)
            termination_condition = 1
        end
    end
    return V_g, S, f_S, termination_status
end

end

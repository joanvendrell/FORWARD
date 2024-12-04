"""
######################## Conditioned Greedy #########################
Select integer variables with greedy algorithm. It is important to
deactivate some constraints in order to compute a solution.
"""
module ConditionedGreedy

include("../forward/Forward.jl")
include("../forward/Auxiliary.jl")
using .Forward
using .Auxiliary
using Graphs
using MetaGraphs

function greedy(G::MetaDiGraph{Int64, Float64},
                Omega::Vector{Int},
                kappa::Int,
                optimizer,
                strategy::String,
                outlev = 0)
    V_g = []; S = nothing; f_S = Inf; termination_status = 0
    while length(V_g)< kappa
        candidates = setdiff(Omega,V_g)
        V_g_t = V_g; S_t = nothing; f_S_t = Inf; termination_status_t = 0
        for x in candidates
            V_g_t_i = convert(Vector{Int}, V_g ∪ x)
            # rearrange graph
            G_i = deepcopy(G)
            for non_candidate in setdiff(candidates,V_g_t_i)
                set_prop!(G_i,non_candidate,:p,0)
            end
            # compute error
            S_i, termination_status_i, f_S_i = Forward.forward(G_i,optimizer,strategy,outlev)
            if ((length(V_g)==kappa-1)&&(termination_status_i==1)&&(f_S_i<=f_S_t))||((length(V_g)<kappa-1)&&(f_S_i<=f_S_t))
                V_g_t = V_g_t_i
                S_t = S_i
                f_S_t = f_S_i
                termination_status_t = termination_status_i
            end
        end
        # update values
        V_g = V_g_t
        S = S_t
        f_S = f_S_t
        termination_status = termination_status_t
    end
    return V_g, S, f_S, termination_status
end

end
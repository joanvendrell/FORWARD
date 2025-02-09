"""
######################## Custom Problem #########################
This file implements a custom optimization problem.
"""
module CustomProblem

include("CustomVariables.jl")
include("CustomConstraints.jl")
include("CustomObjective.jl")
include("CustomWarmUp.jl")
using .CustomVariables
using .CustomConstraints
using .CustomObjective
using .CustomWarmUp
using MetaGraphs
using Graphs
using PowerModelsDistribution
const _PMD = PowerModelsDistribution

# Implementation with warm-start from FORWARD solver
function custom_problem_from_FORWARD(pm::_PMD.AbstractUnbalancedPowerModel, 
                                     k::Int=2, 
                                     S::MetaDiGraph{Int64, Float64}=nothing,
                                     PMDtoFOR::Dict=nothing, 
                                     FORtoPMD::Dict=nothing,
                                     warmup::Bool=false)
    CustomVariables.import_variables_transmission_loss(pm)
    CustomConstraints.import_constraints_transmission_loss(pm,k)
    CustomObjective.custom_objective_transmission_loss_relaxed(pm)
    if warmup
        CustomWarmUp.custom_warmstart_from_FORWARD_loss(pm,S,PMDtoFOR,FORtoPMD)
    end
end

# Implementation with warm-start from MINLP solver
function custom_problem_transmission_loss(pm::_PMD.AbstractUnbalancedPowerModel, k::Int=2, data::Dict=nothing, warmup::Bool=false)
    CustomVariables.import_variables_transmission_loss(pm)
    CustomConstraints.import_constraints_transmission_loss(pm,k)
    CustomObjective.custom_objective_transmission_loss_relaxed(pm)
    if warmup
        CustomWarmUp.custom_warmstart_transmission_loss(pm,data)
    end
end

end

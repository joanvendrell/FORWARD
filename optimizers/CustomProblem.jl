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
using PowerModelsDistribution
const _PMD = PowerModelsDistribution

function custom_problem_transmission_loss(pm::_PMD.AbstractUnbalancedPowerModel, k::Int=2, data::Dict=nothing, warmup::Bool=false)
    CustomVariables.import_variables_transmission_loss(pm)
    CustomConstraints.import_constraints_transmission_loss(pm,k)
    CustomObjective.custom_objective_transmission_loss_relaxed(pm)
    if warmup
        CustomWarmUp.custom_warmstart_transmission_loss(pm,data)
    end
end

end
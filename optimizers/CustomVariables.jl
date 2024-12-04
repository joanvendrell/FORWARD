"""
######################## Custom Power Variables #########################
This file implements a function which imports the optimization 
variables for an specific MILNP defined.
"""
module CustomVariables

using PowerModelsDistribution
import PowerModelsONM
import InfrastructureModels
import JuMP
const _ONM = PowerModelsONM
const _PMD = PowerModelsDistribution
const _IM = InfrastructureModels


function import_variables_transmission_loss(pm::_PMD.AbstractUnbalancedPowerModel)
    _PMD.variable_mc_bus_voltage(pm)                     # creates voltage variables  :va
    _PMD.variable_mc_branch_power(pm)                    # creates power flow variables for all power lines :p
    _PMD.variable_mc_load_indicator(pm, relax=true)      # creates a variable for each load :z_demand
    _PMD.variable_mc_load_power(pm)                      # creates power input variables for all loads :pd
    _PMD.variable_mc_generator_power(pm, bounded=false)  # creates variables for all generators :pg
    _PMD.variable_mc_switch_power(pm)                    # creates power flow variables for all power lines :psw
    _PMD.variable_mc_switch_state(pm, relax=false)       # creates a variable state for each switch :switch_state
    _PMD.variable_mc_gen_indicator(pm, relax=false)      # creates a variable state for each generator :z_gen
    y_variable(pm)
    z_variable(pm) 
    println("Variables defined!")
end

# Auxiliary variable for Line Loss
function y_variable(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Our variable should be the edges ("switch")
    switches = Dict(i => switch["index"] for (i,switch) in ref(pm, nw, :switch))
    # ref(pm,nw,:bus) selects from our PowerModel pm the buses being nw the id of each one of them
    # this is the convention we use to create variables to make accessing them easier
    y = var(pm, nw)[:y] = Dict(i => JuMP.@variable(pm.model,
                                                   [switches[i]], base_name="$(nw)_y_$(i)",
                                                   #start = comp_start_value(ref(pm, nw, :switch, i),
                                                   #["y_start", "y"], switches[i], 1.0)
            ) for i in ids(pm, nw, :switch))
    # syntax used for reporting the variable solution
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :switch, :y, _PMD.ids(pm, nw, :switch), y)
    println("Power loss variable (y_variable) created!")
end

# Auxiliary variable for relaxations
function z_variable(pm::_PMD.AbstractUnbalancedPowerModel; nw::Int=_PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Our variable should be the edges ("switch")
    switches = Dict(i => switch["index"] for (i,switch) in ref(pm, nw, :switch))
    # ref(pm,nw,:bus) selects from our PowerModel pm the buses being nw the id of each one of them
    # this is the convention we use to create variables to make accessing them easier
    z = var(pm, nw)[:z] = Dict(i => JuMP.@variable(pm.model,
                                                   [switches[i]], base_name="$(nw)_z_$(i)",
                                                   #start = comp_start_value(ref(pm, nw, :switch, i),
                                                   #["z_start", "z"], switches[i], 1.0)
                ) for i in ids(pm, nw, :switch))
    # syntax used for reporting the variable solution
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :switch, :z, _PMD.ids(pm, nw, :switch), z)
    println("Relax variable variable (z) created!")
end

end
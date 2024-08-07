module KerrQuasinormalModes

using LinearAlgebra
using NLsolve
using Statistics
using CSV
using StaticArrays

using Parameters
using BSplineKit



include(joinpath(@__DIR__, "ModeSolver" ,"Angular.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"Radial.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"RootSolver.jl"))
include(joinpath(@__DIR__, "ModeSolver/Schwarszchild" ,"Schwarszchild.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"SpinSequenceOptions.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"SpinSequence.jl"))


include(joinpath(@__DIR__, "ModeFunctionInterface" ,"SpinWeightedSphericalLookup.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"Interface.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"ConfluentHeun.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"LinearCombinations.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"Derivative.jl"))


export SpinWeightedSpherical
export HeunConfluentRadial
export qnmfunction
export QuasinormalModeFunction
export ComplexPlot
export CP
export ModeSequence

#For debugging
export HeunConfluentLocalExpansion, ComputeSeries, qnm, Custom
#export importqnm, qnmfunctionnew
end

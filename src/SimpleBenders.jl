module SimpleBenders

using JuMP
const MOI = JuMP.MOI

include("subproblem.jl")
include("master_problem.jl")

end # module

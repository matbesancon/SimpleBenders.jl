module SimpleBenders

using JuMP
const MOI = JuMP.MOI
import Dualization

include("subproblem.jl")
include("master_problem.jl")

end # module

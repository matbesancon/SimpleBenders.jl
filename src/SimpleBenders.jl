module SimpleBenders

using JuMP
import SCIP
import Clp
const MOI = JuMP.MOI

include("subproblem.jl")
include("master_problem.jl")

end # module

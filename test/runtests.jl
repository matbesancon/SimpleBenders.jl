using Test
include("../src/SimpleBenders.jl")
import .SimpleBenders
using JuMP

using Gurobi
using LinearAlgebra



# test from http://www.iems.ucf.edu/qzheng/grpmbr/seminar/Yuping_Intro_to_BendersDecomp.pdf

function test_data()
    c = [2., 3.]
    A = [1 2;2 -1]
    D = zeros(2, 1) .+ [1, 3]
    b = [3, 4]
    return SimpleBenders.SubProblemData(b, D, A, c)
end

function test_result()
    d = test_data()
    m = Model(
        optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0)
    )
    @variable(m, x[1:2] >= 0)
    @variable(m, y[1:1] >= 0)
    @objective(m, Min, d.c'*x + 2y[1])
    @constraint(m, d.A * x + d.D * y .>= d.b)
    optimize!(m)
    return (JuMP.value.(x), JuMP.value.(y), JuMP.objective_value(m))
end

# test2 from Computational Example 3.2 (Subproblem Infeasibility) of Conejo A J, Castillo E, Minguez R, et al. Decomposition techniques in mathematical programming: engineering and science applications[M]. Springer Science & Business Media, 2006.
function test2_data()
    c = vec([-2., -1., 1.])
    A = - [1. 0. 0.; 0. 2. 0.; 0. 0. 1.; 0. 0. 0.]
    D = - [1. 1.; 3. 0.; 0. -7.; -1. 1.]
    b = - vec([3., 12., -16., 2.])
    return SimpleBenders.SubProblemData(b, D, A, c)
end

function test2_result()
    d = test2_data()
    m = Model(
        optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0)
    )
    @variable(m, x[1:3] >= 0)
    @variable(m, y[1:2] >= 0)
    @objective(m, Min, d.c'*x + 3y[1] - 3y[2])
    @constraint(m, d.A * x + d.D * y .>= d.b)
    optimize!(m)
    return (JuMP.value.(x), JuMP.value.(y), JuMP.objective_value(m))
end

@testset "Basic test" begin
    data = test_data()
    f(v) = 2v[1]
    m = Model(
        optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0)
    )
    @variable(m, y[j=1:1] >= 0)
    (m, y, cuts, nopt_cons, nfeas_cons) = SimpleBenders.benders_optimize!(m, y, data, optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0), f)
    (xref, yref, objref) = test_result()
    @test yref[1] ≈ JuMP.value(y[1])
    @test objref ≈ JuMP.objective_value(m)
end

@testset "Test feasibility cut" begin
    (xref, yref, objref) = test2_result()
    data = test2_data()
    f(v) = 3v[1] - 3v[2]
    m = Model(
        optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0)
    )
    @variable(m, 0<=y[j=1:2]<=3)
    (m, y, cuts, nopt_cons, nfeas_cons) = SimpleBenders.benders_optimize!(m, y, data, optimizer_with_attributes(Gurobi.Optimizer, "LogToConsole" => 0), f)
    
    @test yref[1] ≈ JuMP.value(y[1])
    @test yref[2] ≈ JuMP.value(y[2])
    @test objref ≈ JuMP.objective_value(m)
end

# SimpleBenders.jl

[![Build Status](https://github.com/matbesancon/SimpleBenders.jl/workflows/CI/badge.svg)](https://github.com/matbesancon/SimpleBenders.jl/actions)
[![Coverage](https://codecov.io/gh/matbesancon/SimpleBenders.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/matbesancon/SimpleBenders.jl)

A simple implementation of the Benders decomposition method with JuMP.

## Motivation

A good start with this package is the corresponding [blog post](https://matbesancon.github.io/post/2019-05-08-simple-benders/).
The implementation is far from smart since it uses the matrix form of the subproblem.

## Usage

`SimpleBenders.jl` assumes a problem of the form:

```
min f(y) + <c, x>
s.t.
G(x) ∈ S
A x + D y ≥ b
x ≥ 0
```

The projected-out sub-problem requires data stored in a dedicated structure:

```julia
struct SubProblemData
    b::Vector{Float64}
    D::Matrix{Float64}
    A::Matrix{Float64}
    c::Vector{Float64}
end
```

The different bindings correspond to the problem presented above.
Users can build their upper-level model `m` with the `y` variables, and define
the objective on the main variables: `f(y)`. The main function, `benders_optimize!`
can be called with:

```julia
benders_optimize!(m::Model, y::Vector{VariableRef}, sd::SubProblemData, sp_optimizer, f)
```

Which will solve the problem using Benders decomposition.

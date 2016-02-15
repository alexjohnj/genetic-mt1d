# Overview

genetic-mt1d is a 1D magnetotelluric (MT) inversion program that uses a real
valued genetic algorithm to find the best fitting inverse model for a data set
of apparent resistivity and phase responses. It was written as part of a project
for the GEOPH526 class at the University of Alberta. Although I can't make any
guarantees about the correctness of the results, I've had success using it to
reproduce the results of a paper by [Jones & Ferguson (2001)][electric-moho]
where they inverted 1D data to find the electric Moho.

[electric-moho]: http://www.nature.com/nature/journal/v409/n6818/abs/409331a0.html

# Installation

genetic-mt1d is written in [Julia][julia-lang] and tested against the latest
stable build (0.4, at the time of writing), so you'll need to install Julia to
get started. To install genetic-mt1d, run the following from a Julia REPL:

``` julia
Pkg.clone("https://github.com/alexjohnj/genetic-mt1d.git", "MT1DGeneticInversion")
```

Because I wasn't thinking when I named the repository, you have to specify the
package's name as well as its git address.

[julia-lang]: http://julialang.org

# Usage

Let's use a little example to explain things. Say you have some data in a CSV
file that consists of period, apparent resistivity, phase and their associated
errors. You want to invert a 3 layer model from this data. Load this data into
Julia using your preferred method (maybe `readdlm`?)  and wrangle it into a
matrix with 5 columns. The columns should be

1. Period (seconds)
2. Apparent resistivity (Ohm-metres)
3. Phase (degrees)
4. Apparent resistivity uncertainty (Ohm-metres)
5. Phase uncertainty (degrees)

This is your `data` matrix. This is what we're going to invert. Now, in a Julia
script, load genetic-mt1d:

``` julia
using MT1DGeneticInversion
```

The first thing to do is define some constraints for the inversion. We'll want
to place an upper and lower bound on resistivity and the depth to each
layer. Constraints are represented by the `LayerBC` type, a simple struct with a
`min` and `max` field. Each `LayerBC` instance represents the constraints for
one parameter of one layer so we'll need two 3 element arrays of
`LayerBC`s. Here, we'll set the same constraints for each layer:

``` julia
zBounds = [LayerBC(0, 0) LayerBC(0, 1000) LayerBC(0, 1000)] # Depth
rBounds = fill(LayerBC(1, 100), (3, 1)) # Resistivity
```

Note that for the depth constraints, the first layer is constrained to 0 m. This
is a requirement and will produce an error if it's constrained to anything
else. The depth to the other layers is constrained to 0--1000 m. Each layer's
resistivity is constrained to 1--100 Ωm.

Next we need to create an instance of the `Inversion` type using these
constraints. We'll create an `Inversion` with a population size of 100 for the
genetic algorithm.

``` julia
I = Inversion(data, 100, zBounds, rBounds)
```

This creates an `Inversion` instance initialised with a population of 100 random
models, generated within the constraints of `zBounds` and `rBounds`. The number
of layers in the models is derived from the number of constraints
provided. Using this constructor, the following parameters take their default
values:

- Number of elitist clones = 5% of the population size
- Probability of mutation = 0.01
- Probability of selection = 0.65
- Tournament size = 2
- SBX Distribution index (η) = 2

If you want to change these, you can either modify the fields of `I` directly
(see built in documentation) or use the more detailed constructor:

``` julia
I = Inversion(data, 100, zBounds, rBounds, nElitist, probMut, probSel, tournSize, η)
```

Once the `Inversion` instance is created, we can evolve it using:

``` julia
evolve!(I, 2000)
```

This evolves the population of `I` through 2000 generations, changing the state
of `I` in the process. To get the best model after evolution, use the first
element of `I`'s `pop` field:

```julia
bestModel = I.pop[1]
```

The `pop` field is an array of `Model` instances sorted from fittest to
weakest. `Model` instances have two important field. `Model.model` (confusing,
sorry) accesses the matrix representation of the model where each row is a layer
and the first column contains the depth to that layer while the second column
contains its resistivity. `Model.fitness` gives you the root mean squared (RMS)
misfit between the model's forward modelled response and the response in `data`.

That is a (somewhat) long winded (but complete!, kinda) overview of using
genetic-mt1d. Note that there is an example function in
`examples/examples.jl` that inverts some real world data to find the depth
to the electric Moho. You can probably use that as a starting point. Also, all
the code's documented and it integrates with Julia's help system nicely. Try
typing `?Inversion` at the REPL to get some more info on the `Inversion` type,
for example.

# How It Works

The inversion makes use of a real valued genetic algorithm. Models are
represented as Nx2 matrices (N = number of layers) and their apparent
resistivity and phase response is forward modelled using
[Wait (1954)][wait-recursion]. The genetic algorithm tries to minimise the RMS
misfit between the models and the provided data until it is between 0.5 and 1.5.

[wait-recursion]: http://library.seg.org/doi/abs/10.1190/1.1437994

Selection is implemented as a tournament selection with a customisable
tournament size and selection probability. Crossover uses a simulated binary
crossover (SBX) to combine two parents and produce two offsprings. In addition
to the children, a certain number of elitist clones are copied over to the next
generation. Mutation is applied per parameter in each model. For each parameter,
we generate a random number and if it's less than the mutation chance, we pick a
new value for the parameter within the layer's constraints. Because mutation is
applied per parameter, the effect of the mutation rate will depend on the number
of layers you're trying to invert in addition to the population size.


# Performance

Be warned, this isn't written with performance in mind. Nonetheless, the forward
modelling calculation is relatively simple and for smallish populations
(~100--200), performance is acceptable when running on a moderately powerful
laptop (30s--60s for the example script). The genetic algorithm doesn't make use
of multiple cores but one day I'd like to try and make it parallel. First I'll
need to read the docs on parallel computing in Julia though.

The forward modelling code is a straight up port of a MATLAB implementation I'd
written that was highly vectorised. I'm aware that Julia's performance is worse
with vectorised code so one day I'll explore devectorising the forward modelling
code. Some initial benchmarks forward modelling a 3 layer model across 100,000
frequencies showed a large memory allocation, so there might be some performance
to gain from devectorising.

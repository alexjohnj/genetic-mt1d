# genetic-mt1d

This is an inversion program for 1D magnetotelluric (MT) data that uses a
real-valued genetic algorithm. It was written for a report for GEOPH526 at the
University of Alberta. It works and has been tested on some synthetic data and
real world data included in the examples folder.

Currently, documentation for the code is a mixed bag. Some functions are
documented and some aren't. The `testGASnorcle` script in the examples folder
should give you an idea of how to get started with the program, but I'll add
some documentation in the future. At the moment, to change the mutation, elitism
and selection parameters for the algorithm, you need to edit `Chromosome.jl` and
`Population.jl`. In the future, I'm going to make these more easily
configurable.

## The Genetic Algorithm

For the inversion, models are defined in N&times;2 matrices where N is the
number layers. The first column contains the depth to the top of each layer and
the second column contains the layer's resistivity. The GA initialises a
population of random models within a set of provided constraints for depth and
resistivity. These constraints can be global or set on a per-layer basis. The
fitness function the GA is trying to minimise is the root mean squared (RMS)
misfit between the models and the data being inverted. The GA evolves the
population using four steps:

1. Selects parents for crossover using a tournament selection with two
   competitors.
2. Creates elitist clones of the best models and puts them into the next
   generation.
3. Breeds the selected parents using a simulated binary crossover (SBX) to fill
   the next generation to the target size.
4. Mutates the next generation by going through each layer in each model and
   picking a random number. If the number is less than the probability of
   mutation, a new depth or resistivity is picked at random within the initial
   constraints of the GA.

This continues until a certain number of generations is reached or until the
root mean squared misfit between the best model and the data being inverted is
in the range of 0.5 &leq; R &leq; 1.5.

By default, the following parameters are used in the different stages of the GA:

- Probability of selection (P<sub>s</sub>) = 0.6
- Probability of mutation (P<sub>m</sub>) = 0.01
- Elitism rate = 10% of the population
- Distribution parameter (η) = 2

As mentioned, to change them you need to modify the `Population.jl` and
`Chromosome.jl` source files.

## Usage

Running the program requires [Julia][julia-lang]. As previously mentioned, the
`testGASnorlce` script should be enough to get you up and running but basic
usage would go something like this:

[julia-lang]: http://julialang.org

```julia
include("./Population.jl")
data = readdlm("./mt.data")

N = 4 # Number of layers to invert
# Constraints for the depth to the top of each layer. Use one LayerBC type for
# each layer. The first layer must be fixed at 0m.
zBounds = [Layer(BC(0,0)) fill(LayerBC(1,100_000), (N-1,1))]
rBounds = fill(LayerBC(1,50_000), (N,1)) # Same but for resistivity

# Create a population of 100 random models. The models are generated within the
# provided boundaries and are sorted by fitness.
pop = Population(100, data, zBounds, rBounds)

# Evolve the population through 1000 generations, mutating the original pop
# type.
evolve!(pop, 1000)
writedlm("./genetic.model", pop.cs[1].model) # Save the best model to a text file.
```

The data you want to invert needs to be in a delimitered text file with the
following structure:

```
Period (s), Apparent Res (Ωm), Phase (°), Res Uncertainty, Phase Uncertainty
```

You pass it to the `Population` function as a matrix with the same order of
columns.

## The Example

The example script `testGASnorcle` inverts a 4 layer model from the
`SNO-96-106.data` file. This is data from the SNORCLE project, taken over the
Slave Craton in Canada in 1996. If you run the example script, you should
recover a nice four layer model after a few thousand generations. The recovered
model should match an inverted model from the same data set by
[Jones & Ferguson (2001)][jones-ferg]. Notably, there should be a drop in
resistivity around 36 km which corresponds to the Moho discontinuity.

[jones-ferg]: http://www.nature.com/nature/journal/v409/n6818/abs/409331a0.html

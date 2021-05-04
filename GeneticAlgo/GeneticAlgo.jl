## Code showing how genetic algorithms can be used to solve 
## an optimization problem. 
## Here we are finding x and y that maximize the function 
## (see the PDF file assignment_prob1.pdf)

module GeneticAlgo

import Random

##----------------------------------------------------------------------------
##  Types
##----------------------------------------------------------------------------
## While Chromosomes are modeled as strings,
## it is more convenient to work with Arrays. 
const Chromosome = Array{Int64, 1}

## An optimization problem will one or more variables that we have to 
## solve for.
## Varinfo saves 4 characteristics of each variable
struct Varinfo 
    name::Symbol
    minval::Float64
    maxval::Float64
    ## bitlen is the length of the chromosome encoding this specific 
    ## variable. Note this will be generally be different from the whole
    ## chromosome that encodes all variables.
    bitlen::Int64
end


## Parameters specific to the genetic algorithm 
mutable struct AlgParams
    ## popn_size - number of chromosomes in a generation of genetic algorithm
    popn_size::Int64
    trunc_rate::Float64
    mut_rate::Float64
end


## Data structure to store the results.
## Given the random elements in the genetic algorithm, one should 
## run different replications and look for optimum across these 
## different replications.
mutable struct Result
    rseed::Int64
    ## optimal_value of the objective function for the current rseed
    optimal_value::Float64
    ## best_chromosome is the chromosome that yields the optimal_value
    best_chromosome::Chromosome
    ## best_vals is the array of decoded values of actual variables 
    ## corresponding to the best chromosome
    best_vals::Array{Float64, 1}
    ## how many generations was the algorithm run for before stopping it.
    generations::Int64
end

##----------------------------------------------------------------------------
"""
You should not need to change this function.
Function to set the algorithm parameters.
"""
function set_alg_params(; popn_size = 8
                        , trunc_rate = 0.5
                        , mut_rate = 0.2)
    AlgParams( popn_size
             , trunc_rate
             , mut_rate
             )
end

"""
See the PDF file for the objective function
"""
function obj_func(x, y)
    ## The function should implement the objective function
    result = - (y+47)* sin(sqrt(abs((y+47) + (x/2)))) - x * sin(sqrt(abs(x-(y+47))))
    return result
end


"""
** There is no need to change this function **

Function to specify the optimization problem. 
Each of the two variables x and y should be represented using 14 bit chromosomes.
"""
function optproblem()
    Dict( :objfunc => obj_func
        , :vars => [ Varinfo(:x, -512, 512, 14)
                   , Varinfo(:y, -512, 512, 14)
                   ]
        )
end


"""
** There is no need to change this function **

Converts a chromosome with binary digits to an integer. 
"""
function decode_manual(chromosome::Chromosome)
    chrlen = length(chromosome)
    [chromosome[i] * 2^(chrlen - i) for i = 1:chrlen] |> sum
end



"""
** There is no need to change this function **

This method decodes chromosome of a single variable, say x. 

The other method decodes chromosome that contains information about 
all the variables in the problem (both x and y in this case).

When we call the method, Julia knows which method to apply based on the types
of the arguments we pass to it.
"""
function decode_chromosome(var_chromosome::Chromosome, 
                           vardetails::Varinfo)
    ## convert to integer
    intval = decode_manual(var_chromosome)

    ## decode to find the variable value
    varmin = vardetails.minval
    varmax = vardetails.maxval
    ## maxchrom is the maximum possible integer that can be represented
    ## using a chromosome that has the same length as the one passed in 
    ## as an argument.
    maxchrom = [2^(i-1) for i = 1:length(var_chromosome)] |> sum
    var_value = varmin + ((varmax - varmin)/maxchrom) * intval
end

"""
** There is no need to change this function **

Splits the full chromosome into constituent parts, where each part 
corresponds to one variable.

chromosome: full chromosome obtained by joining chromosomes of constituent variables.
chromosome_lengths are lengths of constituent chromosomes (E.g., [8, 8])
"""
function split_chromosome(chromosome::Chromosome, chromosome_lengths::Array{Int64, 1})
    cumul_lengths = cumsum(chromosome_lengths)
    indices = vcat(0, cumul_lengths)
    return [chromosome[indices[i]+1:indices[i+1]] for i = 1:(length(indices)-1)]
end




"""
** There is no need to change this function **

This is the second method that 
decodes the full chromosome containing information about 
all the variables. 
"""
function decode_chromosome(chromosome::Chromosome, problem)
    vars_info = problem[:vars]
    vars_bitlen = [v.bitlen for v in vars_info]
    var_chromosomes = split_chromosome(chromosome, vars_bitlen)
    decoded_vals = [decode_chromosome(var_chromosomes[i], vars_info[i]) 
                        for i = 1:length(vars_info)]
end



"""
This function should calculate the fitness of a chromosome. 
As you can see, you need to first decode it to the get the 
variable values and then pass them to the objective function.

Return value: Float64 (Value of the objective function at the given
                       chromosome)
"""
function calc_fitness(chromosome::Chromosome, problem)
    # Decode the full chromosome
    x,y  = decode_chromosome(chromosome, problem)
    
    ## Use the decoded values to calculate the value of the objective function
    return obj_func(x, y)
end


"""
This function returns a tuple of 4 items.
1. The best trunc_rate proportion of chromosomes. This becomes the 
   parent population for the next generation.
2. The fitness of the best chromosome in the population; this is the 
   current optimum.
3. The best chromosome.
4. The decoded variable values from the best_chromosome.
"""
function trunc_selection(alg_params, popn::Array{Chromosome, 1}, problem)
    popn_with_fitness = [(chrom, calc_fitness(chrom, problem)) for chrom in popn]
    sorted_popn_w_fitness = sort(popn_with_fitness, by = x -> x[2], rev = false)

    best_chromosome, fitness_best = sorted_popn_w_fitness[1]
    decoded_vals = decode_chromosome(best_chromosome, problem)

    num_to_keep = round(Int64, alg_params.trunc_rate * length(popn))
    chroms_to_keep = [sorted_popn_w_fitness[i][1] for i = 1:num_to_keep]
    
    (chroms_to_keep, fitness_best, best_chromosome, decoded_vals)
    # println(best_chromosome)        # test 9
    # println(fitness_best)
    # println(decoded_vals)

end


"""
Creates the necessary number of offsprings from the parent population.
For this application we are using the single-point crossover mechanism 
(see the notes on Genetic Algorithms).

Return value: Array of Chromosomes that are the offsprings. Size of the 
array must be the size of the population - size of the parent population.
"""
function create_offsprings(alg_params, parent_popn::Array{Chromosome, 1})
    num_offsprings = alg_params.popn_size - length(parent_popn)
    offsprings = Array{Chromosome, 1}(undef, 0)
    chromlen = length(parent_popn[1])
    while length(offsprings) < num_offsprings
        parents = rand(parent_popn, 2)
        crossover_point = rand(1:chromlen-1)
        offspring1 = vcat( parents[1][1:crossover_point]
                         , parents[2][crossover_point+1:chromlen]
        )
        offspring2 = vcat( parents[2][1:crossover_point]
                         , parents[1][crossover_point+1:chromlen]
        )
        
        if num_offsprings > length(offsprings)
            push!(offsprings, offspring1)
        end
        if num_offsprings > length(offsprings)
            push!(offsprings, offspring2)
        end
    end
    println(offsprings[1])    
    return offsprings
end


"""
Mutate the chromosomes. See the notes for why mutation is carried out. 

Do not mutate the best chromosome in the population. Others can possibly be mutated 
using the technique we discussed in the class.

Does not return anything. Changes the popn that is one of the inputs.
"""
function mutate_chromosomes!(alg_params, popn::Array{Chromosome, 1}, problem)
    popn_size = alg_params.popn_size
    mut_rate = alg_params.mut_rate
    chromlen = popn[1] |> length
    num_mutations = round(Int64, mut_rate * (length(popn) - 1) * chromlen)

    for i = 1:num_mutations
        ## Since we don't want to mutate the very best of the parent chromosomes 
        ## we ensure that rndrow cannot be 1.
        rndrow = rand(2:popn_size)
        rndcol = rand(1:chromlen)
        popn[rndrow][rndcol] = 1 - popn[rndrow][rndcol]
    end

end


"""
Combines the above function to run one generation of genetic algorithm.
    1. Takes the current population of chromosomes. 
    2. Chooses the best trunc_rate from these to make the parent population. 
    3. Records the fitness of the best chromosome (optimal_val).
    4. Records the best_chromosome and decoded_vals from the best_chromosome.
    4. Create offsprings from the parent population. 
    5. Parent population along with the offsprings creates a 
       new generation (called new_popn).

Returns: (new_popn, optimal_val, best_chromosome, decoded_vals)
"""
function run_generation(alg_params, popn::Array{Chromosome, 1}, problem, last_generation::Bool)
    parents, optimal_val, best_chromosome, decoded_vals = 
        trunc_selection(alg_params, popn, problem)
    offsprings = create_offsprings(alg_params, parents)
    new_popn = vcat(parents, offsprings)
    if !last_generation
        mutate_chromosomes!(alg_params, new_popn, problem)
    end
    return (new_popn, optimal_val, best_chromosome, decoded_vals)
    # println(optimal_val)
    # println(best_chromosome)
    # println(decoded_vals)
end


"""
This function uses the above defined functions to run simulation
for a given random seed until convergence occurs. 

Convergence criterion: The simulation compares the optimal value 
in the new generation with the current optimal value. If it fails 
to get a better optimal value for 25 generations, it gives up (that is,
it considers convergence to have occurred).
Anytime it gets a better optimal value with the new generation, it resets 
tries to zero. 

The `result` struct is returned when convergence has occurred.
"""
function run_simulation(rseed)
    ## Set the random seed (done)
    Random.seed!(rseed)
    ## Set the alg_params and the optimization problem (done)
    alg_params = set_alg_params()
    problem = optproblem()

    vars = problem[:vars]

    ## chromosome_len is the length of full chromosome that encodes 
    ## all variables. Find chromosome_len (done).
    chromsome_len = [v.bitlen for v in vars] |> sum

    ## initial population of chromosomes. 
    ## created randomly from all possible chromosome values (done)
    popn = [rand(0:1, chromsome_len) for i = 1:alg_params.popn_size]

    ## initialize some variables (done)
    gen = 0
    tries = 0
    overall_min = Inf

    ## initializing the result struct. This is where the best will be stored.
    result = Result(rseed, Inf, [], [], 0)

    while (tries < 25) 
        gen = gen + 1
        
        popn, optimal_val, best_chromosome, decoded_vals = run_generation(alg_params, popn, problem, false)
        
        if optimal_val < overall_min 
            overall_min = optimal_val 
            tries = 0
            result = Result(rseed, optimal_val, best_chromosome, decoded_vals, gen)
        else
            tries = tries + 1
        end
    end

    return result
end


"""
** There is no need to change this function **
"""
function main()
    overall_results = Array{Result, 1}(undef, 0)
    for rseed = 1:500
        if rseed % 100 == 0
            println("-"^80)
            println("rseed: $(rseed)")
        end
        result = run_simulation(rseed)
        push!(overall_results, result)
    end

    ## Once you have the result from each random seed, find the 
    ## overall optimum across these random seeds
    best_result = overall_results[1]
    for result in overall_results[2:end]
        if result.optimal_value < best_result.optimal_value
            best_result = result
        end
    end
    return overall_results, best_result
end

end ## end module

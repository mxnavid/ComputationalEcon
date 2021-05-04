## Don't change the module name or the filename.
## Do not import any additional packages.
module Schelling

import Random

import DataFrames
const Df = DataFrames

import StatsBase
const Sb = StatsBase


"""
You should not need to change this function.
""" 
function paramvalues(;threshold = 0.65, propred = 0.45, propblue = 0.45)
    Dict(
        :threshold => threshold, 
        :propred => propred, 
        :propblue => propblue
    )
end


"""
    You should not need to change this function.

    Randomly determine person type for a location, based on the :propred and
    :propblue values given in the parameter values. Used when initializing the grid.

    :R -> Red person
    :B -> Blue person
    :E -> Empty spot
"""
function rand_choice(pv)
    r = rand()
    if r <= pv[:propred]
        :R
    elseif r <= pv[:propred] + pv[:propblue]
        :B
    else
        :E 
    end
end

"""
    You should not need to change this function

    Create and return the initial grid - Two dimensional array

    Use symbol :R for red, :B for blue and :E for empty. 

    Also return an array of empty locations.  

    Each empty location is an array of the form [rownum, colnum]. This will be handy 
    when we want to know the set of possible locations where a person can move.
"""
function initial_grid(numrows, numcols, pv)
    grid = Array{Symbol}(undef, numrows, numcols)
    empty_locs = [] 
    for row in 1:numrows
        for col in 1:numcols
            persontype = rand_choice(pv)
            grid[row, col] = persontype
            if persontype == :E
                push!(empty_locs, [row, col])
            end
        end
    end
    return grid, empty_locs
end

"""
    You should not need to change this function.

    Get neighbor's location on the Schelling grid given 
    relative location (relloc) from a given cell.

    Function to determine the neighbor's row or col position given 
    person's [row, col] and the neighbor's location offset (example, 
    [1, 1] or [-1, 1] etc.). Cannot simply add offset to person location 
    because at the edges the offsets will result in invalid values. 
    For example, if personloc is [1,1], an offset of [-1, -1] will result in 
    neighbor location of [0,0], which is invalid (remember 
    Julia is 1-indexed). 

    This function changes the location to a valid one, if necessary.
    It assumes that the grid wraps around.
"""
function get_nbr_loc(ploc, rel_loc, nrows, ncols)
    local nbrrow, nbrcol   
    
    ## Let us find the neighbor's row
    if ploc[1] + rel_loc[1] <= 0
        nbrrow = nrows
    elseif ploc[1] + rel_loc[1] >= nrows
        nbrrow = 1
    else
        nbrrow = ploc[1] + rel_loc[1]
    end

    ## let us find the neighbor's column 
    if ploc[2] + rel_loc[2] <= 0
        nbrcol = ncols
    elseif ploc[2] + rel_loc[2] >= ncols
        nbrcol = 1
    else
        nbrcol = ploc[2] + rel_loc[2]
    end

    return [nbrrow, nbrcol]
end


"""
    You should not need to change this function

    Return an array consisting of the 8 neighboring locations on the grid.
"""
function get_neighbor_locs(grid, personloc)
    nrows, ncols = size(grid)

    all_rellocs = [[r,c] for r = -1:1 for c = -1:1 if [r,c] != [0,0]]

    [get_nbr_loc(personloc, rel_loc, nrows, ncols) for rel_loc in all_rellocs]
end

"""
    You should not need to change this function.

    Return the identity of the neighbors (removing empty spots :E)
"""
function get_neighbors(grid, personloc)
    nbrs = [grid[nbrloc...] for nbrloc in get_neighbor_locs(grid, personloc)]
    filter(nbr -> nbr != :E, nbrs)
end




"""
    You should not need to change this function.

    Calculate the perc/proportion of own type amongst the neigbors. In the
    Schelling model a person moves when proportion of own type falls below a
    certain threshold.
"""
function calc_perc_own(grid, personloc)
    owntypes = 0
    total = 0

    allnbrs = get_neighbors(grid, personloc)
    persontype = grid[personloc...]

    for nbr in allnbrs 
        if nbr == persontype
            owntypes += 1
        end
        total += 1
    end

    owntypes/total
end

## assumes that person moves to a random spot
"""
    You should not need to change this function.

    Move the person at `personloc` to a random location in 
    `empty_locs`. 

    Update the grid and empty_locs arguments. 
"""
function move_person!(grid, personloc, empty_locs)

    ## find the persontype
    pt = grid[personloc...]

    ## Find a random empty location: rndloc
    rndidx = rand(1:length(empty_locs))
    rndloc = empty_locs[rndidx]

    ## Change that empty location to persontype in the grid
    grid[rndloc...] = pt
    
    ## Change personloc to :E in the grid
    grid[personloc...] = :E

    ## update emptylocs
    empty_locs[rndidx] = personloc
end


"""
    You should not need to change this function.

    Calculate new state of the grid and new empty_locs from the existing 
    grid and empty_locs, by moving any person who is not satisfied at their 
    current location.
"""
function new_state!(grid, empty_locs, pv)
    ## for each square check if the person is satisfied in their current loc
    nrows, ncols = size(grid)
    for row in 1:nrows
        for col in 1:ncols
            if grid[row, col] != :E
                percown = calc_perc_own(grid, [row, col])
                if percown < 1 - pv[:threshold]
                    move_person!(grid, [row, col], empty_locs)
                end
            end
        end
    end
end

##-----------------------------------------------------------------------------
##    Diversity measures (DM2)
##-----------------------------------------------------------------------------
"""
    You should not need to change this function.

    This is the diversity measure that I refer to as DM2 in the my notes (see 
    the PDF file posted on the Schelling model).
"""
function calc_dm2(grid)
    nrows, ncols = size(grid)
    ## average proportion of the other type
    proportions = []
    for row = 1:nrows
        for col = 1:ncols
            if grid[row, col] != :E
                push!(proportions, 1 - calc_perc_own(grid, [row, col]))
            end
        end
    end
    avgprop = Sb.mean(proportions)
end

"""
    **Assignment Q1**

    Write function `prop_dissatisfied` that calculates the proportion of 
    all people that want to move given the schelling grid `grid` (a 
    two dimensional array of symbols) and paramvalues `pv`.

    Make sure that when you calculate the proportion of dissatisfied people you 
    don't count empty spots as people.

    Expected return value: a Float64 value
"""
function prop_dissatisfied(schgrid, pv)
    total_post = 0
    hit = 0
    rows, cols = size(schgrid)
    for curRow = 1:rows
        for curCol = 1:cols
            if schgrid[curRow, curCol] == :E
                continue
            else
                currentPerc = calc_perc_own(schgrid, [curRow, curCol])
                if currentPerc <= 1-pv[:threshold]
                    hit += 1
                end
            end
            total_post += 1
        end
    end
    dissatisfied = hit/total_post
    return dissatisfied
end

"""
    **Assignment Q2**

    Write function `simulate_movement` that takes 5 arguments:
        `pv`        - paramvalues dictionary
        `nrows`     - number of rows in the Schelling grid
        `ncols`     - number of colums in the Schelling grid
        `numsteps`  - integer specifying the number of times movement will be simulated through the entire grid.
        `rseed`     - random seed
        
    (1) Creates initial grid 
        (using the initial_grid function defined above and the given pv, nrows and ncols)
    (2) Calculates and stores 
        (a) `propdisbeg` the proportion of dissatisfied people in the initial grid (using the 
             function you completed above for this purpose).
        (b) diversity measure `dmbeg` using the initial grid (using `calc_dm2` function).
    (3) Creates new states successively building on updated states 
        using the new_state! function numsteps times (the function 
        `new_state!` should be called numsteps times).
    (4) After running `new_state!` numsteps times, uses the resulting grid to calculate and store 
        (a) `propdisend` the proportion of dissatisfied people in this final grid, and 
        (b) the diversity measure `dmend`.
    
    Expected return value:
    The output should be a dictionary with the following keys and values:
        key         |  value
        ---------------------
        :propdisbeg |  Value should be proportion of dissatisfied people in the initial grid.
        :dmbeg      |  Value should be the diversity measure dm2 calculated using the 
                    |  initial grid.
        :propdisend |  Value should be proportion of dissatisfied people at the end (after numsteps)
        :dmend      |  Value should be the diversity measure dm2 calculated using the 
                    |  grid state at the end of numsteps.

"""
function simulate_movement(pv, nrows, ncols, numsteps, rseed)
    Random.seed!(rseed)
    ## Rest of your code here...
    grid, elocs= initial_grid(nrows, ncols, pv)
    propdisbeg = prop_dissatisfied(grid, pv)
    dmbeg = calc_dm2(grid)

    for i in 1:numsteps
        new_state!(grid, elocs, pv)
    end

    propdisend = prop_dissatisfied(grid, pv)
    dmend = calc_dm2(grid)

    return Dict(
        :propdisbeg => propdisbeg, 
        :dmbeg => dmbeg, 
        :propdisend => propdisend,
        :dmend =>dmend
    )
end


"""
    **Assignment Q3**

    Write function run_reps that runs `numreps` replications of the schelling model, 
    and saves data from each replication in a dataframe.

    This function should do the following for `numreps` times. 
        (1) set the random seed
        (2) call `simulate_movement` function you completed above and get its output.
        (3) calculate the `perc_chng_dm` = 100 * ((dmend - dmbeg)/dmbeg)
        (4) store the output of `simulate_movement` along with random seed and `perc_chng_dm` in the dataframe.

    Expected return value:
    Dataframe with 5 columns:
        rseed
        propdisbeg
        dmbeg
        propdisend
        dmend
        perc_chng_dm

"""
function run_reps(numreps, pv, nrows, ncols, numsteps)
    df = Df.DataFrame(
        rseed = Int64[], 
        propdisbeg = Float64[], 
        dmbeg = Float64[], 
        propdisend = Float64[],
        dmend = Float64[],
        perc_chng_dm = Float64[]
    )

    for rseed = 1:numreps
        println("Running replication $(rseed)")
        movement = simulate_movement(pv, nrows, ncols, numsteps, rseed)
        propdisbeg_loc = movement[:propdisbeg]
        dmbeg_loc = movement[:dmbeg]
        propdisend_loc = movement[:propdisend]
        dmend_loc = movement[:dmend]     
        perc_chng_dm = 100 * ((dmend_loc - dmbeg_loc)/dmbeg_loc)
        push!(df,(rseed, propdisbeg_loc, dmbeg_loc, propdisend_loc, dmend_loc, perc_chng_dm))
    end
    return df
end

"""
    **Assignment Q4**

    Write the function `calcavgs` that takes the output of the `run_reps` function 
    and returns a dictionary containing the following keys and values

    Expected return value: Dictionary
        key                    |  value
        ---------------------------------------------------------------------------------
        :avg_propdisbeg        |  Mean of the `propdisbeg` column from the dataframe `df`
        :avg_dmbeg             |  Mean of the `dmbeg` column from the dataframe `df`
        :avg_propdisend        |  Mean of the `propdisend` column from the dataframe `df`
        :avg_dmend             |  Mean of the `dmend` column of the dataframe
        :avg_perc_chng_dm      |  Mean of the `perc_chng_dm` column of the dataframe
"""
function calcavgs(df)
    ## your code here...
    mean_propdisbeg = Sb.mean(df[!,:propdisbeg])
    mean_dmbeg = Sb.mean(df[!, :dmbeg])
    mean_propdisend = Sb.mean(df[!, :propdisend])
    mean_dmend = Sb.mean(df[!, :dmend])
    mean_perc_chng_dm = Sb.mean(df[!, :perc_chng_dm])

    return Dict(
        :avg_propdisbeg => mean_propdisbeg, 
        :avg_dmbeg => mean_dmbeg,
        :avg_propdisend =>avg_propdisend,
        :avg_dmend =>mean_dmend,
        :avg_perc_chng_dm =>mean_perc_chng_dm
    )
end






end ## end module
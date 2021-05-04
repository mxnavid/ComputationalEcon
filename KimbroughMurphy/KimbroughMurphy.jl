module KimbroughMurphy

import Random
import Distributions; const Dist = Distributions
import StatsBase; const Sb = StatsBase
import DataFrames; const Df = DataFrames
# import PrettyTables; const Pt = PrettyTables

"""
You don't need to change this function.

Market demand function: A - b * p.

"""
function paramvalues(; numfirms = 2, 
                       epochlens = Dict(i => 30 for i = 1:numfirms), 
                       initqs = Dict(i => 40 for i = 1:numfirms), 
                       deltas = Dict(i => 3.0 for i = 1:numfirms), 
                       epsilons = Dict(i => 0.7 for i = 1:numfirms), 
                       costs = Dict(i => 0.0 for i = 1:numfirms),  
                       A = 400, 
                       b = 2, 
                       T = 6600)
    

    Dict(
        :numfirms => numfirms,
        :epochlens => epochlens,
        :initqs => initqs,
        :deltas => deltas,
        :epsilons => epsilons,
        :costs => costs, 
        :A => A, :b => b, :T => T
    )
end


"""
You don't need to change this function.

Calculate the Cournot-Nash outcome

Return a dictionary consisting of:
    :firmqtys   - Each firms output in the Cournot Nash equilibrium.
    :Q          - Aggregate output in the Cournot Nash equilibrium.
    :P          - Market price in the Cournot Nash equilibrium.

"""
function cournot_nash(pv)
    A = pv[:A]
    b = pv[:b]
    costs = pv[:costs]
    aggcosts = costs |> values |> sum
    nfirms = length(costs)
    qtynum(i) = A + aggcosts - (nfirms + 1) * costs[i]
    qtyden = (nfirms + 1) * b
    firmqtys = Dict(i => qtynum(i)/qtyden for i = 1:nfirms)
    
    res = Dict{Symbol, Any}(:firmqtys => firmqtys)
    res[:Q] = firmqtys |> values |> sum

    res[:P] = A - b * res[:Q]
    return res
end

"""
Calculate firm's profit given:
    firmidx - one of {1, 2, ..., numfirms}.
    firmqty - qty produced by the firm.
    aggqty  - aggregate quantity in the market.
    pv      - paramvalues

Return value is of type Float64.
"""
function calc_firm_profit(firmidx, firmqty, aggqty, pv)
    ## calculate price using the inverse demand
    ## get firms cost of production form pv using firmidx
    ## calculate and return profit.
    local price::Float64
    local profit::Float64
    local cost::Float64
    price = pv[:A] - pv[:b]*aggqty
    cost = pv[:costs][firmidx]
    return firmqty * (price - cost)
end


"""
Draw the bidQuantity given current quantity (currq) and delta. 
bidQuantity is drawn from uniform distribution U[currq - delta, currq + delta].

Return value is of type Float64
"""
function bidquantity(currq, delta)
    ## Define the uniform distribution with the correct parameters.
    ## return a random value from this uniform distribution
    lo = currq - delta
    hi = currq + delta
    u = rand(Dist.Uniform(lo, hi))
    return u
end

"""
Calculate profits of all firms given the quantities dictionary and paramvalues.
quantity dictionary has firm index (1,2,3...) as keys and the 
quantity (Float64) as values. (firmidx (Int64) => quantity (Float64))

Return value: profits (Dictionary with firm index as keys and the profits (Float64) as values.)
              The keys should go from 1,2,...,numfirms
"""
function calc_profits(quantities, pv)
    ## calculate aggregate quantity using the quantities dictionary.
    ## use the aggregate quantity and quantities dict to calculate and 
    ## return profits of all firms (see the return value described above)
    aggQuantity = quantities |> values |> sum
    # print(aggQuantity)
    profit = Dict(i => calc_firm_profit(i, quantities[i], aggQuantity, pv) for i = 1:pv[:numfirms] )
    # for i = 1:pv[:numfirms]
    #     a = calc_firm_profit(i, quantities[i], aggQuantity, pv)    
    # end    
    return profit
    
end

"""
Calculate the amount by which the current quantity should be changed 
(increased or decreased) by firm `firmidx`, given the `firmreturnsdict` and 
`pv` (paramvalues).
`firmreturnsdict` is a dictionary in the following format
    :up     => Array of returns when firm produced a quantity the same or higher than its current currentQuantity
    :down   => Array of returns when firm produced a quantity lower than its current currentQuantity

Return value is of type Float64.
"""
function calc_currqchange(firmidx, firmreturnsdict, pv)
    ## calculate average profit from the firmreturnsdict[:up] array.
    avg_up = Sb.mean(firmreturnsdict[:up])
    ## calculate average profit from the firmreturnsdict[:down] array.
    avg_down = Sb.mean(firmreturnsdict[:down])
    ## if same or the former is bigger return the epsilon for this firm from pv
    ## if the latter is bigger return -1 * epsilon for this firm from pv
    if avg_up >= avg_down
        pv[:epsilons][firmidx]
    else
        pv[:epsilons][firmidx] * (-1)
    end
end


"""
Generate new currqtys dictionary and returnddict dictionary for all firms.
How? In the following manner:
Loop over each firmidx. For each firm check if it is a new epoch (epoch lengths 
can vary across firms). 

    If it is not a new epoch then keep the currqty the same for the firmidx and keep 
    the returnsdict of the firm the same.

    If it is a new epoch, calculate the new currqty using the current currqty and 
    calc_currqchange function. Also, create new firmreturndict with empty :up and :down 
    arrays.

To check whether it is a new epoch for a firm, calculate (`t` % firm's epoch length). 
It is a new epoch if this values is 1 and t != 1.

"""
function change_currqtys_reset_returnsdict(t, currqtys, returnsdict, pv)
    newcurrqtys = Dict() 
    newreturnsdict = Dict()
    epochlens = pv[:epochlens] 
    for firmidx in keys(currqtys)
        ## Is `t` the start of a new epoch? If yes, 
        if (t % epochlens[firmidx]) == 1 && t != 1
            # V1
            # calculate new currqty 
            newcurrqtys[firmidx] = currqtys[firmidx] + calc_currqchange(firmidx,returnsdict[firmidx], pv)
            # reset firm's retursdict
            newreturnsdict[firmidx] = Dict(:up => [], :down=>[])

            # # V2 
            # newCurrqty = currqtys[firmidx] + calc_currqchange(firmidx, returnsdict, pv)
            # push!(newcurrqtys, firmidx=>newCurrqty)
            # push!(newreturnsdict, firmidx => Dict(:up => [], :down =>[]))
        else
        ## if not,  
            ## add current currqty to newcurrqtys as it is
            newcurrqtys[firmidx] = currqtys[firmidx]
            ## add current firm's returnsdict to newreturnsdict as it is.
            newreturnsdict[firmidx] = returnsdict[firmidx]
        end
    end
    return newcurrqtys, newreturnsdict
    # println(newcurrqtys)      #test-15 => 52.3
end



"""
Runs one episode. 
Returns
    - new currqtys dict 
    - newreturnsdict
    - bidquantity

    For the first two return values, uses the function defined above.

    currqtys is a Dict of type firmidx => firm's current quantity
    returnsdict is Dict of type firmidx => (Dict of :up|:down => returns array)
    quantities is a Dict: firmidx => Array of quantities produced in the epoch
"""
function run_episode(t, currqtys, returnsdict, pv)
    newcurrqtys, newreturnsdict = change_currqtys_reset_returnsdict(t, currqtys, returnsdict, pv)

    bidqs = Dict()
    for firmidx in keys(newcurrqtys)
        ## populate bidqs with each firms bidQuantity for this episode
        bidqs[firmidx] = bidquantity(newcurrqtys[firmidx], pv[:deltas][firmidx])
    end

    ## calculate profits of all firms in this episode.
    # aggqty = 0
    # for i in keys(bidqs)
    #     aggqty = aggqty + bidqs[i]
    # end

    profit =  calc_profits(currqtys, pv)
    # for i in keys(bidqs)
    #     profit[i] = calc_firm_profit(i, bidqs[i], aggqty, pv)
    # end

    for (firmidx, newcurrqty) in newcurrqtys
        ## update the newreturnsdict dictionary appropriately with the 
        ## calculated profits.
        if bidqs[firmidx] >= newcurrqty
            push!(newreturnsdict[firmidx][:up], profit[firmidx])
        else
            push!(newreturnsdict[firmidx][:down], profit[firmidx])
        end
    end
    return (newcurrqtys, newreturnsdict, bidqs)
end

"""
Returns the quantities array for the market simulation (all epochs and all
episodes in one market).
"""
function simulate_mkt(pv)
    quantities = Dict(i => [] for i = 1:pv[:numfirms])
    ## initialize the currqtys dictionary
    currqtys = Dict()
    for i in 1:pv[:numfirms]
        currqtys = pv[:initqs]
    end
    ## initialize returnsdict dictionary for each firm
    newreturnsdict = Dict()
    for i in 1:pv[:numfirms]
        newreturnsdict[i] = Dict(:up => [], :down => [])
    end
    # println(newreturnsdicst)
    
    for t = 1:pv[:T]
        ## obtain new currqtys and new returnsdict and bidqs for episode t
        ## update quantities dictionary using the bidqs data.
        currqtys, newreturnsdict, bidqs = run_episode(t, currqtys, newreturnsdict, pv)
        for x in 1:pv[:numfirms]
            append!(quantities[x], bidqs[x])
        end
    end
    return quantities
end

"""
Return value looks like (given N firms)
    Dict(
        :q1     => Average of last obs_to_use observation from firm 1's quantities array,
        :q2     => Average of last obs_to_use observation from firm 2's quantities array,
        :q3     => Average of last obs_to_use observation from firm 3's quantities array,
        .
        .
        .
        :qN     => Average of last obs_to_use observation from firm N's quantities array, 
        :aggQ   => Sum of values of the keys :q1, :q2, :q3 etc. described above. This is the 
                   average aggregate quantity.
    )
"""
function calc_avg_and_aggavg_qty(quantities, obs_to_use)

    ## Create a dictionary that contains keys such as 
    ## q1, q2, q3 etc. where the number denotes the firmidx. 
    ## The value should be the average of firm's output in the *last* obs_to_use observations. 
    ## find the average aggregate quantity, :aggQ, using the aggregate of firm qty averages.
    ## add this to dictionary above with key :aggQ
    ## return the dictionary
    myDict = Dict(Symbol("q$(k)") => Sb.mean(v[(obs_to_use:end)]) for (k,v) in quantities)
    aggQuantity = values(myDict) |> sum
    myDict[:aggQ] = aggQuantity
    return myDict
end


"""
Should not have to change any code here...
"""
function run_reps(numreps, pv)
    ## This is the empty dataframe to store data
    df = Df.DataFrame( rseed = Int64[], aggQ = Float64[] )
    numfirms = pv[:numfirms]
    for firmidx in 1:numfirms
        df[:, Symbol("q$(firmidx)")] = Float64[]
    end

    ## number of observation used for calculating average
    obs_to_use = 1000
    for rseed in 1:numreps
        # println("Running rep $(rseed)")
        Random.seed!(rseed)
        quantities = simulate_mkt(pv) 
        d = Dict(:rseed => rseed)
        d = merge(d, calc_avg_and_aggavg_qty(quantities, obs_to_use)) 
        push!(df, d)
    end
    return df
end

"""
Should not have to change any code here...
This puts the output in a format that can be passed to pretty_table.
"""
function output(nasheqm, km_mkts, numfirms)
    vars = vcat([:aggQ], [Symbol("q$i") for i = 1:numfirms])
    nashoutput = vcat([nasheqm[:Q]], [nasheqm[:firmqtys][i] for i = 1:numfirms])
    kmoutput = [Sb.mean(km_mkts[:, col]) for col in vars]
    ([:Var, :Cournot, :KM], hcat(vars, nashoutput, kmoutput))
end


function main()
    pv = paramvalues(numfirms = 4)
    pv[:costs][1] = 10.0
    println("paramvalues: $(pv)")
    nasheq = cournot_nash(pv)
    km_mkts_df = run_reps(10, pv)
    # headers, data = output(nasheq, km_mkts, pv[:numfirms])
    # pt = Pt.pretty_table(data, headers)
    # println(pt)
    return nasheq, km_mkts_df
end

end ## end module
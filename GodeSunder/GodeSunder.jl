# Author: Professor Bharat Bhole
#         Navid Nafiuzzaman
# Description: Assignment 2 - GodeSunderAssignmentSp2021

module GodeSunder 

import Distributions; const Dist = Distributions
import StatsBase; const Sb = StatsBase
import DataFrames; const Df = DataFrames
import Random
import CSV
import Plots


## We need a Buyer type
mutable struct Buyer
    id :: Int64
    value :: Float64
    bid :: Float64
end


## We need a Seller type
mutable struct Seller
    id :: Int64
    cost :: Float64
    ask :: Float64
end


"""
**AssignmentQ1**
Write a function paramvalues that returns a dictionary of parameter values.

The function should take the following keyword arguments along with the shown
default values. If you don't remember what keyword arguments are, see the
docs/notes on how to define functions with keyword arguments.

    Name of the argument -> Description for your understanding (Default value)
    --------------------------------------------------------------------
    nb                   -> number of buyers (100)
    ns                   -> number of sellers (100)
    minval               -> minimum value in the buyer population (0)
    maxval               -> maximum value in the buyer population (200)
    mincost              -> minimum cost in the seller population (0)
    maxcost              -> maximum cost in the seller population (200)


Return value ->  Dictionary with keys 
                 :nb, 
                 :ns, 
                 :minval, 
                 :maxval, 
                 :mincost, 
                 :maxcost

Some examples of how the function should work:

    paramvalues()[:nb] = 100
    paramvalues()[:maxcost] = 200
    paramvalues(nb=500)[:nb] = 500
    paramvalues(nb=500)[:maxval] = 200
"""
function paramvalues(;
    nb = 100,
    ns = 100,
    minval = 0.0,
    maxval = 200.0,
    mincost = 0.0,
    maxcost = 200.0)
    dict = Dict(:nb => nb,
                :ns => ns,
                :minval => minval,
                :maxval => maxval,
                :mincost => mincost,
                :maxcost => maxcost
                )
    return dict
    ## function body
end


"""
You don't have to and should not need to modify this function.

Generate and return the population of buyers.

        Args: pv -> paramvalues dictionary
Return value: Array{Buyer, 1}
"""
function gen_buyers(pv)
    unifbuyers = Dist.Uniform(pv[:minval], pv[:maxval])
    rndvals = rand(unifbuyers, pv[:nb])
    [Buyer(id, rndvals[id], -Inf) for id = 1:pv[:nb]]
end

"""
You don't have to and should not need to modify this function.

Generate a population of sellers.

        Args: pv -> paramvalues dictionary
Return value: Array{Seller, 1}
"""
function gen_sellers(pv)
    unifsellers = Dist.Uniform(pv[:mincost], pv[:maxcost])
    rndvals = rand(unifsellers, pv[:ns])
    [Seller(id, rndvals[id], Inf) for id = 1:pv[:ns]]
end


"""
**AssignmentQ2**

Calculate the perfectly competitive equilibrium. 

This is going to be similar to the function we coded in the class. 
Except, you are required to make the following modifications.

1. Change the code to set the equilibrium price to the (average of
buyer valuation and seller cost -- both at the equilibrium quantity).

2. Add code to calculate average buyer value on traded units in the perfectly 
competitive market model.

3. Add code to calculate average seller cost on traded units in the perfectly 
competitive market model.

Return value - Dictionary with keys 
               :eqmqty, 
               :eqmprice,
               :eqmsurplus, 
               :eqmbv (for average equilibrium buyer value), and 
               :eqmsc (for average equilibrium seller cost)
"""
function calc_eqm(buyers::Array{Buyer, 1}, sellers::Array{Seller, 1}, pv)
    buyers_sorted = sort(buyers, by = b -> b.value, rev = true)
    sellers_sorted = sort(sellers, by = s -> s.cost)

    diff_array = [buyers_sorted[i].value - sellers_sorted[i].cost 
                    for i = 1:min(pv[:nb], pv[:ns])]

    diffarr_pos = filter(x -> x >= 0, diff_array)
    eqm_qty = length(diffarr_pos)
    eqm_price = buyers_sorted[eqm_qty].value
    
    
    # buyers_value = Any[]
    # for x in buyers
    #     append!(buyers_value, x.value)
    # end

    # sellers_cost = Any[]
    # for x in sellers
    #     append!(sellers_cost, x.cost)
    # end
    
    # eqm_bv = Sb.mean(buyers_value)
    # eqm_sc = Sb.mean(sellers_cost)    

    eqm_bv = sum([b.value for b in buyers_sorted[1:eqm_qty]])/eqm_qty
    eqm_sc = sum([s.cost for s in sellers_sorted[1:eqm_qty]])/eqm_qty
    eqm_ss = sum(diffarr_pos)

    Dict(:eqmqty => eqm_qty, 
         :eqmprice => eqm_price, 
         :eqmsurplus => eqm_ss,
         :eqmbv => eqm_bv,
         :eqmsc => eqm_sc)    
end

"""
You don't have to and should not need to modify this function.

Buyers randomly decide how much to offer (their bid) 
subject to the constraint that they don't pay over their valuation.
Given buyers array, update the array with their bid.
"""
function gen_bids!(buyers :: Array{Buyer, 1})
    for buyer in buyers
        bid_distn = Dist.Uniform(0, buyer.value)
        buyer.bid = rand(bid_distn)
    end
end

"""
You don't have to and should not need to modify this function.

Sellers randomly decide how much to ask  
subject to the constraint that they don't ask for less than their cost.
Given the sellers array, update the array with their asks.

- Since we deleted the constant MaxVal parameter, we need to give the maximum value 
  as an input. Hence, I have added the `maxval` as the second argument to this function. 

- Any place in the code where you call `gen_asks!` function you will have to ensure that 
  the function call takes this change into account.
"""
function gen_asks!(sellers :: Array{Seller, 1}, maxval)
    for seller in sellers
        ask_distn = Dist.Uniform(seller.cost, maxval)
        seller.ask = rand(ask_distn)
    end
end


"""
You don't have to and should not need to modify this function.

Every buyer purchases at most one unit. 
When a buyer engages in a trade, we remove that buyer from the list of 
potential buyers.
"""
function remove_buyer(buyers :: Array{Buyer, 1}, buyer_to_rem :: Buyer)
    filter(b -> b.id != buyer_to_rem.id, buyers)
end


"""
You don't have to and should not need to modify this function.

Every seller sells at most one unit. 
When a seller engages in a trade, we remove that seller from the list of 
potential sellers.
"""
function remove_seller(sellers :: Array{Seller, 1}, seller_to_rem :: Seller)
    filter!(s -> s.id != seller_to_rem.id, sellers)
end


"""
**AssignmentQ3**

Implement the double-auction market that runs for 0.5 seconds
For each trade, record the buyer value, seller cost and the price.
Return the dataframe consisting of these three columns. This will be 
the same as the function we wrote in the class, except the following 
changes:

- In the class we ran the market for 2 seconds. 
  Change the code to run the market for 0.5 seconds.
- In the class we set the price of a trade equal to buyer value. Now set the 
  price equal to the average of buyer bid and seller ask.
- Add two additional function arguments (I have already done this):
    (i) `numbuyers` (this is to replace the constant Nb)
   (ii) `maxval` (which is used in gen_asks function)

Return value: The function should return a dataframe with three columns: 
              buyervalue, 
              sellercost, and 
              price.
""" 
function market(buyers::Array{Buyer, 1}, 
                sellers::Array{Seller, 1}, 
                numbuyers::Int64, 
                maxval)
    trades = Df.DataFrame(buyervalue = Float64[], 
                          sellercost = Float64[], 
                          price = Float64[])
    
    buyerscopy = deepcopy(buyers)
    sellerscopy = deepcopy(sellers)
    currtime = time()
    
    while time() < currtime + 0.5
        gen_bids!(buyerscopy)
        gen_asks!(sellerscopy, maxval)
        i = 0
        while (i < 5 * numbuyers)               # 5 times the number of buyers before giving up
            i = i + 1
            buyer = Sb.sample(buyerscopy)
            seller = Sb.sample(sellerscopy)
            avg = Sb.mean([buyer.bid, seller.cost])         # average of bid and sale
            if (buyer.bid >= seller.ask)
                push!(trades, [buyer.value, seller.cost, avg])
                buyerscopy = remove_buyer(buyerscopy, buyer)
                sellerscopy = remove_seller(sellerscopy, seller)
                break
            end
        end        
    end
    return trades
end


"""
**AssignmentQ4 - Modified Double Auction Market**

Implement the double-auction market that runs for 0.5 seconds
For each trade, record the buyer value, seller cost and the price.
Return the dataframe consisting of these three columns. This will be 
the same as the above function, except for the following:

    Change in how trades are determined: 
    instead of choosing a buyer-seller pair randomly to trade, 
    choose the buyer with the highest bid and pair them with the seller who has the 
    lowest ask. 

    If the trade is feasible then record buyervalue, sellercost and the price 
    (= average of buyer bid and seller ask), and remove this buyer and seller. 

    Repeat the process starting with new bids and asks. This goes on for 0.5 seconds. 
    
    **Hint: For this market version you can do away with the inner while loop
            that we had in the code we worked through in the class**

Return value: The function should return a dataframe with three columns: 
              buyervalue, sellercost and price.
""" 
function modified_market(buyers::Array{Buyer, 1}, 
                         sellers::Array{Seller, 1}, 
                         maxval)
    trades = Df.DataFrame(buyervalue = Float64[], 
                         sellercost = Float64[], 
                         price = Float64[])
   
    buyerscopy = deepcopy(buyers)
    sellerscopy = deepcopy(sellers)
    currtime = time()
   
    while time() < currtime + 0.5
        gen_bids!(buyerscopy)
        gen_asks!(sellerscopy, maxval)
        buyers_sorted = sort(buyerscopy, by = b -> b.bid, rev = false)      # sort by bid
        sellers_sorted = sort(sellerscopy, by = s -> s.ask)                 # sort by ask
        
        buyer = buyers_sorted[length(buyers_sorted)]                        # pick the highest bidder
        seller = sellers_sorted[1]                                          # pick the lowest ask

        avg = Sb.mean([buyer.bid, seller.cost])         # average of bid and sale
        if (buyer.bid >= seller.ask)
            push!(trades, [buyer.value, seller.cost, avg])
            buyerscopy = remove_buyer(buyerscopy, buyer)
            sellerscopy = remove_seller(sellerscopy, seller)
        end
    end
    return trades
end



"""
**AssignmentQ5**

Calculate the total surplus, quantity and average price given 
`mktdf`: the market outcome as a dataframe with buyervalue, sellercost 
and avgprice as columns.

Modify the code to return two additional variables 
(in addition to the ones that were being returned earlier in the code we
worked through in the class):

    :avgbv (average of the values of the buyers who engage in trade)
    :avgsc (average of the costs of the seller who engage in trade)

Return value -> Dictionary with the following fields
                :quantity
                :avgprice
                :surplus
                :avgbv
                :avgsc

"""
function marketagg(mktdf)
    quantity = size(mktdf)[1]
    avgprice = Sb.mean(mktdf.price)
    surplus = sum(mktdf.buyervalue .- mktdf.sellercost)
    avgbv = Sb.mean(mktdf.buyervalue)
    avgsc = Sb.mean(mktdf.sellercost)
    
    return Dict(:quantity => quantity,
                :avgprice => avgprice,
                :surplus => surplus,
                :avgbv => avgbv,
                :avgsc => avgsc)
    
end


"""
**AssignmentQ6**
Run multiple replications of the simulation

A lot of the code will remain the same as that in the class. But you will need to 
make the following changes in order to accommodate the changes made above.

1. We are storing data on more variables in the eqmdf. Modify the `eqmdf` dataframe accordingly.
2. We are storing data on more variables in gsmktdf. Modify the `gsmktdf` dataframe accordingly.
3. The market function called here now needs `pv` (paramvalues dict) as an additional argument 
   (I have already added it to the list of arguments).
4. We have defined code for two types of market above. Add an argument called `mkt_type` that 
   allows us to specify which version of Gode Sunder market we want to run in replications. 
   Assume that `mkt_type` will be called with one of the following two values:
   :mkt_rand (where buyers and sellers who trade are randomly paired), or
   :mkt_mod (where highest bid buyer is paired with the lowest ask seller).
   Make sure your code calls the correct market function based on the value of `mkt_type`.

Return value -> Array containing two dataframes. 

                The first dataframe is the `eqmdf`. See the output of
                calc_eqm function to figure out the column names for this
                dataframe.

                The second dataframe is the gsmktdf. See the output of the marketagg function to 
                figure out the column names for this dataframe.

"""
function run_reps(numreps :: Int64, pv, mkt_type)
    ## container to save equilibrium outcomes.
    eqmdf = Df.DataFrame(rseed = Int64[],
                        eqmqty = Int64[],
                        eqmprice = Float64[],
                        eqmsurplus = Float64[],
                        eqmbv = Float64[],
                        eqmsc = Float64[])## complete the fields here

    ## container to save Gode-Sunder results.
    gsmktdf = Df.DataFrame(rseed = Int64[],
                            quantity = Int64[], 
                            avgprice= Float64[], 
                            surplus = Float64[],
                            avgbv = Float64[],
                            avgsc = Float64[])## complete the fields here

    ## parameters for usage in the function you will call for gs market
    nb = pv[:nb] 
    maxval = pv[:maxval]

    for rseed in 1:numreps
        ## Setting the random seed to rseed
        Random.seed!(rseed)
        println("Running replication $(rseed)...")

        ## Complete the rest of the code here...
        buyers = gen_buyers(pv)
        sellers = gen_sellers(pv)
        eqm = calc_eqm(buyers, sellers, pv)
        eqm[:rseed] = rseed
        push!(eqmdf, eqm)
        
        if mkt_type == :mkt_rand
            gs_trades = market(buyers, sellers, nb, maxval)
            gsmktagg = marketagg(gs_trades)
            gsmktagg[:rseed] = rseed
        else
            gs_trades = modified_market(buyers, sellers, maxval)
            gsmktagg = marketagg(gs_trades)
            gsmktagg[ :rseed] = rseed
        end
        push!(gsmktdf, gsmktagg)
    end
    return (eqmdf, gsmktdf)
    ## return eqm outcomes container and Gode-Sunder outcomes container
end


"""
**AssignmentQ7**
Write a function to calculate Efficiency of the GS market. 
Efficiency is defined as the ratio GS/ES, where: 
    GS -> Average social surplus realized in Gode-Sunder model across the replications. 
    ES -> Average social surplus realized in Perfectly competitive equilibrium model across
          the replications.

Return value -> A Float64 value (will be between 0 and 1).
"""
function calc_efficiency(eqmdf, gsmktdf)
    ES = Sb.mean(eqmdf.eqmsurplus)
    GS = Sb.mean(gsmktdf.surplus)
    val = GS/ES
    return val  
end

"""
**AssignmentQ8**

Write a function to calculate mean of percentage difference in average prices between the equilibrium model and the GS market.
Specifically, do the following:
    Step 1: Calculate the absolute difference between equilibrium market price (from the perfectly competitive model) and 
           the avgprice in the GS market for each replication.
    Step 2: Divide each difference in Step 1 by the equilibrium price of that replication (this will give the percentage difference).
    Step 3: Calculate the mean of values obtained in Step 2.

Return value -> A Float64 value (will be between 0 and 1).
"""
function calc_perc_price_diff(eqmdf, gsmktdf)
    x = 0
    eqmmarketprice_arr = Any[]      # eqm price for each iteration
    gsmavgprice_arr = Any[]         # gsm price for each iteration
    diff_arr = Any[]                # difference of the above two
    division = Any[]
    for x in eachrow(eqmdf)
        append!(eqmmarketprice_arr, x.eqmprice)
    end

    for x in eachrow(gsmktdf)
        append!(gsmavgprice_arr, gsmktdf.avgprice)
    end

    for s in 1:length(eqmmarketprice_arr)
        difference = abs(eqmmarketprice_arr[s]-gsmavgprice_arr[s])
        append!(diff_arr, difference)
    end

    for i in 1:length(diff_arr)
        divide = diff_arr[i]/eqmmarketprice_arr[i]
        append!(division, divide)
    end   
    
    return Sb.mean(division)
    
end


"""
You don't have to and should not need to modify this function.

This function runs each version of the market model for 100 replications, and 
compares both of them with the equilibrium outcome.
"""
function main()
    pv = paramvalues()
    results = Dict()

    for mkt_type in [:mkt_rand, :mkt_mod]
        println("Running simulations for $(mkt_type)...")
        eqmdf, gsmktdf = run_reps(100, pv, mkt_type)
        eff = calc_efficiency(eqmdf, gsmktdf)
        pricediff = calc_perc_price_diff(eqmdf, gsmktdf)
        results[mkt_type] = Dict(:eff => eff, :pricediff => pricediff)
        println("-"^80)
        println("Results for $(mkt_type)...")
        println("Efficiency: $(eff)")
        println("Perc price difference: $(pricediff)")
        println("-"^80)
    end
    return results
end

end
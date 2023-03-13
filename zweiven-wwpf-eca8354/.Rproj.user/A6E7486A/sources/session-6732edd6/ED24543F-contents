library(here)
library(data.table)

# Load the code to solve
#source(here("code", "SolveBulldozersProbabilistic.R"))
#source(here("code", "SolveBulldozersHeuristic.R"))

############
## Set up problem parameters
##############

# Core ecoregion information
ecoregion.data = fread(here("data", "alldata.csv"))  
# Mammal-based information on complementarity/species overlap across ecoregions
profile.data = fread(here("data", "mammal_based_profiles.csv"))
profile.abundance = profile.data[,filt.abundances]               
profiles = as.matrix(profile.data[,-ncol(profile.data), with=F])

# Baseline fixed parameters
kNumPeriods = 50
kBudgetPerPeriod = 1E9
kSARExponent = 0.2
kCostType = "GDPperm2_partial"
kDefSource = "hansen2018gross"
kDefType = "dozer_rate_2000_2018"
kFAType = "forest_area_2018"
kPFAType = "protected_forest_area_2018"
kCostEl = -6
# Endemic profiles with restoration value at 80% of pristine forest
# (see Newbold et al. 2015); 
kRestMultiplier = 0.8 

# Budget increment for heuristic solutions
kBudgetIncr = 100000

# Calculate area which *could* be reforested and its implications
# for the SAR multiplier
ecoregion.data[,forestable.area.init:=forest_area_2018+18*get(kDefType)]
ecoregion.data[,ae:=exp(log(plant_spcs) - kSARExponent*log(forestable.area.init))]

# Drop regions with missing or clearly bad cost data (i.e. cost should be positive), and
land.setup = ecoregion.data[!is.na(GDPperm2) & GDPperm2>0,
                            .(eco_code = eco_code, 
                              area.forestable.init = forestable.area.init,
                              area.forested.init = get(kFAType),
                              area.reserved.init = get(kPFAType),
                              species = plant_spcs,
                              ce = get(kCostType),
                              def.rate = get(kDefType),
                              ae=ae)] 
land.setup[,ce:=ce/((area.forestable.init-area.reserved.init)^(1/kCostEl))]

# Set up alternative species profiles under assumption of endemism
endemic.profiles = diag(nrow(land.setup))
endemic.abundances = land.setup$species

cost.el = -6


fixed.mc = is.na(cost.el)
land.setup[,forestable.area.available.init:= area.forestable.init-area.reserved.init] # Initially Available Area  
land.setup[,forested.area.available.init:= area.forested.init-area.reserved.init] # Initially Available Area
land.setup[,ze:=log(species/ae)/log(area.forestable.init)] # species area curve exponent




###################################################################################################################

library(Rfast) # for colprods

##############
# Constants
##############
# Tolerance for detecting when land has run out
EPSILON = 1E-3   
# Number of budget reallocation iterations to temporarily ignore an ecoregion
# if it is simply acting as a pass-through for budget transfers from
# one region to another. This occurs if a region is on a steep part
# of the ROI curve and can oscillate between highest and lowest ROI.
# This is just to speed things up; we will always return to the ecoregion.
NIGNOREITER = 10  

# First, define a helper function that characterizes the implications of a particular budget allocation
# Treat all problem parameters as (global) constants -- make this a function of the current candidate budget choice only.
# Gives back a list with the following elements given the candidate budget allocation:
# - cons.end.periods: the end period (race over period) per region
# - area.avail.consendper: the area still available during the race over period per region
# - end.areas.conserved: areas conserved at race end per region 
# - rois: the marginal ROI at the end of the race per region 
# Note the addition of profiles and profile.abundance. 
# - profiles is a matrix in which rows represent types of species, columns represent ecoregions, 
#   and entries indicate presence absence of that species type in that ecoregion.
# - profile.abundances is a vector with length equal to the number of rows in profile.abundance, indicating
#   how many species of that type we think there are.


FindRaceEndConditionsProbabilistic = function(cur.budget, 
                                              land.setup, 
                                              cum.devel.mat, 
                                              profiles, 
                                              profile.abundance, 
                                              cost.el=NA, 
                                              restoration.value=0, 
                                              roi.period) {
  num.regions = nrow(land.setup)
  num.periods = nrow(cum.devel.mat)
  anti.profiles = 1 - profiles
  fixed.mc = is.na(cost.el)
  
  
  
  with(land.setup,
       {
         # Find end period. This is the first period in which, given rates of deforestation and protection,
         # if both sides did not collide, the available land at the end of that period would be non-positive.
         
         # Cumulative land bought per current allocation at end of each period
         # - calculate cumulative spend as of period t, then compute land purchased using that vector
         cum.spend = apply(cur.budget, 1, cumsum)
         initial.forested.mat = matrix(area.forested.init, ncol=num.regions, nrow=num.periods, byrow=T)
         initial.forestable.mat = matrix(area.forestable.init, ncol=num.regions, nrow=num.periods, byrow=T)
         initial.protected.mat = matrix(area.reserved.init, ncol=num.regions, nrow=num.periods, byrow=T)
         if(class(ce)=="list") {
           ce.mat = matrix(unlist(ce), ncol=num.regions, nrow=num.periods, byrow=F)
         } else {
           ce.mat = matrix(ce, ncol=num.regions, nrow=num.periods, byrow=T)
         }
         
         # calculate hypothetical cumulative *new* area reserved by end of each period
         # if not constrained by bulldozers
         if(fixed.mc) {
           res.totals =  apply(cur.budget/t(ce.mat), 1, cumsum) 
         } else {
           if(cost.el!=-1) {
             ep.oneplusep = cost.el/(1+cost.el)
             if(class(ce)=="list") {
               # (potentially) time-varying cost multiplier per ecoregion
               res.totals = initial.forestable.mat - ((initial.forestable.mat-initial.protected.mat)^(1/ep.oneplusep) - apply(cur.budget/(t(ce.mat)*ep.oneplusep), 1, cumsum))^ep.oneplusep - initial.protected.mat
             } else {
               # time-invariant cost multiplier per ecoregion
               res.totals = initial.forestable.mat - ((initial.forestable.mat-initial.protected.mat)^(1/ep.oneplusep) - cum.spend/(ce.mat*ep.oneplusep))^ep.oneplusep - initial.protected.mat
             }
           } else {
             stop("cost elasticity of -1 not implemented")
           }
           # If we end up with NAs/NaNs, it's because there's more than enough spend to conserve everything, and the 
           # marginal and total cost functions aren't well defined once all forestable land has been protected.
           # we'll define both to be infinite beyond that point. Equivalently, we set res.totals to all available forestable land in those cases.
           res.totals[is.na(res.totals)] = (initial.forestable.mat-initial.protected.mat)[is.na(res.totals)]
           
           # We could end up with tiny negative numbers due to apparent numerical precision issues. If cum.spend is zero 
           # for a period and region, raising (initial.forestable.mat-initial.protected.mat) to 1/ep.oneplusep then to ep.oneplusep
           # doesn't result in an entry that exactly equals initial.forestable.mat-initial.protected.mat (but it should).
           # zero those out.
           res.totals[res.totals<0] = 0
         }
         
         # deal with R automatically retyping as numeric for single period case
         if(num.periods==1) {
           res.totals = matrix(res.totals, ncol=num.regions, nrow=num.periods)
         }
         rem.forested.land = forested.area.available.init - t(res.totals + cum.devel.mat) # How much available land is remaining at end of each period given current plan
         
         # Find out when the 'race' is over in each region. Earlier of:
         # a) when conservation meets development: pristine forest runs out by the end of that period
         # b) the end of the planning horizon.
         # we use epsilon as a threshold for detecting when land has run out since numerical precision
         # issues may result in failure to detect that event.
         meeting.periods = apply(rem.forested.land, 1, FUN=function(rem.vec) { min(which(rem.vec<=EPSILON)) })
         cons.end.periods = pmin(meeting.periods, num.periods)
         
         # Next, we want land conserved in each region when the 'race' is over in each region.
         # Because the budget may be more than enough to conserve all remaining land in the 
         # ending period for a particular region, we'll calculate land conserved at the start of
         # that period, calculate how much is conserved in the 'race over' period subject to 
         # feasibility constraints, and add both to the land initially reserved in a region.
         
         # Compute land newly reserved at start of conservation ending period 
         area.res.consendper.start = sapply(1:num.regions, function(r) {
           # find available land at start of period in which race ends.
           race.over.per = cons.end.periods[r]
           if(race.over.per > 1) {
             area.newly.res = res.totals[race.over.per-1, r]
           } else {
             area.newly.res = 0
           }
           return(area.reserved.init[r]+area.newly.res)
         })
         
         # Now find actual land conserved during final conservation period in each region.
         # To do so, we'll need to know constraints: how much is available in the final 
         # period during which conservation takes place
         area.avail.consendper = sapply(1:num.regions, function(r) {
           # find available land at start of period in which race ends.
           cons.end.per = cons.end.periods[r]
           if(cons.end.per > 1) {
             area.avail.consendper = rem.forested.land[r, cons.end.per-1]
           } else {
             area.avail.consendper = forested.area.available.init[r]
           }
           return(area.avail.consendper)
         })
         area.newly.conserved = 
           sapply(1:num.regions, function(r) {  # land added during race
             
             # find available land at start of period in which race ends.
             race.over.per = cons.end.periods[r]
             
             
             # Need to properly handle case in which race ends in first period
             if(race.over.per==1) {
               return(min(area.avail.consendper[r], res.totals[race.over.per,r]))
             } else {
               # In final period, the amount conserved is the minimum of that implied by the budget 
               # and the amount still available. Where the amount implied by the budget is larger, we
               # consider that to be wasted budget rather than an infeasible allocation.
               area.conserved.lastper = min(area.avail.consendper[r], res.totals[race.over.per,r]-res.totals[race.over.per-1,r])
               return(res.totals[race.over.per-1,r] + area.conserved.lastper)
             }
           })
         res.totals.at.consend = area.reserved.init + area.newly.conserved # add in land initially conserved
         # If restoration is allowed, we need to know conditions for when the restoration phase ends.
         # It will be the earlier of 
         # - when all forested land is either conserved or restored
         # - the end of the planning horizon
         # Restoration begins when conservation stops.
         # As specified, note restoration costs the same as protection, but is just less effective per 
         # the restoration.value multiplier. This means we can use protected land totals computed earlier, but just
         # need to know how much is conserved and how much is restored. 
         if(restoration.value > 0) {
           area.restored = res.totals-matrix(area.newly.conserved, nrow=num.periods, ncol=num.regions, byrow=T)
           area.restored[area.restored<0] = 0
           # for restoration, there is no development to consider, so we only need to constrain
           # restoration totals by the total amount of forested land initially available
           # (the sum of conservation and restoration can't exceed initially forested land)
           area.avail.for.rest = matrix(area.forestable.init - res.totals.at.consend,
                                        ncol=num.regions,
                                        nrow=num.periods, 
                                        byrow=T)
           
           # can't restore more than is available for restoration
           area.restored[area.restored>area.avail.for.rest] = area.avail.for.rest[area.restored>area.avail.for.rest]
           
           # Find out when restoration 'finishes' (if it does): when
           # area available for restoration runs out. 
           # We allow for epsilon error in detecting end periods as otherwise numerical precision 
           # issues may lead to a detection failure.
           rest.end.periods = sapply(1:num.regions, FUN=function(r) {
             rest.end = min(which((res.totals.at.consend[r]+area.restored[,r])>=(area.forestable.init[r]-EPSILON)))
           })
           # If we don't actually run out of land to restore in a region,
           # there will be no periods in which the area restored is at least as
           # large as the area available for restoration. The result of the
           # sapply call above will be Inf for that region.
           # In that case, the 'end' period for restoration is the end of the planning horizon,
           # and pmin will do the trick.
           # 
           rest.end.periods = pmin(rest.end.periods, num.periods)
           
           # Find out ending area restored per region
           rest.totals.at.restend = sapply(1:num.regions, FUN=function(r) {
             rest.end = area.restored[rest.end.periods[r],r]
           })
           # now we have a region x period tally of restoration totals (area.restored)
           # and a vector of periods during which restoration finishes (rest.end.periods)
           # and a vector of final restoration areas
         } else {
           rest.end.periods = cons.end.periods
           rest.totals.at.restend = rep(0, num.regions)
         }
         
         # Calculate ROI per region in the specified period, i.e., the marginal species benefit for a dollar spent
         # in each ecoregion in the specified period.
         # To do so in our probabilistic framework, we will compute
         # the expected increase in aggregate species protected across all ecoregions if we marginally
         # increase protection in each ecoregion. 
         # That increase is the sum over types of species in a given ecoregion of the increase in probability
         # that marginal land protection in the ecoregion will newly protect the species AND it is not protected
         # elsewhere.
         # 
         # We need to deal with two discontinuities in the ROI: 
         # 1. All forested land is conserved or developed
         # -- Here the restoration.value multiplier kicks in, which drops the marginal benefit
         # 2. All forestable land is conserved or reforested.
         # -- Here the marginal benefit drops to zero because the probability curve flattens at 1.
         # 
         # Converging on an optimal budget will work better if we treat these discontinuities specially,
         # so that we can get as close to the equimarginal principle defining the simple interior solution.
         # Otherwise we can get to oscillating budget movements where we shift from a region with, e.g., zero 
         # ROI, which then has the highest ROI, we shift budget back so it has zero, etc.
         # 
         # One option is to look at left and right ROIs (approaching from below and above), and then when
         # identifying source regions, we examine the left ROI (what would happen if we removed budget)
         # while we use the right ROI for targets (what if we added $).
         
         
         
         # Calculate probability of protection for each species type in each ecoregion
         # We assume that probability is constant across species in an ecoregion.
         # If a candidate budget allocation would result in more than the initially
         # foreseted area being protected, we cap the probability at one.
         if(restoration.value > 0) {
           p.protect = ((res.totals.at.consend + restoration.value*rest.totals.at.restend)/area.forestable.init)^ze
         } else {
           p.protect = (res.totals.at.consend/area.forestable.init)^ze
         }
         p.protect[p.protect>1]=1
         #print(p.protect)
         
         if(restoration.value > 0) {
           # marginal benefits depend upon whether the marginal unit of area affected by the marginal 
           # unit of budget in the specified period is conserved or restored
           marg.p.protect.right = marg.p.protect.left = ze*((res.totals.at.consend + restoration.value * rest.totals.at.restend)^(ze-1))/(area.forestable.init^ze)
           # Need to adjust marginal probability of protection if the marginal dollar would 
           # go toward restoration. That happens if, by the end of the roi period of interest,
           # remaining forested land that is not protected and not developed is at or below zero.
           roi.period.rest = (rem.forested.land[,roi.period]<0)
           marg.p.protect.left[roi.period.rest] = marg.p.protect.left[roi.period.rest]*restoration.value      
           roi.period.rest = (rem.forested.land[,roi.period]<=0)
           marg.p.protect.right[roi.period.rest] = marg.p.protect.right[roi.period.rest]*restoration.value
           
         } else {
           marg.p.protect.right=marg.p.protect.left = ze*(res.totals.at.consend^(ze-1))/(area.forestable.init^ze)
         }
         
         # And build up terms indicating probability of each species type being unprotected in each ecoregion.
         # this will be a matrix with a column per species type, a row per ecoregion,
         # and the entry representing the probability the species of that type is unprotected
         # in that ecoregion.
         marg.p.terms.left = profiles * matrix(marg.p.protect.left, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         marg.p.terms.right = profiles * matrix(marg.p.protect.right, ncol=nrow(land.setup), nrow=length(profile.abundance), byrow=T)
         # We could get an NaN if a profile entry is zero for a region and the marginal probability of protection is Infinite
         # because there's zero land protected (neither initially nor via the current budget). In those cases 0*Inf is NaN, but
         # we should treat the marginal probability terms as zeros since those species' probability of protection can't be 
         # influenced by protection in that region.
         marg.p.terms.right[is.nan(marg.p.terms.right)] = 0
         marg.p.terms.left[is.nan(marg.p.terms.left)] = 0
         #print(marg.p.terms.left[412,])
         unprot.terms = apply(profiles, 1, function(sp.type.pa.vec) { (1-p.protect*sp.type.pa.vec)})
         
         # now for each ecoregion x species type, we want the probability that
         # additional protected area would newly protect that species type. 
         # We then sum that across species in an ecoregion to get MB, and divide by MC to get ROI.
         
         # we can take full products of unprotected terms across all ecoregions once since that's expensive
         # then in the region-specific loop, we divide by the unprotected term from that focal region to remove its effect.
         fullprods = colprods(unprot.terms)
         print(fullprods)
         mbs.right = sapply(1:num.regions, function(i) {
           sum(profile.abundance*marg.p.terms.right[,i]*(fullprods/unprot.terms[i,]))
         })
         mbs.left = sapply(1:num.regions, function(i) {
           sum(profile.abundance*marg.p.terms.left[,i]*(fullprods/unprot.terms[i,]))
         })
         #print(mbs.left[412])
         #print(mbs.left[413])
         # We'll now compute marginal costs given allocations up to and including roi.period
         if(fixed.mc) {
           period.mcs.right = period.mcs.left = period.mcs = ce.mat[roi.period,]
           
         } else {
           # Need to calculate total reserved land (accounting for constraints of deforestation and forestable land)
           # in each area for comparison with forestable land, which will give us our MC
           roi.period.newly.res.totals = res.totals[roi.period,] 
           if(restoration.value==0) {
             roi.period.newly.res.totals[cons.end.periods==roi.period] =  area.newly.conserved[cons.end.periods==roi.period]
           } else {
             roi.period.newly.res.totals[rest.end.periods==roi.period] = (area.newly.conserved+rest.totals.at.restend)[rest.end.periods==roi.period]
           }
           # MC is a function of forestable land that is not yet protected via either conservation or reforestation
           if(class(ce)=="list") {
             # (potentially) time-varying cost multiplier per ecoregion
             period.mcs = ce.mat[roi.period,]*(area.forestable.init-(roi.period.newly.res.totals+area.reserved.init))^(1/cost.el)
           } else {
             # fixed cost multiplier per ecoregion
             period.mcs = ce*(area.forestable.init-(roi.period.newly.res.totals+area.reserved.init))^(1/cost.el)
           }
           # note again this is not well defined if all available land has been conserved via protection or reforestation, in which case 
           # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
           period.mcs.right = period.mcs.left = period.mcs
           period.mcs.right[res.totals[roi.period,]>=(area.forestable.init-area.reserved.init)] = Inf
           period.mcs.left[res.totals[roi.period,]>=(area.forestable.init-area.reserved.init)] = Inf
         }
         rois.right = mbs.right/period.mcs.right
         rois.left = mbs.left/period.mcs.left
         
         # In the edge case in which both ze and res.totals.at.race.over are zero, this will be NaN.
         # If no land is conserved AND conserving is useless (such that ze=0), then the ROI must 
         # be zero since spending can't increase forested area conserved
         rois.right[is.nan(rois.right)] = 0
         rois.left[is.nan(rois.left)] = 0
         
         
         
         # give back info about ending conditions, ROI, etc.
         if(restoration.value > 0) {
           race.end.periods = rest.end.periods
         } else {
           race.end.periods = cons.end.periods
         }
         
         
         # Calculate two things about the final period during which 
         # either conservation or restoration takes place. 
         # 1. area reserved or restored at start of that period
         area.resrest.raceendper.start = sapply(1:num.regions, function(r) {
           # find available land at start of period in which race ends.
           race.over.per = race.end.periods[r]
           if(race.over.per > 1) {
             area.newly.resrest = res.totals[race.over.per-1, r]
           } else {
             area.newly.resrest = 0
           }
           return(area.reserved.init[r]+area.newly.resrest)
         })
         
         # 2. Area available at the start of that period (end of preceding period).
         
         # note what we consider "available" in the race end period
         # will differ based on whether or not restoration is allowed.
         # - If yes, it's remaining initially forestable land that is not yet conserved.
         # - If no, it's remaining initially forested land that is neither conserved nor developed.
         area.avail.raceendper = sapply(1:num.regions, function(r) {
           # find available land at start of period in which race ends.
           race.end.per = race.end.periods[r]
           if(race.end.per > 1) {
             if(restoration.value > 0) {
               area.avail.raceendper = (forestable.area.available.init[r] - res.totals[race.end.per-1,r]) 
             } else {
               area.avail.raceendper = rem.forested.land[r, race.end.per-1]
             }
           } else {
             if(restoration.value > 0) {
               area.avail.raceendper = forestable.area.available.init[r]
             } else {
               area.avail.raceendper = forested.area.available.init[r]
             }
           }
           return(area.avail.raceendper)
         })
         
         
         return(list(cons.end.periods = cons.end.periods,
                     rest.end.periods = rest.end.periods,
                     race.end.periods = race.end.periods,
                     
                     area.res.consendper.start = area.res.consendper.start,
                     area.resrest.raceendper.start = area.resrest.raceendper.start,
                     
                     area.avail.consendper = area.avail.consendper,
                     area.avail.raceendper = area.avail.raceendper,
                     
                     end.areas.conserved = res.totals.at.consend,
                     end.areas.restored = rest.totals.at.restend,
                     rois.right = rois.right,
                     rois.left = rois.left)) 
         
       })
}




SolveBulldozersProbabilistic =  function(land.setup,
                                         num.periods,
                                         budget.per.period,
                                         profiles, 
                                         profile.abundance,
                                         cost.el=-999999999,
                                         restoration.value = 0,
                                         budget.conv.tol = 100,
                                         conv.tol = 1E-8,
                                         inner.conv.tol = 1E-8,
                                         realloc.tol=1E-12,
                                         budget.nudge = 0.1,
                                         max.fwd.pass=100,
                                         out.file=NULL) {
  
  
  fixed.mc = is.na(cost.el)
  land.setup[,forestable.area.available.init:= area.forestable.init-area.reserved.init] # Initially Available Area  
  land.setup[,forested.area.available.init:= area.forested.init-area.reserved.init] # Initially Available Area
  land.setup[,ze:=log(species/ae)/log(area.forestable.init)] # species area curve exponent
  
  # create some locals from the setup
  lapply(colnames(land.setup), FUN=function(cname) { 
    assign(x=cname, value=land.setup[,get(cname)], pos=1) 
  })
  
  # Compute a few convenience quantities based on problem setup
  num.regions = nrow(land.setup)  
  
  # Assemble matrix of cumulative development (if unchecked).
  # This will be derived from either a single fixed rate per ecoregion OR
  # a vector of rates per ecoregion, in which case def.rate is a list (not unlike cost)
  
  if(class(def.rate)=="list") {
    cum.devel.mat = apply(matrix(unlist(def.rate), nrow=num.regions, ncol=num.periods, byrow=TRUE),
                          1,
                          cumsum)
  } else {
    cum.devel.mat = apply(matrix(def.rate, nrow=num.regions, ncol=num.periods, byrow=FALSE), 
                          1, 
                          cumsum)
  }
  # convert to matrix for special case of 1 period problem since otherwise gets 
  # retyped as numeric            
  if(num.periods == 1) {
    cum.devel.mat = matrix(cum.devel.mat, nrow=num.periods, ncol=num.regions)
  }
  
  
  # Move budget to maximize objective, subject to feasibility
  # Guiding principle: 
  # - iteratively move budget from low to high marginal ROI 
  
  
  # Initial allocation equal across regions in each period
  cur.budget = c(1/num.regions) * matrix(budget.per.period, nrow = num.regions, ncol=num.periods, byrow=T) #Set initial proportion of budget allocation
  
  # Iniitalize other numerical tracking variables
  roi.diff = conv.tol*10
  num.pass.adjustments = 1
  last.race.over.time =ncol(cur.budget)
  fwd.pass = 1
  total.budget.chg = 100000
  

  last.source = last.target = -1
  recent.target = rep(F, num.regions)
  recent.source = rep(F, num.regions)
  t = 1
  # set up copy of budget for convergence checks. We'll look at the difference
  # across the entire time sweep.
  last.budget = cur.budget
  num.pass.adjustments = 0
  total.budget.chg = 0
  target.ignore.ctr = rep(0, nrow(cur.budget))

  print(t)
  # find race over times and ROIs for current allocation.
  # we'll actually do this twice: once considering a small increase in budget everywhere (for right ROIs), 
  # and once with a small decrease in budget everywhere (for left ROIs). This should help deal with 
  # numerical precision issues, especially around discontinuities.
  right.budget = left.budget = cur.budget
  right.budget[,t] = right.budget[,t] + budget.nudge
  left.budget[,t] = pmax(left.budget[,t] - budget.nudge,0)
  
  race.over.info = FindRaceEndConditionsProbabilistic(cur.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
  
  race.over.info.right = FindRaceEndConditionsProbabilistic(right.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)
  race.over.info.left = FindRaceEndConditionsProbabilistic(left.budget, land.setup, cum.devel.mat, profiles, profile.abundance, cost.el, restoration.value, t)

  
  print(race.over.info)
  return(race.over.info)
  }
###################################################################################################################


test = SolveBulldozersProbabilistic(land.setup,
                             num.periods = kNumPeriods,
                             budget.per.period = kBudgetPerPeriod, 
                             profiles = profiles,
                             profile.abundance = profile.abundance,
                             cost.el = kCostEl,
                             restoration.value = kRestMultiplier,
                             budget.conv.tol = 100,
                             conv.tol = 5E-8,
                             inner.conv.tol = 1E-9,
                             realloc.tol = 0.1, 
                             budget.nudge = 0.1,
                             out.file="results/final/main.Rdata")


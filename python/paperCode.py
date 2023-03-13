import pandas as pd
import numpy as np
import pprint as pp
import json
import sys
# Tolerance for detecting when land has run out
EPSILON = 1E-3
# Number of budget reallocation iterations to temporarily ignore an ecoregion
# if it is simply acting as a pass-through for budget transfers from
# one region to another. This occurs if a region is on a steep part
# of the ROI curve and can oscillate between highest and lowest ROI.
# This is just to speed things up; we will always return to the ecoregion
NIGNOREITER = 10    
#计算roi和周期相关参数的help function，主要是为从新分配budgets之后进行ROIs，race over和土地增加等数组的function，最后会return一个包含这些信息的dict
def FindRaceEndConditionsProbabilistic(cur_budget, 
                                    land_setup, 
                                    cum_devel_mat, 
                                    profiles, 
                                    profile_abundance, 
                                    cost_el = None, 
                                    restoration_value = 0,
                                    roi_period = 0):
    
    #roi_period为了R代码保持一致这里使用了R的index方式从1开始，但是python是从0开始，所以在使用roi_period是都使用roi_period - 1
    number_regions = len(land_setup)
    number_periods = len(cum_devel_mat)
    anti_profiles = 1 - profiles
    fixed_mc = cost_el is None

    #cum_spend和现当今budget相关，为每一年相加重复50次
    
    #cur_budget 是一个50*458矩阵 50行458列
    #Cumulative land bought per current allocation at end of each period
    cur_spend = np.cumsum(cur_budget.T,axis=0)
    #第一步生成50*458的矩阵三个和init相关初始矩阵
    #初始化dataframe里的数据，从数据扩展成矩阵
    initial_forested_mat = np.array(land_setup['area_forested_init'])#shape是50*458
    initial_forested_mat = initial_forested_mat.reshape((1, initial_forested_mat.size)).repeat(50, axis=0)

    initial_forestable_mat = np.array(land_setup['area_forestable_init'])#shape是50*458
    initial_forestable_mat = initial_forestable_mat.reshape((1, initial_forestable_mat.size)).repeat(50, axis=0)

    initial_protected_mat = np.array(land_setup['area_reserved_init'])#shape是50*458
    initial_protected_mat = initial_protected_mat.reshape((1, initial_protected_mat.size)).repeat(50, axis=0)

    #扩展ce到50*458维的矩阵
    #这里和R语言写的不太一样主要在于，没有对于'ce'的数据类型进行判断
    ce_mat = np.array(land_setup['ce'])#shape是50*458
    ce_mat = ce_mat.reshape((1,ce_mat.size)).repeat(50,axis=0)
    # calculate hypothetical cumulative *new* area reserved by end of each period if not constrained by bulldozers
    #每个时间节点后reserved面积总量
    if fixed_mc :
        res_totals = np.cumsum(cur_budget.T/ce_mat,axis=0)#shape是50*458
    else:
        if cost_el != -1:
            ep_oneplusep = cost_el/(1 + cost_el)
            res_totals = initial_forestable_mat - np.power(np.power((initial_forestable_mat - initial_protected_mat), 1/ep_oneplusep ) - cur_spend/(ce_mat * ep_oneplusep),ep_oneplusep )  - initial_protected_mat#shape是50*458
        else:
            raise Exception("cost elasticity of -1 not implemented")
        res_totals[np.isnan(res_totals)] = (initial_forestable_mat - initial_protected_mat)[np.isnan(res_totals)]
        res_totals[res_totals < 0] = 0
    
    if number_periods == 1:
        res_totals = np.array(res_totals).reshape(number_regions, number_periods)
    #剩余森林土地的计算
    rem_forested_land_trans = np.array(land_setup['forested_area_available_init']) - (res_totals + cum_devel_mat)#shape为458 * 50
    rem_forested_land = rem_forested_land_trans.T
    # Find out when the 'race' is over in each region. Earlier of:
    # a) when conservation meets development: pristine forest runs out by the end of that period
    # b) the end of the planning horizon.
    # we use epsilon as a threshold for detecting when land has run out since numerical precision
    # issues may result in failure to detect that event.
    #当剩余土地里袁术数值小于等于EPSILON的index都取出来，再取最小之，如果没有的话就去50（因为50年是最大年限）
    cons_end_periods = []
    for r in range(0,number_regions):
        #取rem_forested_land中满足条件的元素的index
        index = np.where(rem_forested_land[r] <= EPSILON)[0]
        if len(index) > 0:
            cons_end_periods.append(min(index))
        else:
            cons_end_periods.append(number_periods-1)
    cons_end_periods = np.array(cons_end_periods)
    
    

    # Next, we want land conserved in each region when the 'race' is over in each region.
    # Because the budget may be more than enough to conserve all remaining land in the 
    # ending period for a particular region, we'll calculate land conserved at the start of
    # that period, calculate how much is conserved in the 'race over' period subject to 
    # feasibility constraints, and add both to the land initially reserved in a region.
    
    # Compute land newly reserved at start of conservation ending period
    #计算newly_reserved
    area_res_consendper_start = []
    for r,rac_over_per in enumerate(cons_end_periods):
        if rac_over_per > 0:
            area_newly_res = res_totals[int(rac_over_per) - 1,r]
        else:
            area_newly_res = 0
        area_res_consendper_start.append(np.array(land_setup['area_reserved_init'])[r]+area_newly_res)
    area_res_consendper_start = np.array(area_res_consendper_start)

    
    # Now find actual land conserved during final conservation period in each region.
    # To do so, we'll need to know constraints: how much is available in the final 
    # period during which conservation takes place
    
    area_avail_consendper = []
    for r in range(0,number_regions):
        cons_end_per = cons_end_periods[r]
        # find available land at start of period in which race ends.
        #判断race over是什么时间
        if cons_end_per > 0:
            area_avail_consendper.append(rem_forested_land[r,int(cons_end_per) - 1])
        else:
            area_avail_consendper.append(np.array(land_setup['forested_area_available_init'])[r])
    area_avail_consendper = np.array(area_avail_consendper)

    area_newly_conserved = []
    for r in range(0,number_regions):
        rac_over_per = cons_end_periods[r]
        # Need to properly handle case in which race ends in first period
        #当race over发生在peroid
        if rac_over_per == 0:
            area_newly_conserved.append(min(area_avail_consendper[r],res_totals[int(rac_over_per),r]))
        else:
            # In final period, the amount conserved is the minimum of that implied by the budget 
            # and the amount still available. Where the amount implied by the budget is larger, we
            # consider that to be wasted budget rather than an infeasible allocation.
            area_conserved_lastper = min(area_avail_consendper[r],res_totals[int(rac_over_per),r] - res_totals[int(rac_over_per) - 1,r])

            area_newly_conserved.append(res_totals[int(rac_over_per) - 1,r] + area_conserved_lastper)
        
    area_newly_conserved = np.array(area_newly_conserved)
    # add in land initially conserved
    #把新conserve的土地加在原来的初始值上
    res_totals_at_consend = np.array(land_setup['area_reserved_init']) + area_newly_conserved

    
    # If restoration is allowed, we need to know conditions for when the restoration phase ends.
    # It will be the earlier of 
    # - when all forested land is either conserved or restored
    # - the end of the planning horizon
    # Restoration begins when conservation stops.
    # As specified, note restoration costs the same as protection, but is just less effective per 
    # the restoration.value multiplier. This means we can use protected land totals computed earlier, but just
    # need to know how much is conserved and how much is restored. 
    #这一部分是为了计算restoration end的时间节点
    if(restoration_value > 0):
        area_restored = res_totals - area_newly_conserved
        area_restored[area_restored < 0] = 0
        area_avail_for_rest = np.array(land_setup['area_forestable_init']) - res_totals_at_consend
        area_avail_for_rest = area_avail_for_rest.reshape((1, area_avail_for_rest.size)).repeat(50, axis=0)
        # can't restore more than is available for restoration
        area_restored[area_restored > area_avail_for_rest] = area_avail_for_rest[area_restored > area_avail_for_rest]
        
        #
        rest_end_periods = []
        for r in range(0,number_regions):
            # mask = (res_totals_at_consend[r] + area_restored[:,r]) >= (np.array(land_setup['area_forestable_init'])[r] - EPSILON)
            #找出符合条件的index
            index = np.where((res_totals_at_consend[r] + area_restored[:,r]) >= (np.array(land_setup['area_forestable_init'])[r] - EPSILON))
            if len(index[0]) > 0:
                rest_end_periods.append(min(index[0]))
            else:
                rest_end_periods.append(number_periods - 1)
        rest_end_periods = np.array(rest_end_periods)

        rest_totals_at_restend = []
        for r in range(0,number_regions):
            rest_end = area_restored[int(rest_end_periods[r]),r]
            rest_totals_at_restend.append(rest_end)
        rest_totals_at_restend = np.array(rest_totals_at_restend)
    else:
        rest_end_periods = cons_end_periods.copy()
        rest_totals_at_restend = np.zeros(number_regions)
        
    
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
    
    #计算ROI的方法，使用到生物的数量'ze'是论文中的z值 z=0.2
    if restoration_value > 0:
        p_protect = np.power(((res_totals_at_consend + restoration_value * rest_totals_at_restend)/np.array(land_setup['area_forestable_init'])),np.array(land_setup['ze']))
    else:
        p_protect = (res_totals_at_consend/np.array(land_setup['area_forestable_init']))**np.array(land_setup['ze'])
    p_protect[p_protect>1] = 1

    if restoration_value > 0:
        # marginal benefits depend upon whether the marginal unit of area affected by the marginal 
        # unit of budget in the specified period is conserved or restored
        # 根据公式计算marg_p_protect，左右没有区别的
        marg_p_protect_right = np.array(land_setup['ze'])*np.power((res_totals_at_consend 
                                + restoration_value * rest_totals_at_restend),(np.array(land_setup['ze']) - 1))/np.power(np.array(land_setup["area_forestable_init"]),np.array(land_setup['ze']))
        marg_p_project_left = np.array(land_setup['ze'])*np.power((res_totals_at_consend 
                                + restoration_value * rest_totals_at_restend),(np.array(land_setup['ze']) - 1))/np.power(np.array(land_setup["area_forestable_init"]),np.array(land_setup['ze']))
        #调整去向restoration地区的roi，当这个peroid结束时remaing forested——land 没有被protect或者developed应该是roi小于或者等于0的
        roi_period_rest = rem_forested_land[:,int(roi_period) - 1] < 0
        marg_p_project_left[roi_period_rest] = marg_p_project_left[roi_period_rest] * restoration_value
        roi_period_rest = rem_forested_land[:,int(roi_period) -  1] <= 0
        marg_p_protect_right[roi_period_rest] = marg_p_protect_right[roi_period_rest] * restoration_value
    else:
        marg_p_protect_right = land_setup['ze']*(res_totals_at_consend**(land_setup['ze'] - 1))/(np.array(land_setup['area_forestable_init'])**land_setup['ze'])
        marg_p_project_left = land_setup['ze']*(res_totals_at_consend**(land_setup['ze'] - 1))/(np.array(land_setup['area_forestable_init'])**land_setup['ze'])
    # And build up terms indicating probability of each species type being unprotected in each ecoregion.
    # this will be a matrix with a column per species type, a row per ecoregion,
    # and the entry representing the probability the species of that type is unprotected
    # in that ecoregion.
    #扩展各种生物的分布与其他的矩阵size匹配，方便后续计算
    marg_p_terms_left = np.array(profiles) *  marg_p_project_left.reshape((1, marg_p_project_left.size)).repeat(len(profiles), axis=0)
    marg_p_terms_right = np.array(profiles) * marg_p_protect_right.reshape((1, marg_p_protect_right.size)).repeat(len(profiles), axis=0)
    #计算时候可能会出现NAN值，把这些地区的值设成0
    marg_p_terms_left[np.isnan(marg_p_terms_left)] = 0
    marg_p_terms_right[np.isnan(marg_p_terms_right)] = 0

    #unprotected terms的计算
    unprot_terms = []
    for row in np.array(profiles):
        specie_p = 1 - p_protect * row.T
        unprot_terms.append(specie_p)
    unprot_terms = np.array(unprot_terms).T
    fullprods = np.cumprod(unprot_terms, axis = 0)[-1]

    # we can take full products of unprotected terms across all ecoregions once since that's expensive
    # then in the region-specific loop, we divide by the unprotected term from that focal region to remove its effect.
    # profile_abundance 是表示各个物种的总数
    mbs_right = []
    for i in range(0,number_regions):
        s_right = np.sum(np.array(profile_abundance)*marg_p_terms_right[:,i]*(fullprods/unprot_terms[i]))
        mbs_right.append(s_right)
    
    mbs_left = []
    for i in range(0,number_regions):
        s_left = np.sum(np.array(profile_abundance)*marg_p_terms_left[:,i]*(fullprods/unprot_terms[i]))
        mbs_left.append(s_left)
    
    mbs_left = np.array(mbs_left)
    mbs_right = np.array(mbs_right)    
    
    
    
    
    if fixed_mc:
        period_mcs_right = period_mcs = ce_mat[int(roi_period) - 1,]
        period_mcs_left = period_mcs = ce_mat[int(roi_period) - 1,]
    else:
        # Need to calculate total reserved land (accounting for constraints of deforestation and forestable land)
        # in each area for comparison with forestable land, which will give us our MC
        roi_period_newly_res_totals = res_totals[int(roi_period) - 1]
        if restoration_value == 0:
            roi_period_newly_res_totals[cons_end_periods == (int(roi_period) - 1)] = area_newly_conserved[cons_end_periods == int(roi_period) - 1]
        else:
            roi_period_newly_res_totals[rest_end_periods == (int(roi_period) - 1)] = (area_newly_conserved+rest_totals_at_restend)[rest_end_periods == (int(roi_period) - 1)]
        
        # MC is a function of forestable land that is not yet protected via either conservation or reforestation
        period_mcs = ce_mat[int(roi_period) - 1] * (np.array(land_setup['area_forestable_init']) - (roi_period_newly_res_totals + np.array(land_setup['area_reserved_init'])))**(1/cost_el)
        
        period_mcs_right = period_mcs.copy()
        period_mcs_left = period_mcs.copy()
        # note again this is not well defined if all available land has been conserved via protection or reforestation, in which case 
        # MC should be infinite and no additional land can be conserved (we correct res.totals earlier).
        period_mcs_right[res_totals[int(roi_period) - 1] >= (np.array(land_setup['area_forestable_init']) - np.array(land_setup['area_reserved_init']))] = np.Inf
        period_mcs_left[res_totals[int(roi_period) - 1] >= (np.array(land_setup['area_forestable_init']) - np.array(land_setup['area_reserved_init']))] = np.Inf
    rois_right = mbs_right/period_mcs_right
    rois_left = mbs_left/period_mcs_left
    #处理NAN值，当ze 和 res.totals.at.race.over 都是0是就会出现NAN
    rois_right[np.isnan(rois_right)] = 0
    rois_left[np.isnan(rois_left)] = 0
    
    
    
    # give back info about ending conditions, ROI, etc.
    if restoration_value > 0:
        race_end_periods = rest_end_periods.copy()
    else:
        race_end_periods = cons_end_periods.copy()
    
    
    # Calculate two things about the final period during which 
    # either conservation or restoration takes place. 
    # 1. area reserved or restored at start of that period
    area_resrest_raceendper_start = []
    for r in range(0,number_regions):
        race_over_per = race_end_periods[r]
        if rac_over_per > 0:
            area_newly_resrest = res_totals[int(race_over_per) - 1, r]
        else:
            area_newly_resrest = 0
        area_resrest_raceendper_start.append(np.array(land_setup['area_reserved_init'])[r] + area_newly_resrest)
    area_resrest_raceendper_start = np.array(area_resrest_raceendper_start)
    
    # 2. Area available at the start of that period (end of preceding period).
    
    # note what we consider "available" in the race end period
    # will differ based on whether or not restoration is allowed.
    # - If yes, it's remaining initially forestable land that is not yet conserved.
    # - If no, it's remaining initially forested land that is neither conserved nor developed.
    area_avail_raceendper = []
    for r in range(0,number_regions):
        race_end_per = race_end_periods[r]
        if race_end_per > 0:
            if restoration_value > 0:
                area_avail_raceendper.append(np.array(land_setup['forestable_area_available_init'])[r] - res_totals[int(race_end_per)-1,r])
            else:
                area_avail_raceendper.append(rem_forested_land[r,int(race_end_per) - 1])
        else:
            if restoration_value > 0:
                area_avail_raceendper.append(np.array(land_setup['forestable_area_available_init'])[r])
            else:
                area_avail_raceendper.append(np.array(land_setup['forested_area_available_init'])[r])    
    area_avail_raceendper = np.array(area_avail_raceendper)
    
    result ={
        'cons_end_periods':cons_end_periods,
        'rest_end_periods':rest_end_periods,
        'race_end_periods':race_end_periods,
        
        'area_res_consendper_start':area_res_consendper_start,
        'area_resrest_raceendper_start':area_resrest_raceendper_start,
        
        'area_avail_consendper':area_avail_consendper,
        'area_avail_raceendper':area_avail_raceendper,
        
        'end_areas_conserved':res_totals_at_consend,
        'end_areas_restored':rest_totals_at_restend,
        'rois_right':rois_right,
        'rois_left':rois_left
        }
        
    return result

# Function to actually solve the race as a function
# of the specified setup. Most arguments are vectors, with one entry
# per region
# Arguments:
# land.setup: data.table with information on ecoregions (one row per ecoregion)在main函数里初始化的，为从dataset中取出来的相关信息
# num.periods: number of periods in optimization (50 in base case in paper)
# budget.per.period: max budget that can be spent per period
# profiles: matrix (rows: profiles, columns: ecoregions) of presence/absence data for species profiles
# profile.abundance: vector of abundances of different species profiles
# restoration.value: multiplier indicating relative value of restored forest vs. intact forest for species protection
# cost.el: cost elasticity of land (as a proxy for more general protection costs)
# budget.conv.tol: tolerance for when budget shifts are sufficiently small to consider the algorithm to have converged
# conv.tol: tolerance for when ROI differences are sufficiently small to consider the algorithm to have converged
# inner.conv.tol: tolerance for when ROI differences are sufficiently small for inner search for budget reallocation to be complete (avoid overshooting on individual reallocations)
# realloc.tol: tolerance for when budget shifts are sufficiently small for inner search for budget reallocation to be complete (avoid overshooting on individual reallocations)
# budget.nudge: small amount by which we shift budget either direction when computing ROI for a budget increase or decrease to deal with potential discontinuities in ROI and numerical imprecision
# max.fwd.pass: limit on forward passes through the time horizon before exiting out
# out.file: name of filename to save results to. If null, results simply returned to calling code.
#主要的计算function，结束后返回一个result的字典
def SolveBulldozersProbabilistic(land_setup,
                                        num_periods,
                                        budget_per_period,
                                        profiles, 
                                        profile_abundance,
                                        cost_el=-999999999,
                                        restoration_value = 0,
                                        budget_conv_tol = 100,
                                        conv_tol = 1E-8,
                                        inner_conv_tol = 1E-8,
                                        realloc_tol=1E-12,
                                        budget_nudge = 0.1,
                                        max_fwd_pass=100,
                                        out_file=None):
    #初始话一些数组，后续计算使用
    fixed_mc = cost_el is None
    land_setup['forestable_area_available_init'] = np.array(land_setup['area_forestable_init']) - np.array(land_setup['area_reserved_init'])
    land_setup['forested_area_available_init'] = np.array(land_setup['area_forested_init']) - np.array(land_setup['area_reserved_init'])
    #ze表示论文中公式使用的z=0.2
    land_setup['ze'] = np.log(np.array(land_setup['species'])/np.array(land_setup['ae']))/np.log(np.array(land_setup['area_forestable_init']))
    
    #一共458个地区
    num_regions = len(land_setup)
    
    # Assemble matrix of cumulative development (if unchecked).
    # This will be derived from either a single fixed rate per ecoregion OR
    # a vector of rates per ecoregion, in which case def.rate is a list (not unlike cost)
    cum_devel_mat = np.array(land_setup['def_rate']).reshape((1,len(land_setup['def_rate']))).repeat(num_periods, axis=0)
    cum_devel_mat = np.cumsum(cum_devel_mat,axis=0)
    
    if num_periods == 1:
        cum_devel_mat = np.array(cum_devel_mat).reshape(num_regions, num_periods)
    #总的budget的初始化，当前为根据region数量平均初始化
    cur_budget = (1/num_regions) * (np.array(budget_per_period * np.ones((num_regions,num_periods))))

    #定义些条件，在循环内会重新初始化
    roi_diff = conv_tol*10
    num_pass_adjustments = 1
    last_race_over_time = cur_budget.shape[1]
    fwd_pass = 1
    total_budget_chg = 100000
    
    # 第一层循环：完全分配次数，总调节数额，最大调节步数
    # max.fwd.pass: limit on forward passes through the time horizon before exiting out
    while num_pass_adjustments > 0 and total_budget_chg > budget_conv_tol and fwd_pass <= max_fwd_pass:
        
        # last_source = -1
        last_target = -1
        recent_target = np.array(False).repeat(num_regions)
        recent_source = np.array(False).repeat(num_regions)
    
        t = 1
        # set up copy of budget for convergence checks. We'll look at the difference
        # across the entire time sweep.
        num_pass_adjustments = 0
        total_budget_chg = 0
        target_ignore_ctr = np.zeros(len(cur_budget))
        
        while t <= last_race_over_time:

            print(t)
            
            #获取target和source地区的budgets等于current budgets
            right_budget = cur_budget.copy()
            left_budget = cur_budget.copy()
            #论文中budget_nudge = 0.1 作用为为了减少计算式出现的不连续性和不精确性，人为制造一定的左右差异
            right_budget[:,t - 1] = right_budget[:,t - 1] + budget_nudge
            left_budget[:,t - 1] = np.maximum(left_budget[:,t - 1] - budget_nudge , 0)
            #带入现在的budgets计算ROI和race over等时间参数
            race_over_info = FindRaceEndConditionsProbabilistic(cur_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)
            race_over_info_right = FindRaceEndConditionsProbabilistic(right_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)
            race_over_info_left = FindRaceEndConditionsProbabilistic(left_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)
            #获得最大的截止时间
            last_race_over_time = max(race_over_info['race_end_periods']) + 1
            #获得左右的ROI值
            rois_right = race_over_info_right['rois_right'].copy()
            rois_left = race_over_info_left['rois_left'].copy()
            #获得当前periods所有地区的budgets
            period_budgets = cur_budget[:,t - 1].copy()
            
            #获取当前ce，这里没有和R一样，没有判断ce的类型，主要是因为ce在计算中一直为一维数据，没有发生改变
            period_ces = np.array(land_setup['ce']).copy()
            

            
            if fixed_mc:
                cost_of_rem_land = race_over_info['area_avail_raceendper'] * period_ces
            else:
                if cost_el != -1:
                    ep_oneplusep = cost_el/(1+cost_el)
                    if restoration_value > 0:
                        #计算剩余土地的cost，isoelastic cost
                        cost_of_rem_land = ep_oneplusep * period_ces*((np.array(land_setup['area_forestable_init']) - race_over_info['area_resrest_raceendper_start']) ** (1/ep_oneplusep))
                    else:
                        cost_of_rem_land = ep_oneplusep * period_ces*(np.array(land_setup['area_forestable_init'])
                                                                      - race_over_info['area_resrest_raceendper_start']) ** (1/ep_oneplusep) - (np.array(land_setup['area_forestable_init']) - race_over_info['area_resrest_raceendper_start'] + race_over_info['area_avail_raceendper'])**(1/ep_oneplusep)
                else:
                    raise Exception("cost elasticity of -1 not implemented")
            # Find out what extra budget we have: beyond whats needed to purchase land remaining in the 
            # final period (which may not be this period -- we include that condition below)这一步的计算可能产生小数点位数的精度问题
            #计算剩余所需的budgets，在下面作为是否归零的判断依据之一
            excess_budget =  period_budgets - cost_of_rem_land

            # The main ROI code doesn't set ROIs to zero when the race is over (because in earlier periods
            # we want to consider it to still be positive due to discreteness of time periods).
            # Set it to zero for regions in which this is the final period and we have more than enough
            # $ to conserve everything.
            #右侧（target）判定条件，余额大于等于0且当前periods还可进行保护设置为0
            mask_rois_right = (race_over_info_right['race_end_periods'] == t - 1) & (excess_budget >= 0)
            rois_right[mask_rois_right] = 0
            #左侧（source）判定条件，余额大于等于0且当前periods还可进行保护设置为0
            mask_rois_left = (race_over_info_left['race_end_periods'] == t - 1) & (excess_budget > 0)
            rois_left[mask_rois_left] = 0

            #已经race over的，表示不需要budgets了
            rois_right[race_over_info_right['race_end_periods'] < t - 1] = 0
            rois_left[race_over_info_left['race_end_periods'] < t - 1] = 0
            
            # Identify source region. Use left ROIs, and only consider regions that
            # have budget to give.
            #获得左侧ROI（source），min寻找最小
            roi_feas_left = rois_left.copy()
            roi_feas_left[period_budgets == 0] = np.nan
            source_region = np.nanargmin(roi_feas_left,axis=0)
            source_roi = roi_feas_left[source_region] 

            # Identify target region. Use right ROIs, and only consider regions where
            # the race is not already over as of this period given current budgets.
            #获得右侧ROI（target），max寻找最大
            roi_feas_right = rois_right.copy()
            roi_feas_right[race_over_info_right['race_end_periods'] < t - 1] = np.nan

            mask = (race_over_info_right['race_end_periods'] == t - 1) & (cost_of_rem_land <= period_budgets)

            roi_feas_right[mask] = np.nan
            roi_feas_right[target_ignore_ctr > 0] = np.nan
            #max寻找最大
            target_region = np.nanargmax(roi_feas_right,axis=0)

            target_roi = roi_feas_right[target_region]
            
            #当有最大和最小的ROI， 且source和target不是一个地区，两个ROI相差超过预设阈值，准备进行重新分配
            if len([source_region]) > 0 and len([target_region]) >  0 and source_region != target_region and source_roi < target_roi and (target_roi - source_roi) > conv_tol:
                num_pass_adjustments = num_pass_adjustments + 1

                print("Finding reallocation for {}->{}: {}, {}".format(source_region, target_region, source_roi, target_roi))
                realloc_amt = cur_budget[source_region,t - 1]
                #当超剩余土地话费大于现阶段budgets是，且超过值小于budgets是，使当前重新分配值为超过值
                if excess_budget[target_region] < 0 and -excess_budget[target_region] < realloc_amt:
                    realloc_amt = -excess_budget[target_region]
                
                #初始三个分配值等于重新分配总值
                last_realloc_delta = realloc_amt
                init_realloc_amt =  realloc_amt
                smallest_flip_realloc = realloc_amt

                
                #进行二分搜索，开始贪心算法重新分配budgets
                while True:
                    
                    trial_budget = cur_budget.copy()
                    #第一步直接把现有分配值全部分配个target地区
                    trial_budget[source_region,t - 1] = trial_budget[source_region,t - 1] - realloc_amt
                    #source地区减去分配值
                    trial_budget[target_region,t - 1] = trial_budget[target_region,t - 1] + realloc_amt

                    
                    right_trial_budget = trial_budget.copy()
                    left_trial_budget = trial_budget.copy()
                    #人为去调整左右budgets确保连续性和准确性
                    right_trial_budget[:,t - 1] = right_trial_budget[:,t - 1] + budget_nudge
                    left_trial_budget[:, t - 1] = np.maximum(left_trial_budget[:, t - 1 ] - budget_nudge , 0)
                    #重新计算ROI
                    trial_race_over_info = FindRaceEndConditionsProbabilistic(trial_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)

                    trial_race_over_info_left = FindRaceEndConditionsProbabilistic(left_trial_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)
                    trial_race_over_info_right = FindRaceEndConditionsProbabilistic(right_trial_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,t)
                    trial_rois_left = trial_race_over_info_left['rois_left'].copy()
                    trial_rois_right = trial_race_over_info_right['rois_right'].copy()
                    #得到当前实验的budegts（为初始的cur_budgets）
                    trial_period_budgets = trial_budget[:,t - 1].copy()
                    
                    if fixed_mc:
                        trial_cost_of_rem_land = trial_race_over_info['area_avil_raceendper'] * period_ces
                    else:
                        if cost_el != -1:
                            ep_oneplusep = cost_el/(1+cost_el)
                            #计算实验的剩余土地cost
                            trial_cost_of_rem_land = ep_oneplusep*period_ces*((np.array(land_setup['area_forestable_init']) - trial_race_over_info['area_resrest_raceendper_start'])**(1/ep_oneplusep)-
                                                                            (np.array(land_setup['area_forestable_init']) - (trial_race_over_info['area_resrest_raceendper_start']+trial_race_over_info['area_avail_raceendper']))**(1/ep_oneplusep))
                        else:
                            raise Exception("cost elasticity of -1 not implemented")
                    #计算实验剩余budgets
                    trial_excess_budget = trial_period_budgets - trial_cost_of_rem_land
                    #判别条件获取target和source的ROI
                    mask_right = (trial_race_over_info_right['race_end_periods'] == t - 1) & (trial_excess_budget >= 0)
                    trial_rois_right[mask_right] = 0
                    mask_left = (trial_race_over_info_left['race_end_periods'] == t - 1) & (trial_excess_budget > 0)
                    trial_rois_left[mask_left] = 0
                    trial_rois_right[trial_race_over_info_right['race_end_periods'] < t -1] = 0
                    trial_rois_left[trial_race_over_info_left['race_end_periods'] < t - 1] = 0
                    
                    trial_roi_target = trial_rois_right[target_region]
                    trial_roi_source = trial_rois_left[source_region]
                    #搜索条件
                    #当调整并计算后的targets ROI 大于 source ROI时：
                    if trial_roi_target > trial_roi_source:
                        if realloc_amt == init_realloc_amt:
                            #初始的调整大小等于调整总数，跳出搜索循环
                            break
                        if last_realloc_delta <= realloc_tol:
                            #重新分配的值等于最小分配阈值，跳出搜索
                            realloc_amt = smallest_flip_realloc
                            break
                        else:
                            #都不满足时，继续添加budgets，搜索补偿为一开始reallmoc_mat的一半
                            realloc_amt = realloc_amt + last_realloc_delta/2
                            last_realloc_delta = last_realloc_delta/2
                    #当调整并计算后的targets ROI 小于等于 source ROI时：
                    if trial_roi_target <= trial_roi_source:
                        #调整最小阈值
                        # What we're really after: shift that makes the target ROI 
                        # lower than the source ROI AFTER the transfer
                        # close to the source ROI AFTER the transfer
                        # higher than the source ROI BEFORE the transfer (so we get contraction of the ROI span)
                        # AND
                        # source ROI AFTER the transfer is no higher than target ROI BEFORE the transfer.
                        # This suggests a contraction of the ROIs among regions having spend, which should
                        # push us toward convergence.
                        if realloc_amt < smallest_flip_realloc:
                            smallest_flip_realloc = realloc_amt
                        if last_realloc_delta <= realloc_tol:
                            #最小阈值小于系统设定值时，跳出搜索
                            break
                        #当不满足继续搜索条件（ROI变化小于0.02左右并且实验target值大于原source值且实验source值小于原target值）跳出搜索
                        if ((trial_roi_source - trial_roi_target) <= inner_conv_tol*realloc_amt) and trial_roi_target >= source_roi and trial_roi_source <= target_roi:
                            #inner_conv_tols是0.1所以一直很小
                            break
                        else:
                            #未达到停止搜索条件，继续减小realloc_amt，继续搜索，搜索步长为原来的一半
                            # We flipped the ranking of the ROIs, but they're not very close. Reallocate less.
                            realloc_amt = realloc_amt - last_realloc_delta/2
                            last_realloc_delta = last_realloc_delta/2

                print("fp: {}; t={}: realloc {} {} {} {} {}".format(fwd_pass, t, source_region, target_region, realloc_amt, source_roi, target_roi))
                #更新当前budgets
                cur_budget[source_region,t - 1] = cur_budget[source_region,t - 1] - realloc_amt
                cur_budget[target_region,t - 1] = cur_budget[target_region,t - 1] + realloc_amt
                
                recent_target[target_region] = True
                recent_source[source_region] = True
                #更新总调整值
                total_budget_chg = total_budget_chg + realloc_amt
                
                
                target_ignore_ctr = np.maximum(target_ignore_ctr - 1,0)
                
                if source_region == last_target:
                    target_ignore_ctr[source_region] = NIGNOREITER
                
                last_source = source_region
                last_target = target_region
            
            else:
                #进入下一个periods
                t = t + 1
                recent_source = np.array(False).repeat(num_regions)
                recent_target = np.array(False).repeat(num_regions)           
        
        fwd_pass = fwd_pass + 1
    
    print('-----------------------------------------Results will be here-----------------------------------------')
    #输出results
    prelim_res = FindRaceEndConditionsProbabilistic(cur_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,1)
    final_rois = {'left':[],'right':[]}
    for i in range(0,num_regions):
        r_res = FindRaceEndConditionsProbabilistic(cur_budget,land_setup,cum_devel_mat,profiles,profile_abundance,cost_el,restoration_value,prelim_res['race_end_periods'][i])
        final_rois['left'].append(r_res['rois_left'][i])
        final_rois['right'].append(r_res['rois_right'][i])

    for ele in prelim_res:
        prelim_res[ele] = prelim_res[ele].tolist()
    result = {'budget':cur_budget.tolist(),'race_end':prelim_res,'land_setup':land_setup.to_json(orient='index',force_ascii=False),'profiles':profiles.tolist(),'profile_abundance':profile_abundance.tolist(),'cost_el':cost_el,'restoration_value':restoration_value}
        
    return result



def main():
    #读取数据集
    ecoregion_data = pd.read_csv('../zweiven-wwpf-eca8354/data/alldata.csv')
    profile_data = pd.read_csv('../zweiven-wwpf-eca8354/data/mammal_based_profiles.csv')
    profile_abundance = np.array(profile_data['filt.abundances'])
    profiles = np.array(profile_data.iloc[:,:len(profile_data.columns)-1])


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

    #计算初始面积
    forestable_area_init = []
    for i,n in enumerate(ecoregion_data['forest_area_2018']):
        forestable_area_init.append(n + 18 * ecoregion_data[kDefType][i])
    # Calculate area which *could* be reforested and its implications
    # for the SAR multiplier
    ecoregion_data['forestable_area_init'] = np.array(forestable_area_init)

    ecoregion_data['ae'] = np.exp(np.log(ecoregion_data['plant_spcs']) - kSARExponent * 
                            np.log(ecoregion_data['forestable_area_init']))

    # Drop regions with missing or clearly bad cost data (i.e. cost should be positive), and
    land_setup = ecoregion_data.dropna(axis=0,subset='GDPperm2')
    land_setup = ecoregion_data.loc[ecoregion_data['GDPperm2_partial'] > 0]
    land_setup = land_setup[['eco_code',
                         'GDPperm2_partial',
                         'forestable_area_init',
                         'ae',
                         'plant_spcs',
                         'forest_area_2018',
                         'protected_forest_area_2018',
                         'dozer_rate_2000_2018']]
    #更改key words名称
    land_setup.rename(columns = {'forestable_area_init' : 'area_forestable_init', 
                                    'GDPperm2_partial' : 'ce',
                                    'forest_area_2018' : 'area_forested_init',
                                    'protected_forest_area_2018' : 'area_reserved_init',
                                    'dozer_rate_2000_2018' : 'def_rate',
                                    'plant_spcs' : 'species'
                                    }, inplace = True)
    #计算花费
    land_setup['ce'] = land_setup['ce']/((np.array(land_setup['area_forestable_init']) - 
                                            np.array(land_setup["area_reserved_init"]))**(1/kCostEl))
    #数据处理，方便后续计算
    for key in land_setup:
        if key == 'eco_code':
            pass
        else:
            land_setup[key] = np.array(land_setup[key],dtype=np.float64)
    
    result = SolveBulldozersProbabilistic(land_setup,kNumPeriods,kBudgetPerPeriod,profiles,profile_abundance,kCostEl,kRestMultiplier,budget_conv_tol=100,conv_tol=5E-8,inner_conv_tol=1E-9,realloc_tol=0.1,budget_nudge=0.1,out_file='result')
    print(result)
    
if __name__ == '__main__':
    main()
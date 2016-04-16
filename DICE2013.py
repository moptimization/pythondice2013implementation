import pandas
import scipy
import pdb
import pyomo.opt
import pyomo.environ as pe
import logging


#using variable and parameter names that are consistent with the GAMS formulation available at the end of 
#http://www.econ.yale.edu/~nordhaus/homepage/documents/Dicemanualfull.pdf

#TODO: some parameters are hard-coded, but should be calculated based on the value of other parameters.  This way they can change appropriately with the 
#value of those other parameters

def getPopulationDict(maxPeriods,pop0,popasym,popadj):
	#create a dictionary of the projected population 
	population = {}
	population[1] = pop0
	for t in xrange(2,int(maxPeriods)+1):
		population[t] = population[t-1]*(popasym/population[t-1])**popadj
	return population

def getGrowthRateOfProductivityDict(ga0,dela,numPeriods,tstep):
	productivityGR = {} #could do this with scipy arrays to be faster, but this is sufficiently small not to matter
	for t in xrange(int(numPeriods)):
		productivityGR[t+1] = ga0 * scipy.exp(-dela*tstep*t)
	return productivityGR
	
def getLevelOfProductivity(a0,ga,numPeriods):
	productivity = {}
	productivity[1] = a0
	for t in xrange(2,int(numPeriods)+1):
		productivity[t] = productivity[t-1]/(1-ga[t-1])
	return productivity
	
def getCumulativeEfficiencyImprovement(gsigma1,dsig,tstep,numPeriods):
	efficiency = {}
	efficiency[1] = gsigma1
	for t in xrange(2,int(numPeriods)+1):
		efficiency[t] = efficiency[t-1]*((1+dsig)**tstep)
	return efficiency
	
def getGrowthRate(sig0,gsig,tstep,numPeriods):
	growthRate = {}
	growthRate[1] = sig0
	for t in xrange(2,int(numPeriods)+1):
		growthRate[t] = growthRate[t-1]*scipy.exp(gsig[t-1]*tstep)
	return growthRate
	
def getBackstopPrice(pback,gback,numPeriods):
	backstopPrice = {}
	for t in xrange(int(numPeriods)):
		backstopPrice[t+1] = pback*(1-gback)**t
	return backstopPrice
	
def getAdjustedCostForBackstop(pbacktime,sigma,expcost2,numPeriods):
	adjCost = {}
	for t in xrange(1,int(numPeriods)+1):
		adjCost[t] = pbacktime[t] * sigma[t] / expcost2 / 1000.
	return adjCost
	
def getEmmissionsFromDeforestation(eland0,deland,numPeriods):
	emmissions = {}
	for t in xrange(1,int(numPeriods)+1):
		emmissions[t] = eland0 * (1 - deland)**(t-1)
	return emmissions
	
def getAverageUtilitySocialDiscountRate(prstp,tstep,numPeriods):
	utilDisc = {}
	for t in xrange(1,int(numPeriods)+1):
		utilDisc[t] = 1./((1+prstp)**(tstep*(t-1)))
	return utilDisc
	
def getExogenousForcingOfOtherGreenhouseGases(fex0,fex1,numPeriods):
	exogForcing = {}
	for t in xrange(1,int(numPeriods)+1):
		if t < 19:
			exogForcing[t] = fex0 + (1./18.)*(fex1-fex0)*(t-1)
		else:
			exogForcing[t] = fex1
	return exogForcing
	
def getFractionOfEmmissionsInControlRegime(periodfullpart,partfractfull,partfract2010,numPeriods):
	partFract = {}
	for t in xrange(1,int(numPeriods)+1):
		if t <= periodfullpart:
			partFract[t] = partfract2010 + (partfractfull-partfract2010)*(t-1)/periodfullpart
		else:
			partFract[t] = partfractfull
	return partFract
	
def getTransientTSCCorrection(c10,c1beta,t2xco2):
	return c10 + c1beta*(t2xco2-2.9)
	
def getCarbonPrice(cprice0,gcprice,tstep,numPeriods):
	carbonPrice = {}
	for t in xrange(1,int(numPeriods)+1):
		carbonPrice[t] = cprice0*(1+gcprice)**(tstep*(t-1))
	return carbonPrice
	
def diceModel2013(**kwargs):
	model = pe.ConcreteModel()
	model.dual = pe.Suffix(direction=pe.Suffix.IMPORT)
	###SETS AND INDICES###
	model.time_periods = pe.Set( initialize=kwargs['rr'].keys()) 
	###PARAMETERS###
	#pdb.set_trace()
	###VARIABLES### 
	#Using the all-caps standard as set in the GAMS model
	def miuBounds(model,t):
		if t == 1:
			return (kwargs['miu0'],kwargs['miu0'])
		elif t < 30:
			return (0.,kwargs['limmiu'])
		else:
			return (0.,kwargs['limmiu']*kwargs['partfract'][t])
			
	def tatmBounds(model,t):
		if t == 1: 
			return (kwargs['tatm0'],kwargs['tatm0'])
		else:
			return (0,40.)
	
	def toceanBounds(model,t):
		if t == 1: 
			return (kwargs['tocean0'],kwargs['tocean0'])
		else:
			return (-1,20.)
			
	def matBounds(model,t):
		if t == 1: 
			return (kwargs['mat0'],kwargs['mat0'])
		else:
			return (10,scipy.inf)
	
	def muBounds(model,t):
		if t == 1: 
			return (kwargs['mu0'],kwargs['mu0'])
		else:
			return (100.,scipy.inf)
	
	def mlBounds(model,t):
		if t == 1: 
			return (kwargs['ml0'],kwargs['ml0'])
		else:
			return (1000.,scipy.inf)
	
	def cBounds(model,t):
		return (2.,scipy.inf)
	
	def kBounds(model,t):
		if t == 1: 
			return (kwargs['k0'],kwargs['k0'])
		else:
			return (1.,scipy.inf)
			
	def cpcBounds(model,t):
		return (.01,scipy.inf)
		
	def sBounds(model,t):
		if t <= kwargs['numPeriods'] - 10:
			return (-scipy.inf,scipy.inf)
		else:
			return (kwargs['optlrsav'],kwargs['optlrsav'])
			
	def ccaBounds(model,t):
		if t == 1:
			return (90,90)
		else:
			return (-scipy.inf,kwargs['fosslim'])
	
	model.MIU = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=miuBounds) #Emission control rate GHGs
	model.FORC = pe.Var(model.time_periods,domain=pe.Reals) #Increase in radiative forcing (watts per m2 from 1900)
	model.TATM = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=tatmBounds) #Increase temperature of atmosphere (degrees C from 1900)
	model.TOCEAN = pe.Var(model.time_periods,domain=pe.Reals,bounds=toceanBounds) #Increase temperatureof lower oceans (degrees C from 1900)
	model.MAT = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=matBounds) #Carbon concentration increase in atmosphere (GtC from 1750)
	model.MU = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=muBounds) #Carbon concentration increase in shallow oceans (GtC from 1750)
	model.ML = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=mlBounds) #Carbon concentration increase in lower oceans (GtC from 1750)
	model.E = pe.Var(model.time_periods,domain=pe.Reals) #Total CO2 emissions (GtCO2 per year)
	model.EIND = pe.Var(model.time_periods,domain=pe.Reals) #Industrial emissions (GtCO2 per year)
	model.C = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=cBounds) #Consumption (trillions 2005 US dollars per year)
	model.K = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=kBounds) #Capital stock (trillions 2005 US dollars)
	model.CPC = pe.Var(model.time_periods,domain=pe.NonNegativeReals,bounds=cpcBounds) #Per capita consumption (thousands 2005 USD per year)
	model.I = pe.Var(model.time_periods,domain=pe.NonNegativeReals) #Investment (trillions 2005 USD per year)
	model.S = pe.Var(model.time_periods,domain=pe.Reals,bounds=sBounds) #Gross savings rate as fraction of gross world product
	model.RI = pe.Var(model.time_periods,domain=pe.Reals) #Real interest rate (per annum)
	model.Y = pe.Var(model.time_periods,domain=pe.NonNegativeReals) #Gross world product net of abatement and damages (trillions 2005 USD per year)
	model.YGROSS = pe.Var(model.time_periods,domain=pe.NonNegativeReals) #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
	model.YNET = pe.Var(model.time_periods,domain=pe.Reals) #Output net of damages equation (trillions 2005 USD per year)
	model.DAMAGES = pe.Var(model.time_periods,domain=pe.Reals) #Damages (trillions 2005 USD per year)
	model.DAMFRAC = pe.Var(model.time_periods,domain=pe.Reals) #Damages as fraction of gross output
	model.ABATECOST = pe.Var(model.time_periods,domain=pe.Reals) #Cost of emissions reductions  (trillions 2005 USD per year)
	model.MCABATE = pe.Var(model.time_periods,domain=pe.Reals) #Marginal cost of abatement (2005$ per ton CO2)
	model.CCA = pe.Var(model.time_periods,domain=pe.Reals,bounds=ccaBounds) #Cumulative industrial carbon emissions (GTC)
	model.PERIODU = pe.Var(model.time_periods,domain=pe.Reals) #One period utility function
	model.CPRICE = pe.Var(model.time_periods,domain=pe.Reals) #Carbon price (2005$ per ton of CO2)
	model.CEMUTOTPER = pe.Var(model.time_periods,domain=pe.Reals) #Period utility
	
	model.UTILITY = pe.Var(domain=pe.Reals) #Welfare function
	
	def obj_rule(model):
		return  model.UTILITY
	model.OBJ = pe.Objective(rule=obj_rule, sense=pe.maximize)
	
	def dummy_rule(model):
		return (model.UTILITY <= 11000000)
	#model.dummy_rule = pe.Constraint(rule=dummy_rule)
	
	def emmissionsEquation(model,t):
		return (model.E[t] == model.EIND[t] + kwargs['etree'][t])
	model.emmissionsEquation = pe.Constraint(model.time_periods,rule=emmissionsEquation)
	
	def industrialEmmissions(model,t): #EIND(t)        =E= sigma(t) * YGROSS(t) * (1-(MIU(t)));
		return (model.EIND[t] == kwargs['sigma'][t] * model.YGROSS[t] * (1 - model.MIU[t]))
	model.industrialEmmissions = pe.Constraint(model.time_periods,rule=industrialEmmissions)
	
	def cumCarbonEmmissions(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else:
			return (model.CCA[t] == model.CCA[t-1] + model.EIND[t-1] * kwargs['tstep'] / 3.666)
	model.cumCarbonEmmissions = pe.Constraint(model.time_periods,rule=cumCarbonEmmissions)
	
	def radiativeForcing(model,t): 
		return (model.FORC[t] == kwargs['fco22x'] * (pe.log10(model.MAT[t]/588.000)/pe.log10(2)) + kwargs['forcoth'][t])
	model.radiativeForcing = pe.Constraint(model.time_periods,rule=radiativeForcing)
	
	def damageFraction(model,t):
		return (model.DAMFRAC[t] == (kwargs['a1'] * model.TATM[t]) + (kwargs['a2'] * model.TATM[t]**kwargs['a3']))
	model.damageFraction = pe.Constraint(model.time_periods,rule=damageFraction)
	
	def damagesConst(model,t):
		return (model.DAMAGES[t] == (model.YGROSS[t] * model.DAMFRAC[t]))
	model.damagesConst = pe.Constraint(model.time_periods,rule=damagesConst)
	
	def abatementCost(model,t):
		return (model.ABATECOST[t] == model.YGROSS[t] * kwargs['cost1'][t] * (model.MIU[t]**kwargs['expcost2']) * (kwargs['partfract'][t]**(1 - kwargs['expcost2'])))
	model.abatementCost = pe.Constraint(model.time_periods,rule=abatementCost)
	
	def mcAbatement(model,t): #constraint is giving an error.  solver or pyomo may not like exponents less than two.  
		myExp = kwargs['expcost2'] - 1
		if myExp >= 2:
			return (model.MCABATE[t] == kwargs['pbacktime'][t] * model.MIU[t]**myExp)
		else:
			myExp += 1
			return (model.MCABATE[t]*model.MIU[t] == kwargs['pbacktime'][t] * model.MIU[t]**myExp)
	#model.mcAbatement = pe.Constraint(model.time_periods,rule=mcAbatement)
	
	def carbonPriceEq(model,t): #constraint is giving an error.  solver or pyomo may not like fractional exponents less than two
		myExp = kwargs['expcost2'] - 1
		if myExp >= 2:
			return (model.CPRICE[t] == kwargs['pbacktime'][t] * (model.MIU[t]/kwargs['partfract'][t])**myExp)
		else:
			myExp += 1
			return (model.CPRICE[t]*model.MIU[t]/kwargs['partfract'][t] == kwargs['pbacktime'][t] * (model.MIU[t]/kwargs['partfract'][t])**myExp)
	#model.carbonPriceEq = pe.Constraint(model.time_periods,rule=carbonPriceEq)
	
	def atmosphericConcentration(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else: #MAT(t+1)       =E= MAT(t)*b11 + MU(t)*b21 + (E(t)*(tstep/3.666));
			return (model.MAT[t] == model.MAT[t-1]*kwargs['b11'] + model.MU[t-1] * kwargs['b21'] + model.E[t-1] * kwargs['tstep'] / 3.666)
	model.atmosphericConcentration = pe.Constraint(model.time_periods,rule=atmosphericConcentration)
	
	def lowerOceanConcentration(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else:
			return (model.ML[t] == model.ML[t-1]*kwargs['b33'] + model.MU[t-1] * kwargs['b23'] )
	model.lowerOceanConcentration = pe.Constraint(model.time_periods,rule=lowerOceanConcentration)
	
	def upperOceanConcentration(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else:
			return (model.MU[t] == model.ML[t-1]*kwargs['b32'] + model.MU[t-1] * kwargs['b22'] + model.MAT[t-1] * kwargs['b12'])
	model.upperOceanConcentration = pe.Constraint(model.time_periods,rule=upperOceanConcentration)
	
	def atmosphericTemperature(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else: 
			return (model.TATM[t] == model.TATM[t-1] + kwargs['c1'] * ((model.FORC[t] - (kwargs['fco22x']/kwargs['t2xco2'])*model.TATM[t-1]) \
				- (kwargs['c3'] * (model.TATM[t-1] - model.TOCEAN[t-1]))))
	model.atmosphericTemperature = pe.Constraint(model.time_periods,rule=atmosphericTemperature)
	
	def oceanTemperature(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else:
			return (model.TOCEAN[t] == model.TOCEAN[t-1] + kwargs['c4'] * (model.TATM[t-1] - model.TOCEAN[t-1]))
	model.oceanTemperature = pe.Constraint(model.time_periods,rule=oceanTemperature)
	
	def grossOutput(model,t):#not sure why this constraint does not produce an error
		#YGROSS(t)      =E= (al(t)*(L(t)/1000)**(1-GAMA))*(K(t)**GAMA);
		coeff = kwargs['al'][t]*(kwargs['l'][t]/1000)**(1-kwargs['gama'])
		return (model.YGROSS[t] == coeff*(model.K[t]**kwargs['gama']))
	model.grossOutput = pe.Constraint(model.time_periods,rule=grossOutput)
	
	def netOutput(model,t):
		return (model.YNET[t] == model.YGROSS[t] * (1-model.DAMFRAC[t]))
	model.netOutput = pe.Constraint(model.time_periods,rule=netOutput)
	
	def outputNetEqn(model,t):
		return (model.Y[t] == model.YNET[t] - model.ABATECOST[t])
	model.outputNetEqn = pe.Constraint(model.time_periods,rule=outputNetEqn)
	
	def consumptionEqn(model,t):
		return (model.C[t] == model.Y[t] - model.I[t])
	model.consumptionEqn = pe.Constraint(model.time_periods,rule=consumptionEqn)
	
	def test2(model,t):
		if t == 1:
			return model.K[t] == kwargs['k0']
		else:
			return pe.Constraint.Skip
	#model.test2 = pe.Constraint(model.time_periods,rule=test2)
	
	def perCapitaConsumption(model,t):
		return (model.CPC[t] == 1000 * model.C[t] / kwargs['l'][t])
	model.perCapitaConsumption = pe.Constraint(model.time_periods,rule=perCapitaConsumption)
	
	def savingsRate(model,t):
		return (model.I[t] ==  model.S[t] * model.Y[t])
	model.savingsRate = pe.Constraint(model.time_periods,rule=savingsRate)
	
	def capitalBalance(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else: #K(t+1)         =L= (1-dk)**tstep * K(t) + tstep * I(t);
			return (model.K[t] == (1 - kwargs['dk'])**kwargs['tstep'] * model.K[t-1] + kwargs['tstep'] * model.I[t-1])
	model.capitalBalance = pe.Constraint(model.time_periods,rule=capitalBalance)
	
	def interestRateEqn(model,t):
		if t == 1:
			return pe.Constraint.Skip
		else:
			return (model.RI[t] == (1 + kwargs['prstp']) * (model.CPC[t]/model.CPC[t-1])** (kwargs['elasmu']/kwargs['tstep']) - 1)
	model.interestRateEqn = pe.Constraint(model.time_periods,rule=interestRateEqn)
	
	def periodUtilityEqn(model,t):
		return (model.CEMUTOTPER[t] ==  model.PERIODU[t] * kwargs['l'][t] * kwargs['rr'][t])
	model.periodUtilityEqn = pe.Constraint(model.time_periods,rule=periodUtilityEqn)
	
	def instUtilityEqn(model,t):
		return (model.PERIODU[t] ==  ((model.C[t] * 1000 / kwargs['l'][t])**(1-kwargs['elasmu'])-1) / (1-kwargs['elasmu']) - 1)
	model.instUtilityEqn = pe.Constraint(model.time_periods,rule=instUtilityEqn)
	
	def utilityCalc(model):
		return (model.UTILITY == kwargs['tstep'] * kwargs['scale1'] * pe.summation(model.CEMUTOTPER) + kwargs['scale2'])
	model.utilityCalc = pe.Constraint(rule=utilityCalc)
	"""

*Utility
 util..               UTILITY        =E= tstep * scale1 * sum(t,  CEMUTOTPER(t)) + scale2 ;

	"""
	model.create()
	solver = pyomo.opt.SolverFactory('ipopt') #options are 'couenne', 'ipopt', and 'bonmin'.  ipopt is continuous and bonmin uses mixed integer
	results = solver.solve(model, tee=True, keepfiles=False, options="max_iter=2000")
	if (results.solver.status != pyomo.opt.SolverStatus.ok):
		logging.warning('Check solver not ok?')
	if (results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal):  
		logging.warning('Check solver optimality?')

	model.load(results)
	print 'Optimal solution value:', model.OBJ()
	return model
	
def createRow(variableName,variable,indices):
	result = [variableName]
	for idx in indices:
		result.append(variable[idx].value)
	return result
	
if __name__ == "__main__":
	paramDF = pandas.read_csv('diceParameters.csv')
	params = dict(zip(paramDF.key,paramDF.value))
	params['l'] = getPopulationDict(params['numPeriods'],params['pop0'],params['popasym'],params['popadj']) 
	params['ga'] = getGrowthRateOfProductivityDict(params['ga0'],params['dela'],params['numPeriods'],params['tstep'])
	params['al'] = getLevelOfProductivity(params['a0'],params['ga'],params['numPeriods'])
	params['gsig'] = getCumulativeEfficiencyImprovement(params['gsigma1'],params['dsig'],params['tstep'],params['numPeriods'])
	params['sigma'] = getGrowthRate(params['sig0'],params['gsig'],params['tstep'],params['numPeriods'])
	params['pbacktime'] = getBackstopPrice(params['pback'],params['gback'],params['numPeriods'])
	params['cost1'] = getAdjustedCostForBackstop(params['pbacktime'],params['sigma'],params['expcost2'],params['numPeriods'])
	params['etree'] = getEmmissionsFromDeforestation(params['eland0'],params['deland'],params['numPeriods'])
	params['rr'] = getAverageUtilitySocialDiscountRate(params['prstp'],params['tstep'],params['numPeriods'])
	params['forcoth'] = getExogenousForcingOfOtherGreenhouseGases(params['fex0'],params['fex1'],params['numPeriods'])
	params['partfract'] = getFractionOfEmmissionsInControlRegime(params['periodfullpart'],params['partfractfull'],params['partfract2010'],params['numPeriods'])
	params['c1'] = getTransientTSCCorrection(params['c10'],params['c1beta'],params['t2xco2'])
	params['cpricebase'] = getCarbonPrice(params['cprice0'],params['gcprice'],params['tstep'],params['numPeriods'])
	model = diceModel2013(**params)
	columns = ['Variable']
	columns.extend(params['l'].keys())
	df = pandas.DataFrame(columns=columns)
	df.loc[0] = createRow('MIU',model.MIU,sorted(params['l'].keys()))
	df.loc[1] = createRow('FORC',model.FORC,sorted(params['l'].keys()))
	df.loc[2] = createRow('TATM',model.TATM,sorted(params['l'].keys()))
	df.loc[3] = createRow('TOCEAN',model.TOCEAN,sorted(params['l'].keys()))
	df.loc[4] = createRow('MAT',model.MAT,sorted(params['l'].keys()))
	df.loc[5] = createRow('MU',model.MU,sorted(params['l'].keys()))
	df.loc[6] = createRow('ML',model.ML,sorted(params['l'].keys()))
	df.loc[7] = createRow('E',model.E,sorted(params['l'].keys()))
	df.loc[8] = createRow('EIND',model.EIND,sorted(params['l'].keys()))
	df.loc[9] = createRow('C',model.C,sorted(params['l'].keys()))
	df.loc[10] = createRow('K',model.K,sorted(params['l'].keys()))
	df.loc[11] = createRow('CPC',model.CPC,sorted(params['l'].keys()))
	df.loc[12] = createRow('I',model.I,sorted(params['l'].keys()))
	df.loc[13] = createRow('S',model.S,sorted(params['l'].keys()))
	df.loc[14] = createRow('RI',model.RI,sorted(params['l'].keys()))
	df.loc[15] = createRow('Y',model.Y,sorted(params['l'].keys()))
	df.loc[16] = createRow('YGROSS',model.YGROSS,sorted(params['l'].keys()))
	df.loc[17] = createRow('YNET',model.YNET,sorted(params['l'].keys()))
	df.loc[18] = createRow('DAMAGES',model.DAMAGES,sorted(params['l'].keys()))
	df.loc[19] = createRow('DAMFRAC',model.DAMFRAC,sorted(params['l'].keys()))
	df.loc[20] = createRow('ABATECOST',model.ABATECOST,sorted(params['l'].keys()))
	df.loc[21] = createRow('MCABATE',model.MCABATE,sorted(params['l'].keys()))
	df.loc[22] = createRow('CCA',model.CCA,sorted(params['l'].keys()))
	df.loc[23] = createRow('PERIODU',model.PERIODU,sorted(params['l'].keys()))
	df.loc[24] = createRow('CPRICE',model.CPRICE,sorted(params['l'].keys()))
	df.loc[25] = createRow('CEMUTOTPER',model.CEMUTOTPER,sorted(params['l'].keys()))
	df.to_csv('results.csv',index=False)
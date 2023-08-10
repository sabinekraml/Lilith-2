# ======================================================================
#			importing libraries
# ======================================================================
import sys, os
import shutil
import re
from scipy.interpolate import griddata
from scipy import optimize 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

print()    
print("".ljust(20,'*')+" STARTING ".ljust(35,"*"))
print()
lilith_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(lilith_dir)
sys.path.append('../..')
import lilith
lilithcalc = lilith.Lilith(verbose=False, timer=False)


print("     Import Lilith's core ".ljust(40,".")+" ok!".rjust(10,"."))

# ====================================================================== 
#			declaring input, scan region, output  
# ======================================================================

smeft_input = "data/smeft.list"  

#fit parameter 
fit_params = {}

fit_params['c3Hq'] = 'p01'
fit_params['cdH']  = 'p02'
fit_params['ceH']  = 'p03'
fit_params['cHe'] = ' (-0.62 * p04)'
fit_params['c1Hl'] = ' 0.78 * p04'
fit_params['c3Hl'] = ' (-0.83 * p05)'
fit_params['cllprime'] = ' 0.55 * p05'
fit_params['cHd'] = '(-0.26 * p06 + 0.24 * p07 )'
fit_params['cHu'] = ' (0.87 * p06 - 0.37 * p07 )'
fit_params['c1Hq'] = '(- 0.42 * p06 - 0.9 * p07)'
fit_params['cHB'] = '( - 0.84 * p08 - 0.31 * p09 + 0.43 * p10 )'
fit_params['cHW'] = '(- 0.27 * p08 + 0.95 * p09 + 0.17 * p10)'
fit_params['cHWB'] = '( 0.47 * p08 - 0.02 * p09 + 0.88 * p10)'
fit_params['cHDD'] = '( - 0.05 * p09 + 0.13 * p10)'
fit_params['cW'] = '(- 0.02 * p08 - 0.01 * p09 + 0.03 * p10)'
fit_params['cuB'] = '(- 0.05 * p08 - 0.04 * p09 + 0.07 * p10)'
fit_params['cuW'] = '(- 0.02 * p08 - 0.01 * p09 + 0.04 * p10)'
fit_params['cHG'] = '( p11 + 0.04 * p12)'
fit_params['cuG'] = '(0.04 * p11 - p12)'
fit_params['cuH'] = '(- 0.09 * p12)'
fit_params['cG'] = '(-0.2 * p13)'
fit_params['c8qd'] = '(-0.05 * p13)'
fit_params['c1qq'] = '(-0.02 * p13)'
fit_params['cqq'] = '(-0.38 * p13)'
fit_params['c3qq'] = '(-0.08 * p13)'
fit_params['c31qq'] = '(-0.78 * p13)'
fit_params['c1qu'] = '(-0.01 * p13)'
fit_params['c8qu'] = '(-0.23 * p13)'
fit_params['c8ud'] = '(-0.05 * p13)'
fit_params['cuu'] = '(-0.02 * p13)'
fit_params['c1uu'] = '(-0.37 * p13)'




#number of grid steps in each of the two dimensions 
grid_subdivisions = 3

# Lilith precision mode
my_precision = "BEST-QCD"

# Higgs mass to test
hmass = 125

# Output files
if (not os.path.exists("results")):
    os.mkdir("results")
output = "results/SMEFT_2d.out"
outputplot = "results/SMEFT_2d.pdf"

officialplot = ''

# print testing
# ~ print("fit-param dict = ",fit_params)
# ~ paramlist = list(fit_params.keys())
# ~ print(paramlist)

# ======================================================================
# 			Loading and import input 
# ======================================================================

lilithcalc.readsmeftinput(smeft_input,fit_params)  
lilithcalc.readexpinput(lilithcalc.smeft_stxs)

# ======================================================================
#		    Optimize for multi-variables minimum   
# ======================================================================

# ----------- optimize minimal --- rolldown method
print()
print("".ljust(10," ")+" Optimizing ".rjust(20,"*")+"".rjust(10,"*"))
print()

def get_ind_param_list(param_dict):
    ind_param_list = []
    for item in param_dict:
        ind_param_list += re.findall('p\w+',param_dict[item])
    ind_param_list = list(dict.fromkeys(ind_param_list))
    return ind_param_list
def get_ini_value_vec(param_list):
    ivalue_list = [0] * len(param_list)
    return ivalue_list
def par_list2dict(value_list):
    param_arg = {}
    for i in range(1,len(value_list)+1):
        d = 'p' + '%02d' %i 
        param_arg[d] = value_list[i-1]
    return param_arg
def compL(ivalue_list):
    param_arg = par_list2dict(ivalue_list)
    lilithcalc.smeftmu_eval(param_arg)
    lilithcalc.computelikelihoodsmeft()
    result_l = lilithcalc.smeft_l
    return result_l

ind_param_list = get_ind_param_list(fit_params)
param_vec_i = get_ini_value_vec(ind_param_list)


# ~ # downhill method
print("".rjust(5," ")+ "Downhill method ".ljust(37,".")+" running", end="\r")
res_downhill = optimize.fmin(compL,param_vec_i, disp=False)
param_res_downhill_dict = par_list2dict(res_downhill)
res_dict_downhill = par_list2dict(res_downhill)
print("".rjust(5," ")+ "Downhill method ".ljust(39,".")+" done!")
    
# ~ # powell method
print("".rjust(5," ")+ "Powell method ".ljust(37,".")+" running", end="\r") 
res_powell = optimize.fmin_powell(compL,param_vec_i, disp=False)
param_res_powell_dict = par_list2dict(res_powell)
res_dict_powell = par_list2dict(res_powell)
print("".rjust(5," ")+ "Powell method ".ljust(39,".")+" done!")

# ~ # BFGS method
print("".rjust(5," ")+ "BGFS method ".ljust(37,".")+" running", end="\r") 
res_bfgs = optimize.fmin_bfgs(compL,param_vec_i, disp=False)
param_res_bfgs_dict = par_list2dict(res_bfgs)
res_dict_bgfs = par_list2dict(res_bfgs)
print("".rjust(5," ")+ "BGFS method ".ljust(39,".")+" done!")

# ~ # ---------- print results --------------------------
print()
print("".ljust(5," ") + "RESULTS") 
print("".ljust(5," ") + "Methods" + "       >>>  " \
     + " downhill".rjust(13," ")\
     + " Powell".rjust(13," ")\
     + " BFGS".rjust(13," ")\
     )
print("".ljust(4," ")+"".rjust(60,"-"))   
print( " min(-2*logL)".rjust(17," ") + "   =   " \
    + " %12.5f"%compL(res_downhill)\
    + " %12.5f"%compL(res_powell)\
    + " %12.5f"%compL(res_bfgs)\
     ) 
print(" at ".rjust(10," "))
for item in ind_param_list:
    print( item.rjust(17," ") + "   =   " \
        + " %12.5f"%res_dict_downhill[item] \
        + " %12.5f"%res_dict_powell[item] \
        + " %12.5f"%res_dict_bgfs[item] \
        )



lilithcalc.smeft_cleanlist()

print()
print("".ljust(22,'-')+" DONE ".ljust(33,"-"))
print()

# ~ # ------------- testing area - 1 sigma calculation -------


#  ----------------------------------------------------------- 

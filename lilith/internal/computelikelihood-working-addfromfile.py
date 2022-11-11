##########################################################################
#
#  This file is part of Lilith
#  v1 (2015) by Jeremy Bernon and Beranger Dumont 
#  v2 (2019) by Sabine Kraml, Tran Quang Loc, Dao Thi Nhung, Le Duc Ninh 
#            converted to Python 3 by Marius Bertrand (Jul/Aug 2020)
#
#  Web page: http://lpsc.in2p3.fr/projects-th/lilith/
#
#  In case of questions email sabine.kraml@lpsc.in2p3.fr 
#
#
#    Lilith is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lilith is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lilith.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

from ..errors import LikelihoodComputationError
import numpy as np

def compute_likelihood(exp_mu, user_mu, user_mode):
    """Computes the likelihood from experimental mu and user mu."""
    likelihood_results = []
    l = 0. # actually -2log(likelihood)
    for mu in exp_mu:
        # compute user mu value scaled to efficiencies
        user_mu_effscaled = {}
        try:
            if mu["dim"] == 1:
                user_mu_effscaled["x"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["x"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["x"] += eff_prod*user_mu[prod,decay]
            elif mu["dim"] == 2:
                user_mu_effscaled["x"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["x"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["x"] += eff_prod*user_mu[prod,decay]

                user_mu_effscaled["y"] = 0.
                for (prod,decay),eff_prod in list(mu["eff"]["y"].items()):
                    if mu["sqrts"] not in ["1.96","7","8","7.","7.0","8.0","8.","7+8"]\
                      and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                        if user_mode == "reducedcouplings":
                            prod = prod + "13"
                    user_mu_effscaled["y"] += eff_prod*user_mu[prod,decay]
            elif mu["dim"] >= 3:
                for i in range(1,mu["dim"]+1):
                    d = "d" + str(i)
                    user_mu_effscaled[d] = 0.
                    for (prod,decay),eff_prod in list(mu["eff"][d].items()):
                        if mu["sqrts"] not in ["1.96","7","8","7.","8.","7.0","8.0","7+8"]\
                          and (prod == "ggH" or prod == "VBF" or prod == "tHq" or prod == "tHW" or prod == "ggZH"):
                            if user_mode == "reducedcouplings":
                                prod = prod + "13"
                        user_mu_effscaled[d] += eff_prod*user_mu[prod,decay]

        except KeyError as s:
            if "s" in ["eff", "x", "y"]:
                # the experimental mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in exp_mu: key "' + str(s) +
                    '" is not found')
            else:
                # the total user mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in user_mu_tot: key "' +
                    str(s) + '" is not found')

        try:
            # likelihood computation in case of a type="normal" (odinary Gaussian approximation)
            if mu["type"] == "n":
                if mu["dim"] == 1:
                    if user_mu_effscaled["x"] < mu["bestfit"]["x"]:
                        unc = mu["param"]["uncertainty"]["left"]
                    else:
                        unc = mu["param"]["uncertainty"]["right"]
                    cur_l = ((mu["bestfit"]["x"] - user_mu_effscaled["x"])**2/unc**2)

                elif mu["dim"] == 2:
                    a = mu["param"]["a"]
                    b = mu["param"]["b"]
                    c = mu["param"]["c"]

                    cur_l = a*(mu["bestfit"]["x"] - user_mu_effscaled["x"])**2
                    cur_l += c*(mu["bestfit"]["y"] - user_mu_effscaled["y"])**2
                    cur_l += (2*b*(mu["bestfit"]["x"] - user_mu_effscaled["x"])
                             * (mu["bestfit"]["y"] - user_mu_effscaled["y"]))

                elif mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])

                    cur_l = mu["param"]["inv_cov_m"].dot(mu_vec).dot(mu_vec.T)

            # likelihood computation in case of a type="variable normal"
            # following "Variable Gaussian 2", Barlow arXiv:physics/0406120v1, Eq. 18
            if mu["type"] == "vn":
                if mu["dim"] == 1:
                    unc_left = abs(mu["param"]["uncertainty"]["left"])
                    unc_right = mu["param"]["uncertainty"]["right"]
                    if unc_left == 0:
                        cur_l = np.abs(user_mu_effscaled["x"]-mu["bestfit"]["x"])/unc_right
                    elif unc_right == 0:
                        cur_l = np.abs(user_mu_effscaled["x"]-mu["bestfit"]["x"])/unc_left
                    else:
                        num = user_mu_effscaled["x"] - mu["bestfit"]["x"]
                        den = unc_left*unc_right + (unc_right - unc_left)*num
                        if den == 0:
                            raise LikelihoodComputationError(
                              'divided by zero in 1D Variable Gaussian case')
                        else:
                            cur_l = num**2/np.abs(den)
                if mu["dim"] == 2:
                    p = mu["param"]["correlation"]
                    sig1p = mu["param"]["uncertainty"]["x"]["right"]
                    sig1m = abs(mu["param"]["uncertainty"]["x"]["left"])
                    sig2p = mu["param"]["uncertainty"]["y"]["right"]
                    sig2m = abs(mu["param"]["uncertainty"]["y"]["left"])
                    z10 = mu["bestfit"]["x"]
                    z20 = mu["bestfit"]["y"]
                    z1 = user_mu_effscaled["x"]
                    z2 = user_mu_effscaled["y"]
                    V1 = sig1p*sig1m
                    V1e = sig1p - sig1m
                    V2 = sig2p*sig2m
                    V2e = sig2p - sig2m
                    V1f = np.abs(V1 + V1e*(z1-z10))
                    V2f = np.abs(V2 + V2e*(z2-z20))
                    cur_l = 1.0/(1-p**2)*((z1-z10)**2/V1f-2*p*(z1-z10)*(z2-z20)/np.sqrt(V1f*V2f)+(z2-z20)**2/V2f)
                elif mu["dim"] >= 3:
                    try:
                       xsbr_input = "validations/ATLAS/HIGG-2021-23-Oct22/Run2-ATLAS-HIGG-2021-23-combine-xsbr.txt"  
                       xsbrdata = np.genfromtxt(xsbr_input,comments="#")                       
                       A_vec = np.array([xsbrdata[:,0]])
                       error_ex_m = np.array([xsbrdata[:,1]])
                       error_ex_p = np.array([xsbrdata[:,2]])
                       B_vec = np.array([xsbrdata[:,3]])
                       error_th_m = np.array([xsbrdata[:,4]])
                       error_th_p = np.array([xsbrdata[:,5]])
#                       print("loading xs*br data... success")
                       xsbrcount = 1
                    except IOError:
                       xsbrcount = 0
#                       print("no xs*br data found")
                    if xsbrcount == 0: 
# original calculation                
                       mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                       for i in range(4,mu["dim"]+1):
                           d = "d"+str(i)
                           mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])

                       unc_sym = np.sqrt(np.abs(mu["param"]["VGau"] + mu["param"]["VGau_prime"]*mu_vec))
                       cov_m = unc_sym*mu["param"]["corr_m"]*unc_sym.T
                       cov_m_tot = cov_m 
                    elif xsbrcount == 1:
                       
## declare data for experimental and theoretical XS.BR for HIGG-2021-23                 
#                    A_vec = np.array([[0.105429,0.0107858,0.00421673,-0.000399613,0.00102226,0.000502212,1.12479,0.123055,0.0796878,0.025986,10.9821,0.849887,\
#                    	0.592278,0.490098,0.206185,2.51406,0.218468,0.126317,0.0502151,0.747513,0.462761,0.118169,27.4805,0.00528577,0.00276746]])
#                    error_ex_m = np.array([[0.0100745,0.00210465,0.0014003,0.000976664,0.00035011,0.000650091,0.123759,0.0399631,0.0498799,0.0170827,1.21176,0.138208,\
#                    	0.267322,0.228647,0.0771942,0.717863,0.039981,0.0748551,0.02753,0.186493,0.104315,0.112019,10.1707,0.00870713,0.00150411]])
#                    error_ex_p = np.array([[0.0104646,0.00235856,0.0015465,0.00109426,0.000371587,0.000813171,0.129834,0.0480131,0.0619694,0.025969,1.25424,0.146516,\
#                    	0.316977,0.314519,0.0817104,0.814827,0.045238,0.0785728,0.0314259,0.200307,0.113702,0.115901,10.5284,0.00876469,0.0015907]])
#                    B_vec = np.array([[0.101598,0.00794477,0.00275964,0.0018069,0.00113455,0.000192426,1.18203,0.0924324,0.0531287,0.0154385,9.63171,\
#                    	0.753178,0.261619,0.171297,0.125799,2.8,0.218954,0.125851,0.0365706,0.7062,0.462391,0.339576,28.0324,0.00983912,0.00119602]])
#                    error_th_p = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])    # use when testing for the case of no theoretical uncertainty
#                    error_th_m = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])    # use when testing for the case of not theoretical uncertainty 
# theoretical error with pdf error            
#                    error_th_p = np.array([[5.3,61,500,140,0.21,2.3,19,6,50,0.13,1.4,31,0.08,8,2.2,20]])
#                    error_th_m = np.array([[5.3,61,500,140,0.21,2.3,19,6,50,0.12,1.4,29,0.11,12,3.1,29]])
# theoretical error without pdf error                    
#                    error_th_m = np.array([[0.00655834,0.000273996,0.0000955456,0.0000841382,0.000117554,0.0000251025,0.0731961,0.002706,0.00152606,\
#                    	0.00138693,0.596397,0.0220437,0.00769851,0.00733883,0.011301,0.174984,0.00667004,0.00376717,0.00329977,0.0187373,\
#                    	0.018916,0.0301977,1.58646,0.000607944,0.0000324972]])
#                    error_th_p = np.array([[0.00654805,0.000278667,0.0000967425,0.0000848243,0.000118148,0.0000134556,0.0726344,0.00269864,\
#                    	0.00152158,0.00136889,0.591797,0.0219799,0.00762871,0.00731812,0.0111538,0.173756,0.00666768,0.00376542,\
#                  	0.00325806,0.0183345,0.0187762,0.0297686,1.57043,0.000602992,0.0000323837]])
#                    f_vec = A_vec/B_vec
#                    f_vec = np.array([[0.96,1.04,1.08,0.96,1.39,2.68,0.59,1.16,3.01,1.09,0.68,1.19,1.1,1.5,1.38,0.79]])
## declare data for HIGG-2018-28
#                    A_vec = np.array([[1.12,0.11,0.075,0.026]])
#                    error_ex_p = np.array([[0.13,0.04,0.06,0.03]])
#                    error_ex_m = np.array([[0.13,0.04,0.05,0.02]])
#                    B_vec = np.array([[1.17,0.092,0.0524,0.0154]])
#                    error_th_p = np.array([[0.08,0.002,0.0027,0.001]])
#                    error_th_m = np.array([[0.08,0.002,0.0049,0.0013]])	
# theoretical error with pdf error 
#                    error_th_p = np.array([[0.06545401, 0.02599244, 0.07375007, 0.10963063]])
#                    error_th_m = np.array([[0.06545401, 0.02599244, 0.13384272, 0.14251982]])

                        f_vec = A_vec/B_vec    	
#                    print("A_vec = ",A_vec.T) 
                        mu_vec = np.array([user_mu_effscaled["d1"] - f_vec[0][0],
                                       user_mu_effscaled["d2"] - f_vec[0][1],
                                       user_mu_effscaled["d3"] - f_vec[0][2]])
                        for i in range(4,mu["dim"]+1):
                           d = "d"+str(i)
                           mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - f_vec[0][i-1]])
                     
                        unc_ex_p = error_ex_p/B_vec
                        unc_ex_m = error_ex_m/B_vec
                        unc_th_p = (error_th_p/B_vec)*f_vec
                        unc_th_m = (error_th_m/B_vec)*f_vec
                        VGau_ex = unc_ex_p*unc_ex_m
                        VGau_ex_prime = unc_ex_p - unc_ex_m
                        Sigma_ex = np.sqrt(np.abs(VGau_ex + VGau_ex_prime*mu_vec))
                        VGau_th = unc_th_p*unc_th_m
                        VGau_th_prime = unc_th_p - unc_th_m
                        Sigma_th = np.sqrt(np.abs(VGau_th + VGau_th_prime*mu_vec))
                    
                        cov_ex = Sigma_ex*mu["param"]["corr_m"]*Sigma_ex.T
                        diag_mat = np.identity(mu["dim"])
                        cov_th = Sigma_th*diag_mat*Sigma_th.T
                        
                        cov_m_tot = cov_ex + cov_th
                    
#                    Binv=1.0/B_vec												
#                    cov_m_tot = Binv*cov_ex*Binv.T + f_vec*(Binv*cov_th*Binv.T)*f_vec.T		# this is wrong 
                    

#                    print("f_vec.T = ",f_vec.T)
#                    print("B_vec.T = ",B_vec.T) 

## the likelihood from covariance
                    inv_cov_m = np.linalg.inv(cov_m_tot)
                    cur_l = (mu_vec.T).dot(inv_cov_m).dot(mu_vec)

            # likelihood computation in case of a type="variable normal 1"
            # following "Variable Gaussian 1", Barlow arXiv:physics/0406120v1, Eq. 15
            if mu["type"] == "vn1":
                if mu["dim"] == 1:
                    unc_left = abs(mu["param"]["uncertainty"]["left"])
                    unc_right = mu["param"]["uncertainty"]["right"]
                    if unc_left == 0:
                        cur_l = 1.0
                    elif unc_right == 0:
                        cur_l = 1.0
                    else:
                        num = user_mu_effscaled["x"] - mu["bestfit"]["x"]
                        den = (2*unc_left*unc_right + (unc_right - unc_left)*num)/(unc_right + unc_left)
                        if den == 0:
                            raise LikelihoodComputationError(
                              'divided by zero in 1D Variable Gaussian 1 case')
                        else:
                            cur_l = (num/den)**2
                if mu["dim"] == 2:
                    p = mu["param"]["correlation"]
                    sig1p = mu["param"]["uncertainty"]["x"]["right"]
                    sig1m = abs(mu["param"]["uncertainty"]["x"]["left"])
                    sig2p = mu["param"]["uncertainty"]["y"]["right"]
                    sig2m = abs(mu["param"]["uncertainty"]["y"]["left"])
                    z10 = mu["bestfit"]["x"]
                    z20 = mu["bestfit"]["y"]
                    z1 = user_mu_effscaled["x"]
                    z2 = user_mu_effscaled["y"]
                    V1s = sig1p + sig1m
                    V1 = 2*sig1p*sig1m/V1s
                    V1e = (sig1p - sig1m)/V1s
                    V2s = sig2p + sig2m
                    V2 = 2*sig2p*sig2m/V2s
                    V2e = (sig2p - sig2m)/V2s
                    V1f = (V1 + V1e*(z1-z10))**2
                    V2f = (V2 + V2e*(z2-z20))**2
                    cur_l = 1.0/(1-p**2)*((z1-z10)**2/V1f-2*p*(z1-z10)*(z2-z20)/np.sqrt(V1f*V2f)+(z2-z20)**2/V2f)
                elif mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])

                    unc_sym = mu["param"]["SGau"] + mu["param"]["SGau_prime"]*mu_vec
                    cov_m = unc_sym*mu["param"]["corr_m"]*unc_sym.T
##                    print("cor_m =",mu["param"]["corr_m"])
##                    print("cov_m =",cov_m)
## full theoretical errors (incl pdf):
#                    error_th_p = np.array([0.06545401, 0.02599244, 0.07375007, 0.10963063])
#                    error_th_m = np.array([0.06545401, 0.02599244, 0.13384272, 0.14251982])
## theoretical errors without pdf:
#                    error_th_p = np.array([0.0629283, 0.00672014, 0.07188083, 0.10060985])
#                    error_th_m = np.array([0.0629283, 0.00672014, 0.13282189, 0.13570322])
#                    mu_th_VGau_sum = error_th_p + error_th_m
#                    mu_th_SGau = 2*error_th_p*error_th_m/mu_th_VGau_sum
#                    mu_th_SGau_prime = (error_th_p - error_th_m)/mu_th_VGau_sum
#                    unc_sym_th = mu_th_SGau + mu_th_SGau_prime*mu_vec
#                    corr_m_th = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
#                    cov_m_th = unc_sym_th*corr_m_th*unc_sym_th.T

## ggF+bbH, VBF, VH (WH+ZH), ttH+tH
#                    error_pdf_p = np.array([0.01800715, 0.0251087, 0.0164991, 0.04354921])
#                    error_pdf_m = np.array([0.01800715, 0.0251087, 0.0164991, 0.04354921])
#                    mu_pdf_VGau_sum = error_pdf_p + error_pdf_m
#                    mu_pdf_SGau = 2*error_pdf_p*error_pdf_m/mu_pdf_VGau_sum
#                    mu_pdf_SGau_prime = (error_pdf_p - error_pdf_m)/mu_pdf_VGau_sum
#                    unc_sym_pdf = mu_pdf_SGau + mu_pdf_SGau_prime*mu_vec
##                    corr_m_pdf = np.array([[1,-1,-1,1],[-1,1,1,-1],[-1,1,1,-1],[1,-1,-1,1]])
#                    corr_m_pdf = np.array([[1,0,0,0],[0,1,1,0],[0,1,1,0],[0,0,0,1]])
#                    cov_m_pdf = unc_sym_pdf*corr_m_pdf*unc_sym_pdf.T

                    cov_m_tot = cov_m # default option
#                    cov_m_tot = cov_m + cov_m_th + cov_m_pdf # include theoretical, pdf errors
##                    print("cor_m_pdf =",corr_m_pdf)
##                    print("cov_m_pdf =",cov_m_pdf)
##                    print("cov_m_tot =",cov_m_tot)
                    inv_cov_m = np.linalg.inv(cov_m_tot)
                    cur_l = inv_cov_m.dot(mu_vec).dot(mu_vec.T)
                    

            # likelihood computation in case of a type="Poisson"
            # following "Generalised Poisson" of Barlow, arXiv:physics/0406120v1, Eq. 10a
            if mu["type"] == "p":
                if mu["dim"] == 1:
                    sigm = abs(mu["param"]["uncertainty"]["left"])
                    sigp = mu["param"]["uncertainty"]["right"]
                    x0 = mu["bestfit"]["x"]
                    x = user_mu_effscaled["x"]

                    gamma = mu["param"]["gamma"]
                    nu = mu["param"]["nu"]
                    alpha = nu*gamma
                    cur_l = -alpha*(x-x0) + nu*np.log(1+alpha*(x-x0)/nu)
                    cur_l = -2.*cur_l

                if mu["dim"] == 2:
                    p = mu["param"]["correlation"]
                    sig1p = mu["param"]["uncertainty"]["x"]["right"]
                    sig1m = abs(mu["param"]["uncertainty"]["x"]["left"])
                    sig2p = mu["param"]["uncertainty"]["y"]["right"]
                    sig2m = abs(mu["param"]["uncertainty"]["y"]["left"])
                    z10 = mu["bestfit"]["x"]
                    z20 = mu["bestfit"]["y"]
                    z1 = user_mu_effscaled["x"]
                    z2 = user_mu_effscaled["y"]

                    gamma1 = mu["param"]["gamma"]["x"]
                    gamma2 = mu["param"]["gamma"]["y"]

                    nu1 = mu["param"]["nu"]["x"]
                    alpha1 = nu1*gamma1
                    nu2 = mu["param"]["nu"]["y"]
                    alpha2 = nu2*gamma2
                    A = mu["param"]["A_corr"]
                    alpha = mu["param"]["alpha_corr"]
                    L2t1 = -alpha1*(z1-z10) + nu1*np.log(1+alpha1*(z1-z10)/nu1)
                    L2t2a = -alpha2*(z2 - z20 + 1/gamma2)*np.exp(alpha*nu1 - A*alpha1*(z1 - z10 + 1/gamma1))
                    L2t2b = -alpha2*(1/gamma2)*np.exp(alpha*nu1 - A*alpha1/gamma1)
                    L2t2c = nu2*np.log(L2t2a/L2t2b)
                    L2t2 = L2t2a - L2t2b + L2t2c
                    cur_l = -2.0*(L2t1 + L2t2)

            # likelihood computation in case of a type="full" (exact likelihood provided in terms of a grid file)
            if mu["type"] == "f":
                if mu["dim"] == 1:
#                    cur_l = mu["Lxy"](user_mu_effscaled["x"]) - mu["LChi2min"]
                    cur_l = mu["Lxy"](user_mu_effscaled["x"])
                elif mu["dim"] == 2:
#                    cur_l = mu["Lxy"](user_mu_effscaled["x"],user_mu_effscaled["y"])[0][0] - mu["LChi2min"]
                    cur_l = mu["Lxy"](user_mu_effscaled["x"],user_mu_effscaled["y"])[0][0]
                if cur_l < 0:
                    cur_l = 0.
        except KeyError as s:
            raise LikelihoodComputationError(
                'there are missing elements in exp_mu: key "' + s +
                '" is not found')

# LDN added a control on the value of loglikelihood
#        if cur_l < 0:
##            if mu["dim"] == 1 and mu["type"] == "f":
##                print("mu =",user_mu_effscaled["x"])
##                print("cur_l min =",mu["LChi2min"])
##                print("cur_l x   =",mu["Lxy"](user_mu_effscaled["x"]))
##                print("cur_l =",cur_l)
#            raise LikelihoodComputationError(
#                'loglikelihood is negative, check the value of mu')
# end


        l += cur_l
        likelihood_results.append(
            {"experiment": mu["experiment"], "source": mu["source"],
            "sqrts": mu["sqrts"], "dim": mu["dim"],
            "type": mu["type"], "eff": mu["eff"], "l": cur_l})
	
    return likelihood_results, l

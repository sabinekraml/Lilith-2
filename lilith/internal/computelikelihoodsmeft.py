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
import os

def compute_likelihood_smeft(exp_mu, smeft_user_mu, smread, smcorr_read):
    """Computes the likelihood from experimental mu and user mu."""
    
    likelihood_results = []
    l = 0. # actually -2log(likelihood)

    for mu in exp_mu:
        # compute user mu value scaled to efficiencies
        user_mu_effscaled = {}
	
	# read SM prediction 
        stxs_th_bin = mu["SMpred"][:,0]
        error_th_m = mu["SMpred"][:,1]
        error_th_p = mu["SMpred"][:,2]
      
        # read SM correlation  
        corr_m_th = mu["SMcorr"] 
	
        try:
            if mu["dim"] >= 3:
                for i in range(1,mu["dim"]+1):
                    d = "d" + str(i)
                    user_mu_effscaled[d] = 0.
                    for (prod,decay),eff_prod in list(mu["eff"][d].items()):
                        if mu["dat"]==0: # signal strength run
                            user_mu_effscaled[d] += eff_prod*smeft_user_mu[i-1]
                        else:            # STXS run    
                            user_mu_effscaled[d] += eff_prod*smeft_user_mu[i-1]*stxs_th_bin[i-1]
            
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
                if mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],
                                       user_mu_effscaled["d2"] - mu["bestfit"]["d2"],
                                       user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])
                        cov_m = np.linalg.inv(mu["param"]["inv_cov_m"])
                    
                    ## include theoretical errors with no correlations:
                    unc_sym_th = np.reshape((error_th_p + error_th_m)/2.0,(1,mu["dim"]))
                    cov_m_th = unc_sym_th*corr_m_th*unc_sym_th.T
                    cov_m_tot = cov_m + cov_m_th    # with theo. correlation 
                    inv_cov_m = np.linalg.inv(cov_m_tot)
                    cur_l = inv_cov_m.dot(mu_vec).dot(mu_vec.T)

            # likelihood computation in case of a type="variable normal"
            # following "Variable Gaussian 2", Barlow arXiv:physics/0406120v1, Eq. 18
            if mu["type"] == "vn":
                if mu["dim"] >= 3:
                    mu_vec = np.array([user_mu_effscaled["d1"] - mu["bestfit"]["d1"],\
                    user_mu_effscaled["d2"] - mu["bestfit"]["d2"],\
                    user_mu_effscaled["d3"] - mu["bestfit"]["d3"]])
                    
                    for i in range(4,mu["dim"]+1):
                        d = "d"+str(i)
                        mu_vec = np.append(mu_vec,[user_mu_effscaled[d] - mu["bestfit"][d]])
                    unc_sym = np.sqrt(np.abs(mu["param"]["VGau"] + mu["param"]["VGau_prime"]*mu_vec))
                    cov_m = unc_sym*mu["param"]["corr_m"]*unc_sym.T
                    
                    ## include theoretical errors with no correlations:
                    mu_th_VGau = error_th_p*error_th_m
                    mu_th_VGau_prime = error_th_p - error_th_m
                    unc_sym_th = np.reshape(np.sqrt(np.abs(mu_th_VGau + mu_th_VGau_prime*mu_vec)),(-1,mu["dim"]))
                    cov_m_th = unc_sym_th*corr_m_th*unc_sym_th.T
                    cov_m_tot = cov_m + cov_m_th 

                    inv_cov_m = np.linalg.inv(cov_m_tot)
                    cur_l = inv_cov_m.dot(mu_vec).dot(mu_vec.T)

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
                    
                    ## include theoretical errors with no correlations:
                    mu_th_SGau = 2.0*error_th_p*error_th_m/(error_th_p + error_th_m)
                    mu_th_SGau_prime = (error_th_p - error_th_m)/(error_th_p + error_th_m)
                    unc_sym_th = np.reshape(np.abs(mu_th_SGau + mu_th_SGau_prime*mu_vec),(-1,mu["dim"]))         
                    cov_m_th = unc_sym_th*corr_m_th*unc_sym_th.T
                    cov_m_tot = cov_m + cov_m_th # with theo. correlation
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

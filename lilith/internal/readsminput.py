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

import sys
try:
    from lxml import etree
except:
    import xml.etree.ElementTree as etree
from ..errors import ExpInputError, ExpInputIOError
from scipy import interpolate
import numpy as np
from scipy.optimize import fsolve
import math
from warnings import warn
#from readexpinput import ReadExpInput

class ReadSMInput:
    """Read and extract the SM prediction values from TXT input."""

    def __init__(self,filename,exp_mu):
    
        self.readsmcount = 0 
        dimcount = exp_mu[0]["dim"]
        try:
            with open(filename) as readtest:
                print("reading SM predictions.............. ok")
                readtest.close()
                self.file = np.genfromtxt(filename)
                if len(self.file) == dimcount:
                    self.readsmcount = 1
                    print("test dim of SM prediction........... ok")
                else:
                    print("test dim of SM prediction........... error!")
                    print("SM prediction data does not have the same dim with experimental data.")    
        except IOError as readerror:
            print("reading SM predictions............... error!")
            print("SM prediction not found!")
        
        if self.readsmcount == 1:
            self.smpredic = np.array(self.file[:,:])
        else:
            self.smpredic = np.zeros((dimcount,3))    

#            self.smcent = np.array([np.zeros(mu["dim"])])
#np.array([SMbin_corr_data[:,:]])
 

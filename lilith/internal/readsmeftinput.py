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
from . import brsm as BR_SM
from warnings import warn
import os
import shutil
import re
from . import readexpinput as ReadEXPinput

class ReadSmeftInput:
    """Read the SMEFT parametrization and extract all information
    Including: 
        + SMEFT parametrization
        + STXS experimental data
    """
    
    def __init__(self):
        
        self.smeft_mu = []
        self.filepath = ""
        self.smeftfilelist = []
        self.rootlist = []
        self.fitparams = []
        
    def get_smeft_mu(self, filepath, fitparams):
        """Read SMEFT parametrization file list and extract user mu"""
        self.filepath = filepath
        # ~ print(self.filepath)
        self.fitparams = fitparams
        self.smeftfilelist = self.get_smeft_filelist(filepath)
        self.rootlist = self.get_rootlist()
        (stxs_source, br_source, acc_source) = self.get_smeft_source()
        self.usermu = self.get_ini_usermu(stxs_source, br_source, acc_source)
        (allparams, paramvalue,self.smeft_mu) = self.get_smeftparam_list(fitparams)
        temp = self.smeft_mu
        return temp
        
# -------------------
        
    def get_smeft_filelist(self,filepath):
        """get list of smeft-parametrization file from .list input"""
        smeft_files = []
        filepath_split = filepath.split("/")
        smeftdata_dir = "/".join(filepath_split[:-1])
        try:
            f = open(filepath,'r')
        except OSError:
            print("          can not open/read file list - ",filepath)
            sys.exit()
        with f: 
            for line in f:
                line = line.split("#")[0].rstrip("\n").strip()
                if line == "": # ignore empty line or comment
                    continue
                if line[0] == "/": # absolute filepath
                    try:
                        temp_reader = open(line,"r")
                        temp_reader.close()
                    except:
                        print("          can not read file '",line.strip(),"'")
                        sys.exit()
                    smeft_files.append(line)
                else: # relative filepath
                    temp_filename = smeftdata_dir+"/"+line
                    try:
                        temp_reader = open(temp_filename,"r")
                        temp_reader.close()
                    except OSError:
                        print("          can not read file '",temp_filename,"'")
                        sys.exit()
                    smeft_files.append(temp_filename)
        print("     Read SMEFT file list ".ljust(40,".")+" ok!".rjust(10,"."))
        print('data to be combined = '.rjust(35," "),len(smeft_files))
        return smeft_files
    
    def get_root(self,filepath):
        """Produce the XML tree with ElementTree"""
        tree = etree.parse(filepath)
        root = tree.getroot()
        return root    
    
    def get_rootlist(self):
        """get list of root from SMEFT XML files"""
        rlist = []
        for f in self.smeftfilelist:
            root = self.get_root(f)
            if root.tag == "smeft":
                rlist.append(root)
            else:
                print("     Check root-tags ".ljust(40,".")+" ERROR!".rjust(10,"."))
                print("          root tag is not correct : ",f)
                print("  >>>> quit process")
                sys.exit()
        print("     Check root-tags ".ljust(40,".")+" ok!".rjust(10,"."))
        return rlist
        
    def get_stxs_list(self):
        """get list of STXS experimental data filepaths and export as a list"""
        stxs_xml_list = "data/smeft-stxs.list"    
        f = open(stxs_xml_list,"w")
        for root in self.rootlist:
            # ~ idx = self.rootlist.index(root)
            for child in root:
                # ~ i = 0
                if child.tag == "stxsexpdata":
                    # ~ i = 1
                # ~ if i==1 and len(child.text)>0:
                    # ~ try:
                        # ~ stxspathcheck = open(child.text,"r")
                        # ~ stxspathcheck.close()
                    # ~ except OSError:
                        # ~ print("     Extract experimental stxs list ".ljust(40,".")+" ERROR!".rjust(10,"."))
                        # ~ print("          Can not open path to exp-data for ",self.smeftfilelist[idx])
                        # ~ print("  >>>> quit process ")
                        # ~ sys.exit()
                        # ~ print()
                    f.write(child.text+"\n")
                # ~ else:
                    # ~ print("     Extract experimental stxs list ".ljust(40,".")+" ERROR!".rjust(10,"."))
                    # ~ print("          can not find exp-data for ",self.smeftfilelist[idx])
                    # ~ print("  >>>> quit process ")
                    # ~ sys.exit()
        f.close()
        print("     Extract experimental stxs list ".ljust(40,".")+" ok!".rjust(10,"."))
        return stxs_xml_list
        
        
    def get_clean(self):
        """remove stxs list after calculate"""
        stxs_xml_list = "data/smeft-stxs.list"
        os.remove(stxs_xml_list)
        
    def get_smeft_source(self):
        """extract initial parametrized STXS bins from sources"""   
        stxs_source = {}
        br_source = {}
        acc_source = {}
        root_counter = -1
        for root in self.rootlist:
            root_counter += 1
            stxs_source[root_counter] = {}
            br_source[root_counter] = {}
            acc_source[root_counter] = {}
            for child in root:
                if child.tag == "smeftsource":
                    try:
                        filepath = self.filepath
                        filepath_split = filepath.split("/")
                        filepath = "/".join(filepath_split[:-2])+child.text
                        fread = open(filepath,"r")
                        lines = fread.readlines()
                        line_count = 0 
                        for line in lines:
                            line_count +=1
                            if line.startswith("#") or line.strip()=="":
                                continue
                            if line.startswith("bin number"):
                                stxs_num = line.split()[2]
                                stxs_text = lines[line_count]
                                stxs_source[root_counter][stxs_num] = stxs_text.strip()
                            elif line.startswith("br number"):
                                br_num = line.split()[2]
                                br_text = lines[line_count]
                                br_source[root_counter][br_num] = br_text.strip()
                            elif line.startswith("acc"):
                                acc_num = line.split()[2]
                                acc_text = lines[line_count]
                                acc_source[root_counter][acc_num] = acc_text.strip()
                    except OSError:
                        print("     Read SMEFT parametrization source ".ljust(40,".")+" ERROR!".rjust(10,"."))
                        print("          can not read source file "+filepath)
                        print("  >>>> quit process")
                        sys.exit()
        print("     Read SMEFT parametrization source ".ljust(40,".")+" ok!".rjust(10,"."))
        return (stxs_source,br_source,acc_source)
         
    def get_ini_usermu(self, ssource, bsource, asource):
        """get initial parametrized expression of user mu"""
        usermu = {}
        root_counter = -1
        sourceref = ""
        for root in self.rootlist:
            root_counter += 1
            usermu[root_counter] = {}
            for child in root:
                if child.tag == "smeftSMsource":
                    sourceref = child.text
                if child.tag == "smeftbin":
                    usermu[root_counter]["dim"] = child.attrib["bindim"]
                    for gchild in child:
                        pro = ""
                        dec = ""
                        acc  = ""
                        for ggchild in gchild:
                            if ggchild.tag == "pro":
                                ggchild_split = ggchild.text.split('+')
                                if len(ggchild_split)==1:
                                    pro = ssource[root_counter][ggchild.text]
                                else:
                                    print("     Warning: combine bin at ",gchild.tag)
                                    try:
                                        f = open(sourceref,'r')
                                        smstxs = {}
                                        for line in f:
                                            line = line.partition('#')[0].rstrip("\n")
                                            line = line.split()
                                            if len(line)==2:
                                                smstxs[line[0]] = float(line[1])
                                        total = 0
                                        pro = ''
                                        for i in ggchild_split:
                                            total += smstxs[i]
                                        for i in ggchild_split:
                                            r = str(smstxs[i]/total)
                                            # ~ cbratio.append(r)
                                            pro += " + "+r+" *( "+ssource[root_counter][i]+" )"                                    
                                    except OSError:
                                        print("     Read parameterized bins ".ljust(40,".")+" ERROR!".rjust(10,"."))
                                        print("          No SM prediction for STXS found")
                                        print("          >>>> bins combined with the same ratio")
                                        combinenum = len(ggchild_split)
                                        cbratio = str(1/combinenum)
                                        for i in ggchild_split:
                                            pro += " + "+cbratio+" *( "+ssource[root_counter][i]+" )"
                            if ggchild.tag == "dec":
                                dec = bsource[root_counter][ggchild.text]
                            if ggchild.tag == "acc":
                                acc = asource[root_counter][ggchild.text]
                        usermu[root_counter][gchild.tag] = "("+pro+")*("+dec+")*("+acc+")"                                    
        print("     Read parameterized bins ".ljust(40,".")+" ok!".rjust(10,"."))
        return usermu
            
    def get_smeftparam_list(self,fitlist):
        bin_list = []
        usermu = self.usermu
        for outkey in usermu:
            del usermu[outkey]["dim"]
            bin_list = bin_list+list(self.usermu[outkey].values())
        params_list = []    # full list of parameters used 
        fitlist = list(self.fitparams.keys())
        smeft_mu = []
        
        for bitem in bin_list:
            params_list += re.findall('c\w+',bitem)
        # ~ print("params_list 0 = ",params_list)
        params_list = list( dict.fromkeys(params_list) )
        # ~ print("params_list 1 = ",params_list)
        params_list = sorted(params_list, key=len, reverse=True)
        # ~ print("params_list 2 = ",params_list)
        param_values = {}
        
        for param in fitlist:
            if param in params_list:
                pass
            else:
                print("     List used SMEFT parameters  ".ljust(40,".")+" ERROR!".rjust(10,"."))
                print("          fitting parameter [",param,"] not in parameter list")
                print("          available parameters are ",params_list)
                print("  >>>> quit process")
                sys.exit()
        for param in params_list:
            if param not in fitlist:
                param_values[param] = '0'
                # ~ print(param," --> ", 0)
            else: 
                # ~ params_value[param] = 'pr'+str(fitlist.index(param))
                param_values[param] = self.fitparams[param]
                # ~ print(param," --> ",param_values[param])
        # ~ print("params_value = ", param_values)
        for outkey in usermu:
            temp_mu = []
            for inkey in usermu[outkey]:
                temp_comp = usermu[outkey][inkey]
                for param in param_values.keys():
                    temp_comp = temp_comp.replace(param,param_values[param])
                    # ~ print(" replace ",param)
                    # ~ print(temp_comp)
                temp_mu.append(temp_comp)
            smeft_mu.append(temp_mu)
        print("     List parameters of interest (POI) ".ljust(40,".")+" ok!".rjust(10,"."))
        print('total parameters = '.rjust(35," "),len(params_list))
        print('POIs = '.rjust(35," "),len(fitlist))
        return (params_list,param_values,smeft_mu)
    

             
        

#######################################################################
# This file is for testing of reading data from expinput  
#######################################################################

from lxml import etree

import numpy as np
import math 

with open("HIGG-2018-28_stxs_vn_dim12-4test.xml","r") as fread:
	tree = etree.parse(fread)

root = tree.getroot()	
#print(root.tag)	
#print(root.attrib)
for child in root:
	if child.tag == "eff":
		eff = child.attrib
		if eff["axis"] == "d3":
#			print(eff["brratio"])
			try:
				print(bool(int(eff["brratio"])))
			except IOError:
				print("no br_ratio")	
#			if int(eff["brratio"]) == 1:
#				print("ok, good")
#				print(bool(int(eff["brratio"])))
			

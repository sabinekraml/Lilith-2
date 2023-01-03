import numpy as np
import matplotlib.pyplot as plt
import matplotlib

f = open("smeftcoef.txt","r")
#for x in f:
#   print(x)
    
cHB = 0
cHWB = 0
cHG = 0
cuH = 0

s = []
for x in f:
    x = x.replace(" ","")
    if not x.startswith("#"):
	    s.append(x)
	
f.close()

s=list(map(str.strip,s))

news = []
i = j = 0
for i in range(1,len(s)):
    if s[i] == "":
        if i-j>1:
            news.append("".join(s[j+1:i]))
        j = i 
        
xpoint = []
y1point = []
y2point = []
y3point = []
for cHW in np.linspace(-6,6,100):
    xpoint.append(cHW)
    y1point.append(eval(news[9]))
    y2point.append(eval(news[9])*eval(news[-2]))
    y3point.append(eval(news[9])*eval(news[-2])*eval(news[-1]))

plt.plot(xpoint,y1point, c="b")
plt.plot(xpoint,y2point,c="r")
plt.plot(xpoint,y3point,c="k")
plt.show()
#plt.figure().savefig("testplot.pdf")



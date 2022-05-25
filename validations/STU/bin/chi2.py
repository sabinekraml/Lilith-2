from scipy import misc 
from scipy import stats

PValueList = [0.9, 0.5, 0.3, 0.32, 0.2, 0.1, 0.05]

global pvalue, dfreedom

def newtons_method(f, x, tolerance=0.0001):
    while True:
        x1 = x - f(x) / misc.derivative(f, x) 
        t = abs(x1 - x)
        if t < tolerance:
            break
        x = x1
    return x

def f(x):
    return 1 - stats.chi2.cdf(x, dfreedom) - pvalue

print(  'df\p' , '| ', PValueList[0], ' | ', PValueList[1], ' | ', PValueList[2], ' | ', \
                      PValueList[3], ' | ', PValueList[4], ' | ', PValueList[5], ' | ', \
                      PValueList[6], ' | ')

for i in range(2,10):
    dfreedom = i 
    Result = []
    for pvalue in PValueList:
        x0 = dfreedom  # x0 approximation
        x = newtons_method(f, x0)
        Result.append(x)
    for i in range(7):
        Result[i] = round(Result[i],3)
    print( dfreedom, ' | ', Result[0], ' | ', Result[1], ' | ', Result[2], ' | ', Result[3], ' | ', \
          Result[4], ' | ', Result[5], ' | ', Result[6], ' | ' )

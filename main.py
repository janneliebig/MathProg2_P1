import matplotlib.pyplot as plt 
import numpy as np
import math

from integrationrules import IntegrationRules
from integrationrules import Function

def plot_Error_M(f ,m):    #b)
    M = np.linspace(1, m, m)
    milne_rule_error = abs(IntegrationRules.milne_rule_error(f, 2**M, 0))
    trapezoid_rule_error = abs(IntegrationRules.trapezoid_rule_error(f, 2**M, 0))
    mid_point_rule_error = abs(IntegrationRules.mid_point_rule_error(f, 2**M, 0))
    simpson_rule_error = abs(IntegrationRules.simpson_rule_error(f, 2**M, 0))

    plt.plot(M, milne_rule_error)
    plt.plot(M, trapezoid_rule_error)
    plt.plot(M, mid_point_rule_error)
    plt.plot(M, simpson_rule_error)

    plt.legend()
    plt.show()    

def plot_Error_n(f, M):   #c)
    M=M+1
    result=np.zeros((M,4))

    for M in range(0,M):
        result[M][3] = abs(IntegrationRules.milne_rule_error(f, 2*2**M, 0))
        result[M][1] = abs(IntegrationRules.trapezoid_rule_error(f, 2*2**M, 0))
        result[M][0] = abs(IntegrationRules.mid_point_rule_error(f, 2*2**M, 0))
        result[M][2] = abs(IntegrationRules.simpson_rule_error(f, 2*2**M, 0))

    for M in range(0,M):
        plt.plot(range(1,5),result[M])
    plt.legend()
    plt.show()    

f = Function(-1, 1)

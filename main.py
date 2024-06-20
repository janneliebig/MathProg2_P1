import matplotlib.pyplot as plt 
import numpy as np
import math

from integrationrules import IntegrationRules
from integrationrules import Function

def plot_Error_M(f: Function,m: int):    #b)
    """
    Plots error depending on number of intervals.

    Args:   
        f: Function -> Function
        m: Number of intervals -> int
    """
    # begin with 0. 2**0 = 1
    M = np.linspace(0, m, m)
    milne_rule_error = abs(IntegrationRules.milne_rule_error(f, 2**M, 0))
    trapezoid_rule_error = abs(IntegrationRules.trapezoid_rule_error(f, 2**M, 0))
    mid_point_rule_error = abs(IntegrationRules.mid_point_rule_error(f, 2**M, 0))
    simpson_rule_error = abs(IntegrationRules.simpson_rule_error(f, 2**M, 0))

    plt.plot(M, milne_rule_error)
    plt.plot(M, trapezoid_rule_error)
    plt.plot(M, mid_point_rule_error)
    plt.plot(M, simpson_rule_error)

    plt.title("Integration error for M = 2$^j$, $j \in \{0, \dots, 8\}$ and $xi = 0.$")

    plt.xscale("log") # optional
    plt.xlabel("function eval.")
    plt.ylabel("intervals")

    plt.legend(["milne", "trapezoid", "mid point", "simpson"], loc="upper right")
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

def plot_Error_n_2(f: Function, start: int, end: int): # c) Andere Art der Darstellung.
    """
    Plots integration error of each rule depending on number of function evaluations.

    Args:
        f: Function -> Function
        start: Minimal number of function calls -> int
        end: Maximal number of function calls -> int
    """
    m = np.linspace(1, end, end)
    milne_rule_error = abs(IntegrationRules.milne_rule_error(f, 5*m, 0))
    trapezoid_rule_error = abs(IntegrationRules.trapezoid_rule_error(f, m*2, 0))
    mid_point_rule_error = abs(IntegrationRules.mid_point_rule_error(f,m , 0))
    simpson_rule_error = abs(IntegrationRules.simpson_rule_error(f, m*3, 0))

    # plot
    plt.plot(5*m, milne_rule_error)
    plt.plot(2*m, trapezoid_rule_error)
    plt.plot(m, mid_point_rule_error)
    plt.plot(3*m, simpson_rule_error)
    
    plt.xlim(right=end, left=start)

    plt.xlabel("func. calls")
    plt.ylabel("error")
    
    plt.title("Integration error depending on number of function evaluations")

    plt.legend(["milne", "trapezoid", "mid point", "simpson"], loc="upper right")
    plt.show()    



f = Function(-1, 1)

plot_Error_M(f, 8)
plot_Error_n_2(f,1, 15)
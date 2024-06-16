import matplotlib.pyplot as plt 
import numpy as np
import math

from integrationrules import IntegrationRules
from integrationrules import Function

f = Function(-1, 1)

M = 8

x = np.linspace(-1, 1, 1000)
milne_rule_error = abs(IntegrationRules.milne_rule_error(f, M, x))
trapezoid_rule_error = abs(IntegrationRules.trapezoid_rule_error(f, M, x))
mid_point_rule_error = abs(IntegrationRules.mid_point_rule_error(f, M, x))
simpson_rule_error = abs(IntegrationRules.simpson_rule_error(f, M, x))

plt.plot(x, milne_rule_error)
plt.plot(x, trapezoid_rule_error)
plt.plot(x, mid_point_rule_error)
plt.plot(x, simpson_rule_error)
plt.legend()
plt.show()
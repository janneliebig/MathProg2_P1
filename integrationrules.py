import numpy as np

class Function:

    def __init__(self, a: float, b: float):
        self.interval = (a,b)

    def at(self, x: float):
        return np.exp(x)

    def derivative_at(self, n: int, x: float):
        return np.exp(x)    
    

class IntegrationRules:
    def mid_point_rule(f: Function, M: int):
        n = 0
        h = (f.interval[1] - f.interval[0])/M
        result: float = 0.0

        for i in range(1, M + 1):
            x_k : float = f.interval[0] +  ((2*i - 1)*h)/2
            result += f.at(x_k)

        return h * result

    def mid_point_rule_error(f : Function, M: int, xi: float):
        # xi must be in (a,b)
        h : float = (f.interval[1] - f.interval[0])/M
        error : float = (f.interval[1] - f.interval[0])/24 * h**2 * f.derivativeat(2, xi)
        return error
    
    def trapezoid_rule(f: Function, M: int):
        n = 1
        b : float = f.interval[1]
        a : float = f.interval[0]

        h = (b - a) / M
        result : float = 0.0

        result += f.at(a)

        for i in range(1, M):
            x_k : float = a + i*h
            result += 2 * f.at(x_k)

        result += f.at(b)
        return (h/2) * result
    
    def simpson_rule(f : Function, M: int):
        n = 2
        b : float = f.interval[1]
        a : float = f.interval[0]

        h = (b - a) / M

        area : float = f.at(a)

        for i in range(1, M):
            x_2k : float = a + (2*i*h)/2
            area += 2*f.at(x_2k)

        for i in range(0, M):
            x_2k1 : float = a + ((2*i + 1)*h / 2)
            area += 4 *f.at(x_2k1)

        area += f.at(b)
        return (h/6)*area
    
    def simpson_rule_error(f : Function, M : int, xi: float):
        # xi must be in (a,b)
        a : float = f.interval[0]
        b : float = f.interval[1]

        h : float = (b-a)/M

        error : float = -(b-a)/180 * (h/2)**4 * f.derivative_at(4, xi)
        return error

f = Function(-1, 1)
print("int_a^b exp(x) dx: " + str(np.exp(1) - np.exp(-1)))
print("mid point rule:    " + str(IntegrationRules.mid_point_rule(f, 100)))
print("trapezoid rule:    " + str(IntegrationRules.trapezoid_rule(f, 100)))
print("simpson rule:      " + str(IntegrationRules.simpson_rule(f, 100)))

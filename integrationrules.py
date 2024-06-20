import numpy as np

class Function:

    def __init__(self, a: float, b: float):
        self.interval = (a,b)

    def at(self, x: float) -> float:
        """
        Returns value of function in x.

        Args:
            x: point -> float
        Returns:
            Returns exp(x) as default.
        """
        return np.exp(x)

    def derivative_at(self, n: int, x: float) -> float:
        """
        Returns derivative of function in x.

        Args:
            x: point -> float
        Returns:
            Returns exp(x) as default as exp is its derivative.
        """
        return np.exp(x)    
    

class IntegrationRules:
    def mid_point_rule(f: Function, M: int) -> float:
        """
        Returns the definite numeric integral of a function using the mid point rule.

        Args:
            f: Function -> Function
            M: number of intervals -> int
        Returns:
            Approximate definite integral of f using mid point rule.
        """
        b : float = f.interval[1]
        a : float = f.interval[0]
        h = (b - a)/M
        result: float = 0.0

        for i in range(1, M + 1):
            x_k : float = a +  ((2*i - 1)*h)/2
            result += f.at(x_k)

        return h * result

    def mid_point_rule_error(f : Function, M: int, xi: float) -> float:
        """
        Returns the numeric error of mid point rule.

        Args:
            f: Function -> Function
            M: Number of intervals -> int
            xi: Point in interval
        
        Returns:
            Numeric error between integral and its approximation using the mid point rule. 
        """
        b : float = f.interval[1]
        a : float = f.interval[0]

        if xi <= a or xi >= b:
            raise Exception("Xi must be between a and b.")

        h : float = (b - a)/M
        error : float = (b - a)/24 * h**2 * f.derivative_at(2, xi)
        return error
    
    def trapezoid_rule(f: Function, M: int) -> float:
        """
        Returns the definite numeric integral of a function using the trapezoid rule.

        Args:
            f: Function -> Function
            M: number of intervals -> int
        Returns:
            Approximate definite integral of f using trapezoid rule.
        """
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
    
    def trapezoid_rule_error(f: Function, M: int, xi : float) -> float:
        """
        Returns the numeric error of trapezoid rule.

        Args:
            f: Function -> Function
            M: Number of intervals -> int
            xi: Point in interval
        
        Returns:
            Numeric error between integral and its approximation using the trapezoid rule. 
        """
        b : float = f.interval[1]
        a : float = f.interval[0]
        h = (b-a)/M

        if xi <= a or xi >= b:
            raise Exception("Xi must be between a and b.")

        error = -(b-a)/12 * (h**2) * f.derivative_at(2, xi)
        return error
    
    def simpson_rule(f : Function, M: int) -> float:
        """
        Returns the definite numeric integral of a function using the simpson rule.

        Args:
            f: Function -> Function
            M: number of intervals -> int
        Returns:
            Approximate definite integral of f using the simpson rule. 
        """
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
    
    def simpson_rule_error(f : Function, M : int, xi: float) -> float:
        """
        Returns the numeric error of simpson rule.

        Args:
            f: Function -> Function
            M: Number of intervals -> int
            xi: Point in interval
        
        Returns:
            Numeric error between integral and its approximation using the simpson rule. 
        """
        a : float = f.interval[0]
        b : float = f.interval[1]

        if xi <= a or xi >= b:
            raise Exception("Xi must be between a and b.")

        h : float = (b-a)/M

        error : float = -(b-a)/180 * (h/2)**4 * f.derivative_at(4, xi)
        return error
    
    def milne_rule(f: Function, M : int) -> float:
        """
        Returns the definite numeric integral of a function using the milne rule.

        Args:
            f: Function -> Function
            M: number of intervals -> int
        Returns:
            Approximate definite integral of f using the milne rule. 
        """
        a : float = f.interval[0]
        b : float = f.interval[1]
        h = (b - a)/(M)
        area : float = 0.0
        for i in range(0, M):
            area += 7*f.at(a + i*h) + 32 * f.at(a + (i+ 1/4)*h) + 12 * f.at(a + (i+1/2)*h) + 32 * f.at(a + (i+3/4)*h) + 7 * f.at(a + (i+1)*h)

        return (b-a)/(90 * M) * area
    
    def milne_rule_error(f: Function, M: int, xi: float) -> float:
        """
        Returns the numeric error of milne rule.

        Args:
            f: Function -> Function
            M: Number of intervals -> int
            xi: Point in interval
        
        Returns:
            Numeric error between integral and its approximation using the milne rule. 
        """
        a : float = f.interval[0]
        b : float = f.interval[1]
        
        if xi <= a or xi >= b:
            raise Exception("Xi must be between a and b.")
        
        h = (b - a)/M

        error = (-2) * (b-a)/945 * ((h/4)**6) * f.derivative_at(6, xi)

        return error        

f = Function(-1, 1)
m = 100
print("int_a^b exp(x) dx: " + str(np.exp(1) - np.exp(-1)))
print("mid point rule:    " + str(IntegrationRules.mid_point_rule(f, m)))
print("mid point error:   " + str(IntegrationRules.mid_point_rule_error(f, m, 0)))
print("trapezoid rule:    " + str(IntegrationRules.trapezoid_rule(f, m)))
print("trapezoid error:   " + str(IntegrationRules.trapezoid_rule_error(f, m, 0)))
print("simpson rule:      " + str(IntegrationRules.simpson_rule(f, m)))
print("simpson error:     " + str(IntegrationRules.simpson_rule_error(f, m, 0)))
print("milne rule:        " + str(IntegrationRules.milne_rule(f, m)))
print("milne error:       " + str(IntegrationRules.milne_rule_error(f, m, 0)))

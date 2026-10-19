def scalePolynomial(poly: list, scaler):
    for i in range(len(poly)):
        poly[i] *= scaler

def multiplyPolynomial(x: list, y: list):
    n = len(x)
    m = len(y)
    poly: list = [0] * (n + m - 1)
    for i in range(n):
        for j in range(m):
            poly[i+j] += x[i]*y[j]
            
    return poly

def dividePolynomial(x: list, y: list):
    n = len(x)
    m = len(y)
    if n < m:
        return [0]                 # divisor has higher degree: quotient is 0

    rem = x[:]                     # copy, so we don't mutate the input
    quot = [0] * (n - m + 1)

    for k in range(n - m, -1, -1):
        quot[k] = rem[k + m - 1] / y[-1]
        for j in range(m):
            rem[k + j] -= quot[k] * y[j]

    return quot                    # list only, remainder dropped

def multiplyNum(x: list):
    n = len(x)
    poly = [1]
    for i in range(n):
        poly = multiplyPolynomial(poly, [-1*x[i], 1])
    
    return poly
    
def lagrange(x: list, y: list):
    n = len(x)
    m = len(y)
    if(n != m):
        # error
        return
    
    fullNumerator = multiplyNum(x)
     
    numeratorProduct = [0] * n
    denumeratorProduct = [0] * n
        
    for i in range(n):
        numeratorProduct[i] = dividePolynomial(fullNumerator, [-1*x[i], 1])
        pro = 1
        for j in range(n):
            if(j!=i):
                pro *= x[i] - x[j]
        denumeratorProduct[i] = pro
        
    result = [0] * n
        
    for i in range(n):
        scalePolynomial(numeratorProduct[i], y[i]/denumeratorProduct[i])
        for k in range(n):
            result[k] += numeratorProduct[i][k]
    
    return result

def evalPoly(coeffs: list, t):
    result = 0
    for c in reversed(coeffs):
        result = result * t + c
    return result

def closeLists(a: list, b: list, tol=1e-9):
    return len(a) == len(b) and all(abs(p - q) < tol for p, q in zip(a, b))

def testLagrange():
    # (xs, ys, expected coefficients, lowest power first)
    cases = [
        ([2],            [5],            [5]),                  # constant
        ([0, 1],         [1, 3],         [1, 2]),               # 1 + 2x
        ([0, 1, 2],      [1, 3, 2],      [1, 3.5, -1.5]),       # the example from earlier
        ([-1, 0, 1],     [1, 0, 1],      [0, 0, 1]),            # x^2
        ([1, 2, 3],      [2, 3, 4],      [1, 1, 0]),            # collinear: top coeff is 0
        ([0, 1, 2, 3],   [0, 1, 8, 27],  [0, 0, 0, 1]),         # x^3
    ]

    for xs, ys, expected in cases:
        result = lagrange(xs, ys)
        assert closeLists(result, expected), f"xs={xs}, ys={ys}: got {result}, expected {expected}"

    # Property 1: the polynomial passes through every given point
    xs = [-2, 0, 1, 3, 4]
    ys = [5, -1, 2, 0.5, 7]
    poly = lagrange(xs, ys)
    assert len(poly) == len(xs)            # degree <= n-1
    for xi, yi in zip(xs, ys):
        assert abs(evalPoly(poly, xi) - yi) < 1e-9, f"misses point ({xi}, {yi})"

    # Property 2: point order must not change the result
    order = [3, 0, 4, 2, 1]
    poly2 = lagrange([xs[i] for i in order], [ys[i] for i in order])
    assert closeLists(poly, poly2), "result depends on point order"

    print("lagrange: all tests passed")

def main():
    testLagrange()
    pass

if __name__ == "__main__":
    main()
from math import comb, ceil

def uniform_distribution(eta):
    """
    Probability density function of the uniform distribution for -eta to eta 
    :param eta: integer
    :returns: A dictionnary {x:p(x) for x in {-eta .. eta}}
    """

    D = {}
    for i in range(-eta, eta + 1):
        D[i] = 1 / (2*eta + 1)

    return D

def centered_binomial_distribution(eta):
    """
    Probability density function of the centered binomial distribution for -eta to eta
    :param eta: integer
    :returns: A dictionary {x:p(x) for x in {-eta .. eta}}
    """
    
    return {k: comb(2*eta, eta +k) * (0.5)**(2*eta) for k in range(-eta, eta+1)}

def convolution_law(A, B):
    """
    Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionary)
    :param B: second input law (dictionary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a + b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def clean_dist(A):
    """ 
    Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B


def iter_convolution_law(A, i):
    """ 
    compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = convolution_law(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = convolution_law(D, A)
            D = clean_dist(D)
    return D


def product_law(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionary)
    :param B: second input law (dictionary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = abs(max(D.keys()))

    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s
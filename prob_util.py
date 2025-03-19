from math import comb

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

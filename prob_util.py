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
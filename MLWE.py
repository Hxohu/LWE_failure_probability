import prob_util
from math import log, comb

def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q

def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D

def centered_binomial_distribution(eta):
    """
    Probability density function of the centered binomial distribution for -eta to eta
    :param eta: integer
    :returns: A dictionary {x:p(x) for x in {-eta .. eta}}
    """
    return {k: comb(2*eta, eta +k) * (0.5)**(2*eta) for k in range(-eta, eta+1)}

class Kyber_768():
    q = 3329
    n = 256
    k = 3
    eta_1 = 2
    eta_2 = 2
    d_u = 10
    d_v = 4

    e = None
    r = None
    s = None
    e_1 = None
    e_2 = None
    compress_u = None
    compress_v = None
    compress_t = None

    def __init__(self):
        self.e = prob_util.centered_binomial_distribution(self.eta_1)
        self.r = prob_util.centered_binomial_distribution(self.eta_1)
        self.s = prob_util.centered_binomial_distribution(self.eta_1)
        self.e_1 = prob_util.centered_binomial_distribution(self.eta_2)
        self.e_2 = prob_util.centered_binomial_distribution(self.eta_2)
        self.compress_u = build_mod_switching_error_law(self.q, 2**10)
        self.compress_v = build_mod_switching_error_law(self.q, 2**4)
        self.compress_t = build_mod_switching_error_law(self.q, 2**12)

    def failure_prob(self):
        compress_s = prob_util.convolution_law(self.s, self.compress_t)
        compress_s = prob_util.product_law(compress_s, self.e_1)
        B1 = prob_util.product_law(self.r, compress_s)
        C1 = prob_util.iter_convolution_law(B1, self.n * self.k)

        compress_e = prob_util.convolution_law(self.e, self.compress_u)
        B2 = prob_util.product_law(self.s, compress_e)
        C2 = prob_util.iter_convolution_law(B2, self.n * self.k)

        C = prob_util.convolution_law(C1, C2)
        F = prob_util.convolution_law(self.compress_v, self.e_2)
        D = prob_util.convolution_law(C, F)

        return prob_util.tail_probability(D, self.q / 4)
    
if __name__ == "__main__":
    kyber_768 = Kyber_768()
    # print(kyber_768.compress_v)
    f = kyber_768.failure_prob()
    print ("failure: %.1f = 2^%.1f"%(f, log(f + 2.**(-300))/log(2)))
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

def convolution_remove_dependency(A, B, q, p):
    normalizer = {}
    maxa = q / p
    for a in A:
        normalizer[a % maxa] = normalizer.get(a % maxa, 0) + A[a]

    C = {}
    for a in A:
        for b in B:
            c = a - b
            if (c % maxa == 0):
                C[c] = C.get(c, 0) + A[a] * B[b] / normalizer[a % maxa]
    
    return C

class Saber():
    l = 3
    n = 256
    q = 2**13
    p = 2**10
    T = 2**4
    mu = 8 // 2

    def failure_prob(self):
        s = prob_util.centered_binomial_distribution(self.mu)    # pdf secret
        e = build_mod_switching_error_law(self.q, self.p)    	# pdf error

        se = prob_util.product_law(s, e)
        se2 = prob_util.iter_convolution_law(se, self.n * self.l)
        se2 = convolution_remove_dependency(se2, se2, self.q, self.p)

        e2 = build_mod_switching_error_law(self.q, self.T)
        
        D = prob_util.convolution_law(se2, e2)
        
        pro = prob_util.tail_probability(D, self.q/4.)
        return pro
    
if __name__ == "__main__":
    saber = Saber()
    f = saber.failure_prob()
    print ("failure: %.1f = 2^%.1f"%(f, log(f + 2.**(-300))/log(2)))
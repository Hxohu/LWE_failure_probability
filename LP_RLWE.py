import prob_util
from math import log, comb

def LAC_phi_2():
    D = {1 : 1/8, 0 : 3/4, -1 : 1/8}
    return D

class LAC_192():
    n = None
    q = None
    l_t = 8   # can correct 8 errors
    l_n = 256 + 72  # total bit nubmber

    s_Gen = None
    s_Enc = None
    e_Gen = None
    e_Enc = None
    e_Enc_ = None

    def __init__(self, n, q):
        self.n = n
        self.q = q

    def select_distribution(self):
        self.s_Gen = LAC_phi_2()
        self.s_Enc = LAC_phi_2()
        self.e_Gen = LAC_phi_2()
        self.e_Enc = LAC_phi_2()
        self.e_Enc_ = LAC_phi_2()

    def single_bit_failure_prob(self):
        self.select_distribution()

        es_ = prob_util.product_law(self.e_Gen, self.s_Enc)
        e_s = prob_util.product_law(self.e_Enc, self.s_Gen)

        es_ = prob_util.iter_convolution_law(es_, self.n)
        e_s = prob_util.iter_convolution_law(e_s, self.n)

        D = prob_util.convolution_law(e_s, es_)
        F = prob_util.convolution_law(D, self.e_Enc_)

        return prob_util.tail_probability(F, int(self.q / 4))
    
    def dec_failure(self):
        bit_failure = self.single_bit_failure_prob()
        bit_failure_ = 1 - bit_failure
        pro = 0.0

        for i in range(self.l_t + 1, self.l_n + 1):
            pro += comb(self.l_n, i) * (bit_failure**i) * (bit_failure_ ** (self.l_n - i))

        return pro
    
if __name__ == "__main__":
    lac_192 = LAC_192(1024, 251)

    f = lac_192.single_bit_failure_prob()
    print ("failure: %.2f = 2^%.2f"%(f, log(f + 2.**(-300))/log(2)))

    dec_f = lac_192.dec_failure()
    print ("failure: %.2f = 2^%.2f"%(f, log(dec_f + 2.**(-300))/log(2)))
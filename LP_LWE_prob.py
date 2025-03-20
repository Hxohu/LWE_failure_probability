import prob_util
from math import log

class LP_LWE():
    q = None
    n1 = None
    n2 = None
    m_len = None

    R1 = None
    R2 = None
    e1 = None
    e2 = None
    e3 = None

    def __init__(self, q, n1, n2, m_len):
        self.q = q
        self.n1 = n1
        self.n2 = n2
        self.m_len = m_len

    def select_distribution(self):
        self.R1 = prob_util.centered_binomial_distribution(int(input("请输入R1的中心二项分布概率参数:")))
        self.R2 = prob_util.centered_binomial_distribution(int(input("请输入R2的中心二项分布概率参数:")))
        self.e1 = prob_util.centered_binomial_distribution(int(input("请输入e1的中心二项分布概率参数:")))
        self.e2 = prob_util.centered_binomial_distribution(int(input("请输入e2的中心二项分布概率参数:")))
        self.e3 = prob_util.centered_binomial_distribution(int(input("请输入e3的中心二项分布概率参数:")))

    def failure_prob(self):
        self.select_distribution()

        e2R2 = prob_util.product_law(self.e2, self.R2)
        e1R1 = prob_util.product_law(self.e1, self.R1)

        e2R2 = prob_util.iter_convolution_law(e2R2, self.n2)
        e1R1 = prob_util.iter_convolution_law(e1R1, self.n1)

        D = prob_util.convolution_law(e1R1, e2R2)
        F = prob_util.convolution_law(D, self.e3)

        return prob_util.tail_probability(F, int(self.q / 4)) * self.m_len
    
if __name__ == "__main__":
    lp_lwe = LP_LWE(3329, 640, 640, 128)
    f = lp_lwe.failure_prob()

    print ("failure: %.1f = 2^%.1f"%(f, log(f + 2.**(-300))/log(2)))

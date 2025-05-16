import prob_util
from math import log

def frodo1344_dis():
    pro = [18286, 14320, 6876, 2023, 364, 40, 2]
    D = {}
    for i in range(-6, 7):
        D[i] = pro[abs(i)] / (2**16)
    
    return D

class Frodo_1344():
    q = None
    n1 = None
    n2 = None

    R1 = None
    R2 = None
    e1 = None
    e2 = None
    e3 = None

    def __init__(self, q, n1, n2):
        self.q = q
        self.n1 = n1
        self.n2 = n2
        
    def select_distribution(self):
        self.R1 = frodo1344_dis()
        self.R2 = frodo1344_dis()
        self.e1 = frodo1344_dis()
        self.e2 = frodo1344_dis()
        self.e3 = frodo1344_dis()

    def failure_prob(self):
        self.select_distribution()

        e2R2 = prob_util.product_law(self.e2, self.R2)
        e1R1 = prob_util.product_law(self.e1, self.R1)

        e2R2 = prob_util.iter_convolution_law(e2R2, self.n2)
        e1R1 = prob_util.iter_convolution_law(e1R1, self.n1)

        D = prob_util.convolution_law(e1R1, e2R2)
        F = prob_util.convolution_law(D, self.e3)

        return prob_util.tail_probability(F, int(self.q / 32)) 
    
if __name__ == "__main__":
    frodo_1344 = Frodo_1344(2**16, 1344, 1344)
    f = frodo_1344.failure_prob()

    print ("failure: %.1f = 2^%.1f"%(f, log(f + 2.**(-300))/log(2)))

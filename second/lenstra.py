from elliptic_curve import EllipticCurve
from point import Point
from ideal import Ideal
from gmpy2 import invert
from gmpy2 import f_mod
from gmpy2 import mpz_random
from random import randint
from gmpy2 import random_state
from gmpy2 import gcd
from gmpy2 import fac

from gmpy2 import is_prime
from gmpy2 import next_prime
import math
import copy

class Lenstra(object):
    def create_curve(self):
        """
        :return: random Elliptic curve  and init Point
        """
        for i in range(10):
            try:
                self.X = self._generate_random_number()
                self.Y = self._generate_random_number()
                self.A = self._generate_random_number()
                self.B = f_mod((self.Y * self.Y - self.X * self.X * self.X - self.A * self.X), self.N)
                curve = EllipticCurve(self.N, self.A, self.B)
                return curve
            except Exception as e:
                print("Error when curve is generated %s" % str(e))
                raise Exception

    def __init__(self, N, t=0):
        self.N = N
        self.t = t if t != 0 else None
        self.curve = self.create_curve()
        self.point = Point(self.curve, self.X, self.Y)

    def _generate_random_number(self):
        """
        :return: random number  
        """
        rstate = random_state()
        r = randint(40, 100)
        return f_mod(mpz_random(rstate, 2 << r).bit_set(0).bit_set(r), self.N)

    def is_not_point(self, p):
        return (isinstance(p, Ideal) == False and isinstance(p, Point) == False)

    def factorECM(self, b1, max_round):
        # build k such that k is the product of all the primes < b1
        k = 1
        for i in range(2, b1):
            if is_prime(i):
                t = 1
                while t < b1:
                    k = k*i
                    t = t*i
        # define second smoothness bound
        b2 = int(math.floor(60 * b1)) ###
        # estimated number of primes <= b2
        num_primes = math.floor(b2/math.log(b2))
        q = next_prime(b1)
        delta_max = 2

        # create delta vector with the distance between primes in [b1, b2]
        delta_vector = []
        while True:
            q1 = next_prime(q + 1)
            if q1 > b2:
                break
            dd = q1 - q
            q = q1
            if dd > delta_max:
                delta_max = dd
            delta_vector.append(dd)
        
        num_primes = len(delta_vector)
        
        print("greatest prime minor than b2 is: ", q1)

        ### matrix
        bb = []
        num_rows = int(delta_max/2)
        for row in range(num_rows):
            cols = []
            for c in range(2):
                cols.append(0)
            bb.append(cols)

        for j in range(max_round):
            self.create_curve() # generate new random curve c and point p
            e = self.curve
            p = self.point
            print("phase 1")
            p = self.mul(k, p)
            if self.is_not_point(p):
                if p is None:
                    raise BaseException('prime number')
                print("found factor here", p)
                return

            print("phase 2")
            
            
            prime_i = next_prime(b1)
            while prime_i < b2:
                q_i = self.mul(prime_i, p)
                if self.is_not_point(q_i):
                    if q_i is None:
                        raise BaseException('prime number')
                    raise BaseException('p2 is not a point, found factor: ', q_i)
                prime_i = next_prime(prime_i)
                
            print("phase two failed")
            
#             R_i_d = {}
#             R_i = []
#             for i in range(2, delta_max + 1, 2):
#                 p2 = self.mul(i, p)
#                 R_i_d[i] = p2
#                 R_i.append(p2)
#                 if self.is_not_point(p2):
#                     raise BaseException('p2 is not a point')
#             
#             #prime_i_minus_1 = next_prime(b1)
#             prime_i = next_prime(b1)
#             while prime_i < b2:
#                 prime_i_minus_1 = prime_i
#                 prime_i = next_prime(prime_i)
#                 if prime_i <= b2:
#                     R_i.append(prime_i - prime_i_minus_1)
#                     prime_i = next_prime(prime_i)
#                 else:
#                     break
#             
#                 
#             prime_i = next_prime(b1)
#             q1 = self.mul(prime_i, p)
#             if self.is_not_point(q1):
#                 raise BaseException('p2 is not a point, found factor: ', q1)
#             prime_i = next_prime(prime_i)
#             Q_i = []
#             Q_i.append(q1)
#             print(len(Q_i))
#             i = 1
#             while prime_i < b2:
#                 i += 1
#                 q_i_minus_1 = Q_i[i - 2]
#                 r_i_minus_1 = R_i[i - 2]
#                 q_i = self.partial_addition(q_i_minus_1, r_i_minus_1)
#                 if (self.is_not_point(q_i)):
#                     raise BaseException('p2 is not a point, found factor: ', q1)
#                 Q_i.append(q_i)
#                 prime_i = next_prime(prime_i)
              
              
              
              # ecm.gp              
#             q = next_prime(b1)
#             p2 = self.mul(q, p)
#             if self.is_not_point(p2):
#                 print("found factor ", p2)
#                 return
#             pp = self.partial_addition(p, p)
#             p3 = copy.deepcopy(pp)
#             for i in range(num_rows):
#                 bb[i][0] = p3.X
#                 bb[i][1] = p3.Y
#                 p3 = self.partial_addition(p3, pp)
#  
#             q1 = q
#             for i in range(num_primes):
#                 kk = int(delta_vector[i]/2)
#                 #print("kk", kk, " and bb len", len(bb))
#                 q1 = q1 + delta_vector[i]
#                 p2 = self.partial_addition(p2, Point(self.curve, bb[kk-1][0], bb[kk-1][1]))
#             print("phase 2 failed")

            

    def partial_addition(self, P, Q):
        """
        :param P: One point on EC 
        :param Q: Second Point on EC
        :return: Partial addition point P + Q according to Lenstra publication
        """
        if P.ideal() and not Q.ideal():
            return Q
        elif not P.ideal() and Q.ideal():
            return P
        elif not P.ideal() and not Q.ideal():
            N = P.curve.P
            d = gcd(P.X - Q.X, N)
            if 1 < d < N:
                return d # non trivial factor
            if d == 1: # x1 = x2
                iv = invert(f_mod(P.X - Q.X, N), N)
                delta = f_mod((P.Y - Q.Y) * iv, N)
                x_3 = f_mod(delta * delta - P.X - Q.X, N)
                y_3 = f_mod(delta * (P.X - x_3) - P.Y, N)
                return Point(P.curve, x_3, y_3)
            else:
                d = gcd(P.Y + Q.Y, N)
                if 1 < d < N:
                    return d # non trivial factor
                if d == N:
                    return Ideal() # non trivial factor
                if d == 1: # x1 != x2
                    iv = invert(f_mod(P.Y + Q.Y, N), N)
                    delta = f_mod((3 * P.X * P.X + P.curve.A) * iv, N)
                    x_3 = f_mod(delta * delta - P.X - Q.X, N)
                    y_3 = f_mod(delta * (P.X - x_3) - P.Y, N)
                    return Point(P.curve, x_3, y_3)

    def mul(self, k, p):
        """
        :param k: how many times 
        :param p: Point
        :return: k*P
        """
#         q = Point(p.getCurve(), p.getX(), p.getY())
#         try:
#             i = 1
#             while i != k:
#                 p = p.__add__(q)
#                 i += 1
#         except BaseException:
#             print("exception")
#             return i
        
        e = bin(k)[2:]  # transform to binary
        e = e[::-1]  # reverse it
        Q = p
        R = Ideal() if e[0] == '0' else p # ideal if k is even, p otherwise
        for i in e:
            Q = self.partial_addition(Q, Q)
            if self.is_not_point(Q): # neither ideal nor point
                return Q
            if i == '1':
                R = self.partial_addition(R, Q)
                if self.is_not_point(R):
                    return R
        return R

    def factor(self):
        """
        factor a number, just do it! 
        """
        k = randint(1000000, 100000000)
        Q = self.point
        x = 2
        while x < k:
            x += 1
            Q = self.mul(x, Q)
            if self.is_not_point(Q):
                return Q
        return None

if __name__ == "__main__":
    p = next_prime(10**15)
    q = next_prime(10**50)
       
    
    l = Lenstra(100000000000003700000000000000000000000000000000000000000000000000000000000000129000000000004773)
    l = Lenstra(100000000380000000361)
    # 70 target
    l = Lenstra(10000000000000000000000000000000000000000000000000000000000000000000033)
    
    l = Lenstra(1965161657812314865184863151684561324886451235479813157)
    l = Lenstra(100000000000003700000000000000000000000000000000151000000000005587)
    l = Lenstra(1000000000000037)
    l.factorECM(10000, 100)
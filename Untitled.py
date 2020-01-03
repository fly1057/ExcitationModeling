'''
方程组的矩阵形式
PE=A*E+B*U
Y=C*E
U=D*Y


'''

import numpy as np
import matplotlib.pyplot as plt


class Model():
    def __init__(self):
        super().__init__()
        self.dt = 0.01
        self.t0 = 0.0
        self.y0 = 3.0
        self.t = self.t0
        self.y = self.y0
        self.tseq = []
        self.yseq = []
        self.tseq_original = np.arange(0, 1, 0.01)
        self.yseq_original = self.tseq_original**3 - self.tseq_original**2 + 3
        self.Tq0p = 0.9
        self.Td0p = 8
        self.Tq0pp = 0.04
        self.Td0pp = 0.06
        self.a = 1
        self.b = 0.192
        self.n = 6.246
        self.ud = [0,]
        self.uq = [0.95,]
        self.Efd = [0.95,]
        self.KG0 = self.KG(self.a,self.b,self.n,self.uq[0])
        self.Edp = [0,]
        self.Eqp = [self.Efd[0]/self.KG0,]
        self.Edpp = [0,]
        self.Eqpp = [self.Efd[0]/self.KG0,]
        
        self.E = np.array([self.Edp[-1],self.Eqp[-1],self.Edpp[-1],self.Eqpp[-1]]).reshape(-1,1)#[Edp,Eqp,Edpp,Eqpp]
        self.pE = np.array([0,0,0,0]).reshape(-1,1)#-1表示我懒得计算该填什么数字，由python通过原数组和其他的值3推测出来。
        self.Y = np.array([0,0]).reshape(-1,1)#ud,uq
        self.A = np.array([[-1/self.Tq0p , 0,0,0],\
                           [0 , -self.KG(self.a,self.b,self.n,self.Eqp[-1])/self.Td0p,0,0],\
                           [(1/self.Tq0pp-1/self.Tq0p),0,-1/self.Tq0pp , 0],\
                           [0 , (1/self.Td0pp-self.KG(self.a,self.b,self.n,self.Eqp[-1])/self.Td0p),0,-1/self.Td0pp ]
                           ])
        self.B = np.array([0,1/self.Td0p,0,1/self.Td0p]).reshape(-1,1)
        self.C = np.array([1,1]).reshape(-1,1)

        print(self.A)
        print(self.B)
        print(self.C)
        print(self.E)

    def KG(self,a,b,n,Eqpp):
        return 1+b/a*Eqpp**(n-1)
        


    def rk4(self, y, f, dt, t):
        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        ynext = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        return ynext

    def f(self, t, y):
        #y=t^3-t^2+3
        #y'=3t^2-2t
        return 3 * t * t - 2 * t

    def calculate(self):
        while self.t < 1:
            self.tseq.append(self.t)
            self.yseq.append(self.y)
            self.y = self.rk4(self.y, self.f, self.dt, self.t)
            self.t = self.t + self.dt

    def plot(self):
        plt.plot(self.tseq, self.yseq, '+-')
        plt.plot(self.tseq_original, self.yseq_original, '-')
        plt.show()


if __name__ == "__main__":
    model = Model()
    model.calculate()
    model.plot()

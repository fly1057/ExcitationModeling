'''
1 方程组的矩阵形式
PE=A*E+B*U   (1-1)
Y=C*E        (1-2)
U=D*Y        (1-3)

2 步骤：
2.1 首先定义并且进行初始化
E = [ Edp , Eqp , Edpp ,Eqpp]
U = Efd
Y = [ud , uq ]
A、B是可求的，C是确定的，但是D中有两个不定参数K和uref。先给出Y根据（1-2）可求E，根据（1-1）可求U，根据U、Y可求D。
(1)给出Y
(2)求出E = inv(C)*Y
(3)求出U = invoice(B)*A*E
(4)假设D=PID(K,uref,ud,uq),设定K，则可求出uref

2.2 计算流程
(1)对E,U,Y初始化
(2)对式（1-1）调用龙格库塔计算E，调用（1-2）计算Y，调用（1-3）更新控制量U
(3)使用龙格库塔法的前提是必须将所有变量均表示成为E的函数。
因此可表示成为
PE=A*E+B*D*C*E
只有这种形式才能调用龙格库塔法，相当于只有一个变量E

'''

import numpy as np
import matplotlib.pyplot as plt


class Model():
    def __init__(self):
        super().__init__()
        self.DefaultParametersetting()

    def DefaultParametersetting(self):
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
        self.ud = [0]
        self.uq = [0.95]
        self.KG0 = self.KG(self.a, self.b, self.n, self.uq[0])
        self.Efd = [self.uq[0]*self.KG0]
        self.Edp = [0]
        self.Eqp = [self.uq[0]]
        self.Edpp = [0]
        self.Eqpp = [self.uq[0]]
        self.uref = [self.ParameterInitial()]

        self.E = np.array(
            [self.Edp[-1], self.Eqp[-1], self.Edpp[-1],
             self.Eqpp[-1]]).reshape(-1, 1)  #[Edp,Eqp,Edpp,Eqpp]
        self.pE = np.array([0, 0, 0, 0]).reshape(
            -1, 1)  #-1表示我懒得计算该填什么数字，由python通过原数组和其他的值3推测出来。
        self.Y = np.array([0, 0]).reshape(-1, 1)  #ud,uq
        self.A = np.array([[-1/self.Tq0p , 0,0,0],\
                           [0 , -self.KG(self.a,self.b,self.n,self.Eqp[-1])/self.Td0p,0,0],\
                           [(1/self.Tq0pp-1/self.Tq0p),0,-1/self.Tq0pp , 0],\
                           [0 , (1/self.Td0pp-self.KG(self.a,self.b,self.n,self.Eqp[-1])/self.Td0p),0,-1/self.Td0pp ]
                           ])
<<<<<<< HEAD
        self.B = np.array([0, 1 / self.Td0p, 0, 1 / self.Td0p])
        print(self.A)
        print(self.B)
=======
        self.B = np.array([0, 1 / self.Td0p, 0, 1 / self.Td0p]).reshape(-1, 1)
        self.C = np.array([1, 1]).reshape(-1, 1)

    def  ParameterInitial(self):
        self.Kp = 500
        self.uref = self.Efd[-1]/self.Kp+np.sqrt(self.ud[-1]**2 + self.uq[-1]**2)
        print(self.uref)
>>>>>>> f73bca7eaa453f9092d91d3c630116c93678b2c8

    def PID(self, Kp , uref, ud, uq):
            return Kp* (uref - np.sqrt(ud**2 + uq**2))

        # print(self.A)
        # print(self.B)
        # print(self.C)
        # print(self.E)

    def KG(self, a, b, n, Eqpp):
        return 1 + b / a * Eqpp**(n - 1)

    def rk4(self, y, f, dt, t):
        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        ynext = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        return ynext

    def test_f(self, t, y):
        #y=t^3-t^2+3
        #y'=3t^2-2t
        return 3 * t * t - 2 * t

    def f(self,t ,E):
        return self.A*self.E+self.B

    def test_calculate(self):
        while self.t < 1:
            self.tseq.append(self.t)
            self.yseq.append(self.y)
            self.y = self.rk4(self.y, self.test_f, self.dt, self.t)
            self.t = self.t + self.dt

    def test_plot(self):
        plt.plot(self.tseq, self.yseq, '+-')
        plt.plot(self.tseq_original, self.yseq_original, '-')
        plt.show()


if __name__ == "__main__":
    model = Model()
    model.test_calculate()
    model.test_plot()

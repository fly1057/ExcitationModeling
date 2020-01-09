'''
1 方程组的矩阵形式
PE=A*E+B*U   (1-1)

2 步骤：
2.1 首先定义并且进行初始化
E = [ Edp , Eqp , Edpp ,Eqpp , Efd1 ,Efd2 ,Efd3 ,Efd4]
 Efd1 ,Efd2 ,Efd3 ,Efd4为由于励磁PID导致的状态变量


2.2 计算流程
(1)对E,U初始化
(2)对式（1-1）调用龙格库塔计算E，所有的变量都是E或者E中的量，这样才能进行龙格库塔计算，因为只有E一个变量

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['font.sans-serif'] = ['SimHei']  # 步骤一（替换sans-serif字体）
plt.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）


class Model():
    def __init__(self):
        super().__init__()
        self.DefaultParametersetting()

    def DefaultParametersetting(self):
        self.dt = 0.01
        self.t = 0.0
        self.y = 3.0
        self.tseq = [0.0]
        self.yseq = [3.0]
        self.tseq_original = np.arange(0, 1, 0.01)
        self.yseq_original = self.tseq_original**3 - self.tseq_original**2 + 3

        #控制参数
        self.Tr = 0.01
        self.TA = 0.01
        self.K = 22*7.84
        self.T1 = 1.0
        self.T2 = 4.0
        self.T3 = 1.0
        self.T4 = 1.0
        #发电机参数
        self.Tq0p = 0.2
        self.Td0p = 9.45
        self.Tq0pp = 0.198   #0.198
        self.Td0pp = 0.091
        self.a = 1.0
        self.b = 0.192
        self.n = 6.246
        #需求变量定义及初始化
        self.ud0 = 0
        self.uq0 = 0.95
        self.Edp0 = 0
        self.Eqp0 = 0.95
        self.Edpp0 = 0
        self.Eqpp0 = 0.95
        self.KG0 = 1 + self.b / self.a * self.Eqpp0**(self.n - 1)
        self.Efd0 = self.uq0 * self.KG0
        self.uref0 = self.Efd0/self.K+self.uq0
        self.deltaU = 0.0488
        self.Tstart = 1
        self.Tend = 6
        self.Efd10 = self.uref0-self.uq0
        self.Efd20 = self.Efd0
        self.Efd30 = self.Efd0
        self.Efd40 = self.Efd0
        
        self.dtvector = np.array([0.001]*8).reshape(-1,1)
        self.tvector = np.array([0]*8).reshape(-1,1)
        self.tmatrix = np.array([0]*8).reshape(-1,1)

        self.Evector = np.array([self.Edp0, self.Eqp0, self.Edpp0,self.Eqpp0,self.Efd10,self.Efd20,self.Efd30,self.Efd40]).reshape(-1, 1)  #[Edp,Eqp,Edpp,Eqpp]
        self.Ematrix = self.Evector

        self.pEvector = np.array([0]*8).reshape(-1, 1)  #-1表示我懒得计算该填什么数字，由python通过原数组和其他的值3推测出来。
        self.pEmatrix = self.pEvector

        self.A = np.array([[-1/self.Tq0p , 0,0,0,0,0,0,0],\
                           [0 , -(1 + self.b / self.a * float(self.Evector[3])**(self.n - 1))/self.Td0p,0,0,0,0,0,1/self.Td0p],\
                           [(1/self.Tq0pp-1/self.Tq0p),0,-1/self.Tq0pp ,0, 0,0,0,0],\
                           [0 , (1/self.Td0pp-(1 + self.b / self.a * float(self.Evector[3])**(self.n - 1))/self.Td0p),0,-1/self.Td0pp,0,0,0,1/self.Td0p],\
                           [0,0,0,0,-1/self.Tr,0 ,0,0],\
                           [0,0,0,0,self.K/self.TA,-1/self.TA ,0,0],\
                           [0,0,0,0,self.T1/self.T2*self.K/self.TA,1/self.T2-self.T1/(self.T2*self.TA),-1/self.T2,0],\
                           [0,0,0,0,self.T1/self.T2*self.T3/self.T4*self.K/self.TA ,self.T3/self.T4*(1/self.T2-self.T1/(self.T2*self.TA)),(1/self.T4-self.T3/(self.T4*self.T2)),-1/self.T4]
                           ])
        self.B = np.array([0,0,0,0,1,0,0,0]).reshape(-1, 1)

        print("A = \n",self.A)
        print("B = \n",self.B)
        print("E = \n",self.Ematrix)

    def rk4(self, y, f, dt, t):

        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        ynext = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        return ynext

    def uref_fun(self,t):

        if t< self.Tstart :
            temp = self.uref0
        elif t<self.Tstart+0.02:
            temp = self.uref0 + self.deltaU*(t-self.Tstart)/0.02
        else :#这边瞬间上升有点问题
            temp = self.uref0 + self.deltaU
        return temp

    def duref_fun(self,t):
        
        if t< self.Tstart :
            temp = 0
        elif t<self.Tstart+0.02:
            temp = self.deltaU/0.02
        else :#这边瞬间上升有点问题
            temp = 0
        return temp

    def test_f(self, t, y):
        # y=t^3-t^2+3
        # y'=3t^2-2t
        return 3 * t * t - 2 * t

    def test_f2(self, t, Evector):

        return self.A@Evector+self.B*((self.uref_fun(self.tvector[0])-float(np.sqrt(Evector[2]**2+Evector[3]**2)))/self.Tr+self.duref_fun(self.tvector[0]) )

    def test_calculate2(self):
        try:
            while self.tvector[0] < self.Tend:
                self.Evector = self.rk4(self.Evector, self.test_f2, self.dtvector, self.tvector)
                self.Ematrix = np.hstack((self.Ematrix,self.Evector))
                self.tvector = self.tvector + self.dtvector
                self.tmatrix = np.hstack((self.tmatrix , self.tvector))
        except Exception as e:
            print(e)
            print(e.__traceback__.tb_frame.f_globals["__file__"])  # 发生异常所在的文件
            print(e.__traceback__.tb_lineno)  # 发生异常所在的行数

    def test_calculate(self):
        try:
            while self.t < 1:
                self.yseq.append( self.rk4(self.yseq[-1], self.test_f, self.dt, self.tseq[-1]) )
                self.t = self.t + self.dt
                self.tseq.append(self.t)
        except Exception as e:
            print(e)
            print(e.__traceback__.tb_frame.f_globals["__file__"])  # 发生异常所在的文件
            print(e.__traceback__.tb_lineno)  # 发生异常所在的行数

    def test_plot(self):
        plt.plot(self.tseq, self.yseq, '+-')
        plt.plot(self.tseq_original, self.yseq_original, '-')
        plt.show()

    def test_plot2(self):
        plt.plot(self.tmatrix[1,:],self.Ematrix[1,:], '-')
        plt.show()


if __name__ == "__main__":
    #读取实测采样波形
    df = pd.read_csv('C:/Users/ll/Desktop/zaoshistep.csv')
    #构造仿真数据储存格式
    df2=pd.DataFrame
    meas_t = df["t"]
    meas_ug = df["UAB2"]

    model = Model()
    model.test_calculate2()
    # model.test_plot2()
    #计算后仿真结果保存csv
    temp = np.vstack((model.tmatrix[1,:],model.Ematrix))#列向量横向拼凑
    #重新构造dataframe的列名
    df2=pd.DataFrame(temp.transpose(),columns=['t','Edp','Eqp','Edpp','Eqpp','Efd1','Efd2','Efd3','Efd4'])
    #dataframe保存时将索引去掉
    df2.to_csv('C:/Users/ll/Desktop/zaoshistepsimulate考虑Tr.csv',index=0)

    #绘图
    plt.plot(meas_t,meas_ug)
    plt.plot(model.tmatrix[1,:],model.Ematrix[1,:], '-')
    plt.legend(["实测录波","python仿真"])
    plt.grid()
    plt.show()
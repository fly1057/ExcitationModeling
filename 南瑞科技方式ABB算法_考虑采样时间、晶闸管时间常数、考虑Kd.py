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
(3)南瑞科技的控制模型里面仍然使用cos=K(uref-ug)/ug  里面的量并未涉及励磁系统的定标，串联传递函数也未涉及
励磁系统的定标，那么也就是说厂家用什么定标实际上对仿真没有影响。仿真实际上做的是仿真的框架。
'''

import os 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#plt.rcParams['font.sans-serif'] = ['SimHei']  # 步骤一（替换sans-serif字体）
#plt.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）


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
########################################################################################
        #控制参数  五强溪3数据
        self.Tr = 0.02
        self.TA = 0.01
        self.K = 1
        self.Kv = 0
        self.Kp = 60
        self.Ki = 20
        self.Kd = 1
        self.Td = 0.01
        #将并联PID转化为串联PID
        self.ax = (self.Kp*self.Td + self.Kd )/self.Ki
        self.bx = self.Td + self.Kp/self.Ki
        self.cx = 1
        self.x1 = (-self.bx + np.sqrt((self.bx)**2-4*self.ax*self.cx))/(2*self.ax)
        self.x2 = (-self.bx - np.sqrt((self.bx)**2-4*self.ax*self.cx))/(2*self.ax)
        self.T1 = -1/self.x1
        self.T2 = 1/self.Ki
        self.T3 = -1/self.x2
        self.T4 = self.Td

        #发电机参数 五强溪3数据
        self.Tq0p = 0.001
        self.Td0p = 8.37
        self.Tq0pp = 0.06   
        self.Td0pp = 0.05
        self.a = 1.0
        self.b = 0.081
        self.n = 9.966

        # #控制参数 湘潭3数据
        # self.Tr = 0.02
        # self.TA = 0.01
        # self.K = 500
        # self.T1 = 1.0
        # self.T2 = 8.33
        # self.T3 = 1.0
        # self.T4 = 1.0
        # #发电机参数 湘潭3数据
        # self.Tq0p = 0.95
        # self.Td0p = 9.32
        # self.Tq0pp = 0.069   
        # self.Td0pp = 0.045
        # self.a = 1.0
        # self.b = 0.186
        # self.n = 8.357
#########################################################################################
        #需求变量定义及初始化
        self.ud0 = 0
        self.uq0 = 0.94765
        self.Edp0 = 0
        self.Eqp0 = self.uq0
        self.Edpp0 = 0
        self.Eqpp0 = self.uq0
        self.KG0 = 1 + self.b / self.a * self.Eqpp0**(self.n - 1)
        self.Efd0 = self.uq0 * self.KG0
        self.uref0 = self.uq0  #关键
        self.deltaU = 0.05
        self.Tstart = 1
        self.Tend = 8
        self.Tstepdelay = 0.02 #不可为0，否则将除以0 
        self.Efd10 = 0 #关键
        self.Efd20 = 0 #关键
        self.Efd30 = self.Efd0 #关键
        self.Efd40 = self.Efd0 #关键
        
        self.dtvector = np.array([0.005]*8).reshape(-1,1)
        self.tvector = np.array([0]*8).reshape(-1,1)
        self.tmatrix = np.array([0]*8).reshape(-1,1)

        self.Evector = np.array([self.Edp0, self.Eqp0, self.Edpp0,self.Eqpp0,self.Efd10,self.Efd20,self.Efd30,self.Efd40]).reshape(-1, 1)  #[Edp,Eqp,Edpp,Eqpp]
        self.Ematrix = self.Evector

        self.pEvector = np.array([0]*8).reshape(-1, 1)  #-1表示我懒得计算该填什么数字，由python通过原数组和其他的值3推测出来。
        self.pEmatrix = self.pEvector

        #只需要改变A阵即可
        self.A = np.array([[-1/self.Tq0p , 0,0,0,0,0,0,0],\
                           [0 , -(1 + self.b / self.a * float(self.Evector[3])**(self.n - 1))/self.Td0p,0,0,0,0,0,1/self.Td0p],\
                           [(1/self.Tq0pp-1/self.Tq0p),0,-1/self.Tq0pp ,0, 0,0,0,0],\
                           [0 , (1/self.Td0pp-(1 + self.b / self.a * float(self.Evector[3])**(self.n - 1))/self.Td0p),0,-1/self.Td0pp,0,0,0,1/self.Td0p],\
                           [0,0,0,0,-1/self.Tr,0 ,0,0],\
                           [0,0,0,0,self.K/self.TA,-1/self.TA ,0,0],\
                           [0,0,0,0,self.T1/self.T2*self.K/self.TA,1/self.T2-self.T1/(self.T2*self.TA),-self.Kv/self.T2,0],\
                           [0,0,0,0,self.T1/self.T2*self.T3/self.T4*self.K/self.TA ,self.T3/self.T4*(1/self.T2-self.T1/(self.T2*self.TA)),(1/self.T4-self.Kv*self.T3/(self.T4*self.T2)),-1/self.T4]
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
        elif t<self.Tstart+self.Tstepdelay:
            temp = self.uref0 + self.deltaU*(t-self.Tstart)/self.Tstepdelay
        else :
            temp = self.uref0 + self.deltaU
        return temp

    def duref_fun(self,t):
        
        if t< self.Tstart :
            temp = 0
        elif t<self.Tstart+self.Tstepdelay:
            temp = self.deltaU/self.Tstepdelay
        else :
            temp = 0
        return temp

    def test_f(self, t, y):
        # y=t^3-t^2+3
        # y'=3t^2-2t
        return 3 * t * t - 2 * t

    def test_f2(self, t, Evector):#即为rk4中需要的导数函数
        Edp = self.Evector[0]
        Eqp = self.Evector[1]
        Edpp = self.Evector[2]
        Eqpp = self.Evector[3]
        Efd = self.Evector[7]
        ug = float(np.sqrt(Edpp**2+Eqpp**2))
        uref = self.uref_fun(self.tvector[0])
        duref = self.duref_fun(self.tvector[0])
        
        temp = (uref*ug-ug**2)/self.Tr
        + ug*duref + \
            uref/ug*( \
            (1/self.Tq0pp - 1/self.Tq0p)*Edp - \
            1/self.Tq0p*Edpp + \
            (1/self.Td0pp -  (1 + self.b / self.a * float(Eqpp)**(self.n - 1)))/self.Td0p*Eqp- \
            1/self.Td0pp*Eqpp+1/self.Td0p*Efd)
        return self.A@Evector+self.B*(temp)

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
    df = pd.read_csv((os.getcwd()).replace("\\","/")+'/wuqiangxi3step2.csv')
    #df = pd.read_csv('C:/Users/ll/Desktop/xiangtan3.csv')
    #构造仿真数据储存格式
    df2=pd.DataFrame
    meas_t = df["t"]
    meas_ug = df["UAB"]

    model = Model()
    model.test_calculate2()
    # model.test_plot2()
    #计算后仿真结果保存csv
    temp = np.vstack((model.tmatrix[1,:],model.Ematrix))#列向量横向拼凑
    #重新构造dataframe的列名
    df2=pd.DataFrame(temp.transpose(),columns=['t','Edp','Eqp','Edpp','Eqpp','Efd1','Efd2','Efd3','Efd4'])
    #dataframe保存时将索引去掉
    df2.to_csv((os.getcwd()).replace("\\","/")+'/wuqiangxi3stepsimulate考虑TrTA 南瑞科技方式.csv',index=0)

    #绘图
    plt.plot(meas_t,meas_ug, linewidth = '1', label = "test1", linestyle='-', marker='')
    plt.plot(model.tmatrix[1,:],model.Ematrix[3,:], linewidth = '2', label = "test2",  linestyle='--')
    plt.legend(["Measurement","python  Simulation"])
    plt.title("NES5100 device")
    plt.grid()
    plt.show()
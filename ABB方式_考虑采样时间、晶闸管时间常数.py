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
import  os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['font.sans-serif'] = ['SimHei']  # 步骤一（替换sans-serif字体）
plt.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）


class Model():
    def __init__(self):
        super().__init__()
        self.DataInput()

    def DataInput(self):

        # #控制参数  皂市1数据
        # self.Tr = 0.01
        # self.TA = 0.01
        # self.K = 22*7.84
        # self.T1 = 1.0
        # self.T2 = 4.0
        # self.T3 = 1.0
        # self.T4 = 1.0
        # #发电机参数 皂市1数据
        # self.Tq0p = 0.2
        # self.Td0p = 9.45
        # self.Tq0pp = 0.198   #0.198
        # self.Td0pp = 0.091
        # self.a = 1.0
        # self.b = 0.192
        # self.n = 6.246
        #励磁控制参数 湘潭3数据
        self.Tr = 0.02
        self.TA = 0.01
        self.K = 500
        self.Kv = 1
        self.ModeFlag = "EXC9000"
        #如果Kv等于0，那么状态变量需要重新初始化，另外T1之类需要重新计算，其余和串联PID一致
        self.T1 = 1.0
        self.T2 = 8.33
        self.T3 = 1.0
        self.T4 = 1.0
        #发电机参数 湘潭3数据
        self.Tq0p = 0.95
        self.Td0p = 9.32
        self.Tq0pp = 0.069   
        self.Td0pp = 0.045
        self.a = 1.0
        self.b = 0.186
        self.n = 8.357
        # #控制参数  鸟儿巢1数据
        # self.Tr = 0.02
        # self.TA = 0.01
        # self.K = 18*10.6 #4.69*36 #
        # self.T1 = 1.0
        # self.T2 = 20.0
        # self.T3 = 0.2
        # self.T4 = 0.05
        # #发电机参数 鸟儿巢1数据
        # self.Tq0p = 0.0001
        # self.Td0p = 4.25
        # self.Tq0pp = 0.035   
        # self.Td0pp = 0.042
        # self.a = 1.0
        # self.b = 0.118
        # self.n = 7.1652
#########################################################################################
        #######需求变量定义及初始化###########
        #发电机变量的初始化
        self.ud0 = 0
        self.uq0 = 0.9493  
        #从0.95开始封脉冲计算Td0p不对，必须从线性段封脉冲才行比如0.4左右，从0.7开始误差也比较大
        self.Edp0 = 0
        self.Eqp0 = self.uq0
        self.Edpp0 = 0
        self.Eqpp0 = self.uq0
        self.KG0 = 1 + self.b / self.a * self.Eqpp0**(self.n - 1)
        self.Efd0 = self.uq0 * self.KG0

        #励磁状态变量初始化
        self.Efd40 = self.Efd0
        self.Efd30 = self.Efd40
        self.Efd20 = self.Kv*self.Efd30
        self.Efd10 = self.Efd20/self.K

        if self.ModeFlag =="ABB":
            self.uref0 = self.Kv*self.Efd0/self.K + self.uq0 #EXC9000 初始化关键

        elif self.ModeFlag =="EXC9000":
            self.uref0 = self.Kv*self.Efd0/(self.K*self.uq0) + self.uq0 #ABB 初始化关键

        elif self.ModeFlag =="NES5000":
            self.uref0 = self.Kv*self.Efd0/(self.K*self.uq0) + self.uq0 #ABB 初始化关键

        #阶跃变量初始化
        self.deltaU = 0.0498  #不影响封脉冲
        self.Tstart = 1
        self.Tend = 20
        self.Tstepdelay = 0.02 #不可为0，否则将除以0

        #快速形成列向量的方法
        self.dtvector = np.array([0.01]*8).reshape(-1,1)
        self.tvector = np.array([0]*8).reshape(-1,1)
        self.tmatrix = np.array([0]*8).reshape(-1,1)

        self.Evector = np.array([self.Edp0, self.Eqp0, self.Edpp0,self.Eqpp0,self.Efd10,self.Efd20,self.Efd30,self.Efd40]).reshape(-1, 1)  #[Edp,Eqp,Edpp,Eqpp]
        self.Ematrix = self.Evector

        self.pEvector = np.array([0]*8).reshape(-1, 1)  #-1表示懒得计算该填什么数字，由python通过原数组和其他的值推测出来。
        self.pEmatrix = self.pEvector

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

    #uref的定义和duref的定义也很有意思
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

        if self.ModeFlag =="ABB":
            tempABB = (uref*ug-ug**2)/self.Tr + \
                ug*duref + \
                uref/ug*( \
                (1/self.Tq0pp - 1/self.Tq0p)*Edp - \
                1/self.Tq0pp*Edpp + \
                (1/self.Td0pp - (1 + self.b / self.a * float(Eqpp)**(self.n - 1))/self.Td0p) *Eqp - \
                1/self.Td0pp*Eqpp + \
                1/self.Td0p*Efd)
            temp = tempABB

        elif self.ModeFlag =="EXC9000":
            tempEXC9000 = (uref-ug)/self.Tr+duref
            temp = tempEXC9000

        elif self.ModeFlag =="NES5000":
            tempNES5000 = (uref-ug)/self.Tr+duref
            temp = tempNES5000

        #用于模拟封脉冲，让所有的控制不起作用，仅让电机本体起作用，在这里B阵乘的应该为0
        return self.A@Evector+self.B*(temp)

    def test_calculate2(self):#实现利用rk4进行计算，然后更新各向量或者矩阵
        try:
            while self.tvector[0] < self.Tend:
            #rk4的格式为rk4( y, f, dt, t)
                ############################用于模拟封脉冲
                # self.Evector[4] = 0
                # self.Evector[5] = 0
                # self.Evector[6] = 0
                # self.Evector[7] = 0
                ############################用于模拟封脉冲
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
    #这里采用了系统还是来得到当前文件的路径，从中将默认的分隔符进行了替换
    df = pd.read_csv((os.getcwd()).replace("\\","/")+'/xiangtan3step.csv')
    #df = pd.read_csv((os.getcwd()).replace("\\","/")+'/niaoerchao1step.csv')
    #构造仿真数据储存格式
    df2=pd.DataFrame
    meas_t = df["t"]
    meas_ug = df["UAB"]

    model = Model()
    model.test_calculate2()
    # model.test_plot2()
    #计算后仿真结果保存csv
    #这里利用了vstack这个函数，很方便的进行列矩阵的构造
    temp = np.vstack((model.tmatrix[1,:],model.Ematrix))#列向量横向拼凑
    #重新构造dataframe的列名
    df2=pd.DataFrame(temp.transpose(),columns=['t','Edp','Eqp','Edpp','Eqpp','Efd1','Efd2','Efd3','Efd4'])
    #dataframe保存时将索引去掉
    df2.to_csv((os.getcwd()).replace("\\","/")+'/xiangtan3stepsimulate考虑TrTA ABB方式.csv',index=0)
    #df2.to_csv((os.getcwd()).replace("\\","/")+'/niaoerchao1stepsimulate考虑TrTA ABB方式.csv',index=0)

    #绘图
    plt.plot(meas_t,meas_ug, linewidth = '1', label = "test1", linestyle='-', marker='')
    plt.plot(model.tmatrix[1,:],model.Ematrix[1,:], linewidth = '2', label = "test2",  linestyle='--')
    #plt.plot(model.tmatrix[1,:],model.Ematrix[7,:], linewidth = '2', label = "test2",  linestyle='--')
    plt.legend(["Measurement Ug","Simulation Ug","Simulation Efd"])
    plt.title("湘潭3 ABB5000 ")
    plt.xlabel("t/s")
    plt.ylabel("Ut/p.u.")
    plt.grid()
    plt.show()
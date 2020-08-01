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

黔东的计算不对，有两个原因一种是采用了MF卡的方程，另外一种是采用了MG卡的方程

'''
import  os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#plt.rcParams['font.sans-serif'] = ['SimHei']  # 步骤一（替换sans-serif字体）
#plt.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）


class Model():
    def __init__(self):
        super().__init__()
        
    def CalculateInitial(self,InitialUg,Stepdt):
        #######需求变量定义及初始化###########
        #发电机变量的初始化
        self.ud0 = 0
        self.uq0 = InitialUg #0.9493  
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

        if self.ExciterType =="ABB":
            self.uref0 = self.Kv*self.Efd0/(self.K*self.uq0) + self.uq0

        elif (self.ExciterType =="EXC9000") or (self.ExciterType =="NES5000"):
            self.uref0 = self.Kv*self.Efd0/self.K + self.uq0 #ABB 初始化关键

        #快速形成列向量的方法
        self.dtvector = np.array([Stepdt]*8).reshape(-1,1)  # Stepdt=0.01
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
    
    def ExciterInitial(self,Tr,TA,K,Kv,ExciterType,T1,T2,T3,T4):
        # 励磁控制参数 湘潭3数据 
        # self.Tr = 0.02
        # self.TA = 0.01
        # self.K = 500
        # self.Kv = 1
        # self.ExciterType = "EXC9000"
        # #如果Kv等于0，那么状态变量需要重新初始化，另外T1之类需要重新计算，其余和串联PID一致
        # self.T1 = 1.0
        # self.T2 = 8.33
        # self.T3 = 1.0
        # self.T4 = 1.0
        #励磁控制参数 湘潭3数据 
        self.ExciterType = ExciterType
        self.Tr = Tr
        self.TA = TA
        self.K = K
        self.Kv = Kv
        if (self.ExciterType =="ABB") or (self.ExciterType =="EXC9000"):
        #如果Kv等于0，那么状态变量需要重新初始化，另外T1之类需要重新计算，其余和串联PID一致
            self.T1 = T1
            self.T2 = T2
            self.T3 = T3
            self.T4 = T4


        elif self.ExciterType =="NES5000":
            self.Kp = T1  #仍然采用串联PID的赋值方式
            self.Ki = T2
            self.Kd = T3
            self.Td = T4

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
            


    def GeneratorInitial(self,Tq0p,Td0p,Tq0pp,Td0pp,a,b,n):
        #发电机参数 湘潭3数据
        # self.Tq0p = 0.95
        # self.Td0p = 9.32
        # self.Tq0pp = 0.069   
        # self.Td0pp = 0.045
        # self.a = 1.0
        # self.b = 0.186
        # self.n = 8.357
        #发电机参数 湘潭3数据
        self.Tq0p = Tq0p
        self.Td0p = Td0p
        self.Tq0pp = Tq0pp   
        self.Td0pp = Td0pp
        self.a = a
        self.b = b
        self.n = n

    def DisturbanceInitial(self,DeltaU,Tstart,Tend,Tstepdelay):
        # #阶跃变量初始化
        # self.DeltaU = 0.0497  #不影响封脉冲
        # self.Tstart = 1
        # self.Tend = 20
        # self.Tstepdelay = 0.02 #不可为0，否则将除以0
        #阶跃变量初始化
        self.DeltaU = DeltaU  #不影响封脉冲
        self.Tstart = Tstart
        self.Tend = Tend
        self.Tstepdelay = Tstepdelay #不可为0，否则将除以0

    def rk4(self, y, f, dt, t):

        k1 = f(t, y)
        k2 = f(t + dt / 2, y + dt / 2 * k1)
        k3 = f(t + dt / 2, y + dt / 2 * k2)
        k4 = f(t + dt, y + dt * k3)
        ynext = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        return ynext

    #uref的定义和duref的定义也很有意思
    def uref_fun(self,t):

        if t<= self.Tstart :
            temp = self.uref0
        elif (t> self.Tstart) and (t<=self.Tstart+self.Tstepdelay) :
            temp = self.uref0 + self.DeltaU*(t-self.Tstart)/self.Tstepdelay
        else :
            temp = self.uref0 + self.DeltaU

        return temp

    def duref_fun(self,t):
        
        if t<= self.Tstart :
            temp = 0
        elif (t> self.Tstart) and (t<=self.Tstart+self.Tstepdelay) :
            temp = self.DeltaU/self.Tstepdelay
        else :
            temp = 0

        return temp

    def DifferentialEquation(self, t, Evector):#即为rk4中需要的导数函数
        Edp = self.Evector[0]
        Eqp = self.Evector[1]
        Edpp = self.Evector[2]
        Eqpp = self.Evector[3]
        Efd = self.Evector[7]
        ug = float(np.sqrt(Edpp**2+Eqpp**2))
        uref = self.uref_fun(self.tvector[0])
        duref = self.duref_fun(self.tvector[0])

        if self.ExciterType =="ABB":
            temp = (uref*ug-ug**2)/self.Tr + \
                ug*duref + \
                uref/ug*( \
                (1/self.Tq0pp - 1/self.Tq0p)*Edp - \
                1/self.Tq0pp*Edpp + \
                (1/self.Td0pp - (1 + self.b / self.a * float(Eqpp)**(self.n - 1))/self.Td0p) *Eqp - \
                1/self.Td0pp*Eqpp + \
                1/self.Td0p*Efd)

        elif (self.ExciterType =="EXC9000") or (self.ExciterType =="NES5000"):
            temp = (uref-ug)/self.Tr+duref

        #用于模拟封脉冲，让所有的控制不起作用，仅让电机本体起作用，在这里B阵乘的应该为0
        return self.A@Evector+self.B*(temp)

    def Calculate(self):#实现利用rk4进行计算，然后更新各向量或者矩阵
        try:
            while self.tvector[0] < self.Tend:
            #rk4的格式为rk4( y, f, dt, t)
                ############################用于模拟封脉冲
                # self.Evector[4] = 0
                # self.Evector[5] = 0
                # self.Evector[6] = 0
                # self.Evector[7] = 0
                ############################用于模拟封脉冲
                self.Evector = self.rk4(self.Evector, self.DifferentialEquation, self.dtvector, self.tvector)
                self.Ematrix = np.hstack((self.Ematrix,self.Evector))
                self.tvector = self.tvector + self.dtvector
                self.tmatrix = np.hstack((self.tmatrix , self.tvector))
        except Exception as e:
            print(e)
            print(e.__traceback__.tb_frame.f_globals["__file__"])  # 发生异常所在的文件
            print(e.__traceback__.tb_lineno)  # 发生异常所在的行数

if __name__ == "__main__":
    #读取实测采样波形
    #这里采用了系统还是来得到当前文件的路径，从中将默认的分隔符进行了替换
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/xiangtan3step.csv')
    #df = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/niaoerchao1step.csv')
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/wanmipo3step.csv')
    # df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/wuqiangxi3step.csv')
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/qiandong1step.csv')
    # df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/zaoshi1step.csv')
    # df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/baishi2step.csv')
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/dayuandu2step.csv')
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/sanbanxi1step.csv')
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/yueyanglaji1step.csv')  
    #df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/leiyang3step.csv')
    df_measurement = pd.read_csv((os.getcwd()).replace("\\","/")+'/data/yiyang1step.csv')
    #构造仿真数据储存格式
    meas_t = df_measurement["t"]
    meas_ug = df_measurement["UAB"]

    model = Model()

    #################声明################
    #def ExciterInitial(self,Tr,TA,K,Kv,ExciterType,T1,T2,T3,T4)
    #def GeneratorInitial(self,Tq0p,Td0p,Tq0pp,Td0pp,a,b,n)
    #def DisturbanceInitial(self,DeltaU,Tstart,Tend,Tstepdelay)
    #def CalculateInitial(self,InitialUg,Stepdt)
    #####################################
    ############湘潭3号机#####################
    #############碗米坡3号机####################
    # model.ExciterInitial(0.01,0.01,500,1,"ABB",0.01,0.0025,1.5,15)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.001,5.36,0.068,0.0466,1,0.141,7.386)
    # model.DisturbanceInitial(0.05,1,6,0.02)#0.04985
    # model.CalculateInitial(0.95,0.001)
    #############碗米坡3号机####################
    #############五强溪3号机####################
    # model.ExciterInitial(0.01,0.01,1.2,0,"NES5000",60,20,1,0.01)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.0001,8.37,0.06,0.05,1,0.081,9.966)
    # model.DisturbanceInitial(0.0498,1,5,0.01)#0.04985
    # model.CalculateInitial(0.95,0.001)  
    #############五强溪3号机####################
    # ###############皂市1号机####################
    # model.ExciterInitial(0.01,0.01,22*7.84,1,"ABB",1,4,1,1)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.001,9.45,0.198,0.091,1,0.192,6.246)
    # model.DisturbanceInitial(0.0481,1.018,5,0.01)#0.04985
    # model.CalculateInitial(0.95,0.005)
    # ##############皂市1号机####################
    # #############凌津滩5号机####################
    # model.ExciterInitial(0.01,0.01,500,1,"ABB",0.01,0.0025,1.5,15)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.001,3.1,0.03,0.02,1,0.245,6.225)
    # model.DisturbanceInitial(0.05,1.02,5,0.01)#0.04985
    # model.CalculateInitial(0.95,0.005)
    ##############凌津滩5号机####################
    # #############黔东1号机####################
    # model.ExciterInitial(0.01,0.01,1,0,"NES5000",60,10,0,0.01)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.94,7.47,0.069,0.045,1,0.0605,13.646)
    # model.DisturbanceInitial(0.0498,1,5,0.01)#0.04985
    # model.CalculateInitial(0.95,0.001)
    # #############黔东1号机####################
    #############白市2号机####################
    # model.ExciterInitial(0.01,0.01,1,0,"NES5000",90,30,0,0.01)#TA和Tr用0.02误差很大，用0.01,0.01行
    # model.GeneratorInitial(0.1,7.24,0.121,0.055,1,0.129,6.89)
    # model.DisturbanceInitial(0.0498,1,5,0.01)#0.04985
    # model.CalculateInitial(0.9518,0.001)
    #############白市2号机####################
    #############大源渡2号机####################
    #model.ExciterInitial(0.01,0.01,22*4,1,"EXC9000",1,4,1.0,1.0)
    #model.GeneratorInitial(0.0001,3.5,0.024,0.06,1,0.113,8.445)
    #model.DisturbanceInitial(0.05,1,6,0.02)
    #model.CalculateInitial(0.95,0.0005)
    #############################################
    ##############三板溪1号机####################
    #model.ExciterInitial(0.01,0.01,300/1.12,1,"ABB",1,6.667 ,1.0,1.0)
    #model.GeneratorInitial(0.001,11.2,0.04,0.06,1,0.075,9.472)
    #model.DisturbanceInitial(0.050,1,6,0.02)
    #model.CalculateInitial(0.95,0.005)
    #############################################
    #############岳阳垃圾1号机####################
    #model.ExciterInitial(0.01,0.01,50*8.62,1,"ABB",1,8.33 ,1.0,1.0)
    #model.GeneratorInitial(0.95,10.1,0.04,0.06,1,0.12,8.58)
    #model.DisturbanceInitial(0.0503,1,6,0.02)
    #model.CalculateInitial(0.95,0.005)
    ############################################# 
    ###############耒阳3号机####################  
    #model.ExciterInitial(0.01,0.01,500,1,"ABB",1.8, 11.25 ,1.0,1.0)
    #model.GeneratorInitial(0.94,8.19,0.079,0.045,1,0.135,10.4)
    #model.DisturbanceInitial(0.0493,1,6,0.01)
    #model.CalculateInitial(0.95,0.002)
    ##############益阳1号机####################  
    model.ExciterInitial(0.01,0.01,40*1.875*9.3,1,"EXC9000",3, 30 ,1.0,1.0)
    model.GeneratorInitial(0.94,9.2,0.079,0.035,1,0.129,8.124)
    model.DisturbanceInitial(0.05,1.02,10,0.01)
    model.CalculateInitial(0.95,0.002)
    #################声明################
    #def ExciterInitial(self,Tr,TA,K,Kv,ExciterType,T1,T2,T3,T4)
    #def GeneratorInitial(self,Tq0p,Td0p,Tq0pp,Td0pp,a,b,n)
    #def DisturbanceInitial(self,DeltaU,Tstart,Tend,Tstepdelay)
    #def CalculateInitial(self,InitialUg,Stepdt)
    #####################################
    
    model.Calculate()
    #计算后仿真结果保存csv，这里利用了vstack这个函数，很方便的进行列矩阵的构造
    temp = np.vstack((model.tmatrix[1,:],model.Ematrix))#列向量横向拼凑
    #重新构造dataframe的列名
    df_simulation=pd.DataFrame(temp.transpose(),columns=['t','Edp','Eqp','Edpp','Eqpp','Efd1','Efd2','Efd3','Efd4'])
    #dataframe保存时将索引去掉
    df_simulation.to_csv((os.getcwd()).replace("\\","/")+'/stepsimulatetemp.csv',index=0)

    #绘图
    plt.plot(meas_t,meas_ug, linewidth = '1', label = "test1", linestyle='-', marker='')
    plt.plot(model.tmatrix[1,:],model.Ematrix[1,:], linewidth = '2', label = "test2",  linestyle='--')
    #plt.plot(model.tmatrix[1,:],model.Ematrix[7,:], linewidth = '2', label = "test2",  linestyle='--')
    plt.legend(["Measurement Ug","Simulation Ug","Simulation Efd"])
    plt.title("Measurement and Simulation ")
    plt.xlabel("t/s")
    plt.ylabel("Ut/p.u.")
    plt.grid()
    plt.show()

    
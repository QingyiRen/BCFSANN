%VERSION of May 13, 2010,by phD Zhijun Zhang and MS Jun Lee
clc;
%clear all;
format long;
global alpha r q0 T lambda sigma aa alpha_i a_i d_i D sa ca s2a c2a pen_long K k;

%circle data
T=5.0;             %ͼ����ʾʱ��
r=0.45;
alpha=pi/6;
sigma=pi/1;
aa=0.05; 

%%%p560 kenimatic data: with a long tool
%a2=0.0;a3=0.0;d3=0.0;d4=0.0;d6=0.0;        %P560���Ȳ���

%Jaco��е��ģ�Ͳ���
%ĩ��ָ��
% D = [0.2755, 0.4100, 0.2073, 0.07433, 0.07433, 0.1687, 0.0098]; %D1~D6;E2
%ĩ��Ϊ��
D = [0.2755, 0.4100, 0.2073, 0.07433, 0.07433, 0.13, 0.0098]; %D1~D6;E2
pen_long=0.075;                                 %�ʳ�
aaa = 11.0*pi/72;
ca = cos(aaa);
sa = sin(aaa);
c2a = cos(2*aaa);
s2a = sin(2*aaa);
d4b = D(3)+sa/s2a*D(4);
d5b = sa/s2a*D(4)+sa/s2a*D(5);
d6b = sa/s2a*D(5)+D(6);
alpha_i = [pi/2, pi, pi/2, 2*aaa, 2*aaa, pi];
a_i = [0, D(2), 0, 0, 0, 0];
d_i = [D(1), 0, -D(7), -d4b, -d5b, -d6b];

%selection of an initial state q0
qa=[0;0;0;0;0;0];                                           %��е�۳�ʼλ�õ�
qb=[0;-pi/4;0;pi/2;-pi/4;0];
qc=[0;+pi/4;0;pi/2;-pi/4;0];
qd=[0;-pi/4;0;2*pi/3;-pi/4;0];
qe=[0;+pi/4;0;2*pi/3;-pi/4;0];
qf=[0;-pi/4;0;pi/2;-pi/4;0];
qg=[0;-pi/4;0;pi/2;-pi/4;0];

qh=[1.675;1.343;-3.716;4.187;-1.71;1.308];          %jaco������ʼ�ؽڽ�
qi=[1.675;2.843;-3.216;4.187;-1.71;-2.29];          %д�ֵĳ�ʼ�ؽ�
qj=[1.675;2.843;-3.216;4.187;-1.71;-2.65];          %������д�ֵĳ�ʼ�ؽ�
q0=qj;%an initial state q0

global mu_p mu_v qP qM qDp qDm qDDp qDDm myInf n m Pinfty Minfty gamma0 beta Time a2;
n=6;                                    %�ؽ���Ŀ��Ҳ��6���ɶ�
m=3;                                    %ĩ��ִ�����ѿ����ռ�ά�ȣ�Ҳ��3��λ������
gamma0=10;                              %��Ӧ�����е�DNNϵ��alpha
beta=4;%may need adjust
lambda=4;
myInf=1e10;
mu_p=2;
K=5*eye(m,m);
qP=[10000;223;71;10000;10000;10000];    %�ؽڼ��ޣ��ϼ���
qM=[-10000;-137;-289;-10000;-10000;-10000];      %�ؽڼ��ޣ��¼���
rad=180/pi;                             %����ת��
qP=qP/rad;                              %����ת��
qM=qM/rad;                              %����ת��
qDp=ones(n,1)*1.5;qDm=-qDp;%not real, just a test           %�ٶȲ㼫��
% qDDp=ones(6,1)*4;qDDm=-qDDp;%not real, just a test
Pinfty=myInf*ones(m,1);
Minfty=-Pinfty;

Time=0;                                 %��ʼ���ظ��˶�����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fuzzy_control
zz=10e-4;

a = newfis('fuzzy tank');

a = addvar(a,'input','e',[0,10]*zz);
a = addmf(a,'input',1,'ZO','gaussmf',[1.5,0]*zz);
a = addmf(a,'input',1,'PS','gaussmf',[1.5,5]*zz);
a = addmf(a,'input',1,'PB','gaussmf',[1.5,10]*zz);

zzz=500;

a = addvar(a,'output','v',[0,10]*zzz);
a = addmf(a,'output',1,'ZO','trimf',[-5,0,5]*zzz);
a = addmf(a,'output',1,'PS','trimf',[0,5,10]*zzz);
a = addmf(a,'output',1,'PB','trimf',[5,10,15]*zzz);
 
%����ģ������
rulelist=[1 1 1 1;
         2 2 1 1;
         3 3 1 1];
a = addrule(a,rulelist);

a1 = setfis(a,'DefuzzMethod','centroid');
writefis(a1,'tank');
a2 = readfis('tank');

figure(1);
plotfis(a2);
figure(2);
plotmf(a,'input',1);
figure(3);
plotmf(a,'output',1);
  
showrule(a);
ruleview('tank');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





dq0=zeros(n,1);                             %��Ҫ�������ݣ���Ӧ�����е� Theta dot
u0=[dq0;zeros(m,1)];                        %��Ӧ��������Ҫ����u
% u0=(rand(n+m,1)-0.5)*10
init=[q0;u0];%acutally, abitrary initial y0
options=odeset('Mass', @mvn_ZNN_matrixALL, 'MaxStep',0.05,'RelTol',1e-6,'AbsTol',1e-8*ones(n+n+m,1));          %����ODE����ز����������������������;�������
% options=odeset('Mass', @mvn_ZNN_matrixALL, 'RelTol',1e-6,'AbsTol',1e-8*ones(n+n+m,1));          %����ODE����ز����������������������;�������
tic;   %������ʱ
N=1;
j=1;
[t,y]=ode15s(@mvn_ZZJLJrebi1_net,[0,N*T],init,options);%ode15s much better than ode45            %ODE���
cpu_time=toc                                                        %��ȡ��ʱ
size(t)
size(y)
%save cpu_time to "txt" file
fid=fopen('mvn_cost_time.txt','a');                                 %���ļ���¼ʱ��
fprintf(fid,'Computing cost time: %g\n',cpu_time);
fclose(fid);                                                        %�ر��ļ�

%����Ϊ�����ļ�
save mvn_SYSdata u0 aa beta sigma r alpha lambda alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T mu_p mu_v qP qM qDp qDm qDDp qDDm m Pinfty Minfty gamma0 n myInf Time k K a2;
save mvn_INITdata t y;

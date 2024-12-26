format long;

global aa alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T mu_p mu_v qP qM qDp qDm qDDp qDDm m Pinfty Minfty gamma  n K;
load mvn_SYSdata;
load mvn_SNDdata;

Time=0;                             %Time清零
t_last=0;
t_time=0;
[z,n]=size(qAll);
q_Time=zeros(z,1);

[Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q0);          %计算雅可比矩阵
ix0=Ppx(7);iy0=Ppy(7);iz0=Ppz(7);
for jj=1:length(t)
   qjj=qAll(jj,:)';%simulated end effector position
 [Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qjj); 
   j0px(jj,1)=Ppx(1);
   j1px(jj,1)=Ppx(2);
   j2px(jj,1)=Ppx(3);
   j3px(jj,1)=Ppx(4);
   j4px(jj,1)=Ppx(5);
   j5px(jj,1)=Ppx(6);
   posx(jj,1)=Ppx(7);%joint6
   j0py(jj,1)=Ppy(1);
   j1py(jj,1)=Ppy(2);
   j2py(jj,1)=Ppy(3);
   j3py(jj,1)=Ppy(4);
   j4py(jj,1)=Ppy(5);
   j5py(jj,1)=Ppy(6);
   posy(jj,1)=Ppy(7);%joint6
   j0pz(jj,1)=Ppz(1);
   j1pz(jj,1)=Ppz(2);
   j2pz(jj,1)=Ppz(3);
   j3pz(jj,1)=Ppz(4);
   j4pz(jj,1)=Ppz(5);
   j5pz(jj,1)=Ppz(6);
   posz(jj,1)=Ppz(7);%joint6
   %-------------------------------------------------
   dqjj=dqAll(jj,:)';%simulated end effector velocity
  [J,DJ]=ZYNjdj(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qjj,dqjj); 
   %-----------------
   dpos=J*dqjj;
   dposx(jj,1)=dpos(1);
   dposy(jj,1)=dpos(2);
   dposz(jj,1)=dpos(3);
   
   ddqjj=ddqAll(jj,:)';
   ddpos=J*ddqjj+DJ*dqjj;
   ddposx(jj,1)=ddpos(1);
   ddposy(jj,1)=ddpos(2);%simulated end-effector acceleration
   ddposz(jj,1)=ddpos(3);

t_time=t(jj)-Time*T;                    %每个周期内的相对于周期起点的运行时间，对该变量初始化计算出处于本周期中的哪一个时刻
if t_time>T
   q_Time(jj,1)=1;
   Time=Time+1;
   t_time=t_time-T;
end
   


phi=2*pi*(sin(pi*t(jj)/2/T))^2;
phiDot=pi^2/T*sin(pi*t(jj)/T);
phiDotDot=pi^3/T/T*cos(pi*t(jj)/T);
% phi=2*pi/T*(t(jj)-T/2/pi*sin(2*pi*t(jj)/T));
% phiDot=2*pi/T*(1-cos(2*pi*t(jj)/T));
% phiDotDot=4*pi*pi/T/T*sin(2*pi*t(jj)/T);
rx(jj,1)=aa*(2*cos(t(jj))-cos(2*t(jj)))+ix0-0.1;
ry(jj,1)=aa*(2*sin(t(jj))-sin(2*t(jj)))+iy0;
rz(jj,1)=iz0;
drx(jj,1)=aa*(-2*sin(t(jj))+2*sin(2*t(jj)));
dry(jj,1)=aa*(2*cos(t(jj))-2*cos(2*t(jj)));
drz(jj,1)=0;
ddrx(jj,1)=-aa*(-2*cos(t(jj))+4*cos(2*t(jj)));
ddry(jj,1)=-aa*(-2*sin(t(jj))+4*sin(2*t(jj))); 
ddrz(jj,1)=0;


rjj=[rx(jj,1);ry(jj,1);rz(jj,1)];
[fx,fy,fz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qjj);
fjj=[fx(7),fy(7),fz(7)]';
znn_W=norm(qjj-q0)^2*eye(n,n);
znn_Q=[znn_W, J'; J, zeros(m,m)];                 %A(t)
znn_d=[drx(jj,1);dry(jj,1);drz(jj,1)];

ejj=norm(rjj-fjj);

znn_u=[-k*(qjj-q0); znn_d+K*(rjj-fjj)];
% error(jj)=norm(rjj-fjj);
error_function=znn_Q*uAll(jj,:)'-znn_u;
% error(jj)=norm(error_function,2);
if(jj>10)
error(jj)=norm(error_function,2);
else
ejj=0;    
end
t(jj)


end
%--Errors-
if (jj>2)
erposx=rx-posx;
erposy=ry-posy;
erposz=rz-posz;
averposx=mean(erposx);
averposy=mean(erposy);
averposz=mean(erposz);
erdposx=drx-dposx;
erdposy=dry-dposy;
erdposz=drz-dposz;
erddposx=ddrx-ddposx;
erddposy=ddry-ddposy;
erddposz=ddrz-ddposz;
else
erposx=0;
erposy=0;
erposz=0;
averposx=0;
averposy=0;
averposz=0;
erdposx=0;
erdposy=0;
erdposz=0;
erddposx=0;
erddposy=0;
erddposz=0;
end
%均方根
% rmse_erposx=rms(erposx);
% rmse_erposy=rms(erposy);
% rmse_erposz=rms(erposz);
% rmse_erpos=rms(sqrt(erposx.^2+erposy.^2+erposz.^2))

save mvn_SRDdata1 t qAll dqAll ddqAll uAll posx posy posz dposx dposy dposz ddposx ddposy ddposz;
save mvn_SRDdata2 erposx erposy erposz erdposx erdposy  erdposz erddposx erddposy erddposz error averposx averposy averposz;
save mvn_SRDdata3 j0px j1px j2px j3px j4px j5px j0py j1py j2py j3py j4py j5py j0pz j1pz j2pz j3pz j4pz j5pz rx ry rz;

%以下为数据筛选和输出
%等间隔数据筛选
[m,n]=size(qAll);
step_long=0.000005;
jjj=1;
k=1;
q_step(1,:)=qAll(1,:);
for iii=2:m
    deta=(posx(iii)-posx(k))*(posx(iii)-posx(k))+(posy(iii)-posy(k))*(posy(iii)-posy(k))+(posz(iii)-posz(k))*(posz(iii)-posz(k));
    if deta>step_long||q_Time(iii)==1
        jjj=jjj+1;
        k=iii;
        q_step(jjj,:)=qAll(iii,:);
        t_step(jjj)=t(iii);
        q_deta(jjj,1)=deta;
        q_finish(jjj,1)=q_Time(iii);
    end
end
q_deta(1,1)=step_long;
%数据输出，由于Jaco建模和实际的Jaco有初始0度偏角以及旋转方向的的差异，在这里首先对偏角和旋转方向进行修正，再保存数据
q_tran=[2*pi*ones(jjj,1),(1/2*pi)*ones(jjj,1),(3/2*pi)*ones(jjj,1),0*ones(jjj,1),(pi)*ones(jjj,1),0*ones(jjj,1),0*ones(jjj,1),0*ones(jjj,1)];
q_Real=[-q_step(:,1),q_step(:,2),q_step(:,3),q_step(:,4),q_step(:,5),-q_step(:,6), q_deta, q_finish];
q_Real=q_Real+q_tran;
save mvn_Datarun q_Real;


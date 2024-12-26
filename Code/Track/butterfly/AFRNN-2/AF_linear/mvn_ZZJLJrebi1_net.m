function doty=mvn_ZZJLJrebi1_net(t,y)

global aa alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T m gamma0 n Time k K a2 dr znngamma;

persistent ix0 iy0 iz0 ix1 iy1 iz1 t_last;

% figure(1);
% plotfis(a2);



   q=y(1:n);                         %读取初始状
   dq=y(n+1:n+n);                        %读取Theta dot，对应文章中的Theta dot
   yy=y(n+1:n+n+m);


if t==0
    [x,y,z]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q0);  %计算初始位置
    ix0=x(7);iy0=y(7);iz0=z(7);
end

t_time=t-Time*T;                    %每个周期内的相对于周期起点的运行时间，对该变量初始化计算出处于本周期中的哪一个时刻

if t_time>T
   Time=Time+1;
   t_time=t_time-T;
end

phi_sin=2*pi*sin(0.5*pi*t/T);
phi=phi_sin*sin(0.5*pi*t/T);
fac=exp(cos(phi))-2*cos(4*phi)-sin(phi/12)^5;%factor

rx=aa*sin(phi)*fac+ix0;%*T/2/pi;
ry=aa*(cos(phi)*fac)+iy0-0.07182818;%T/2/pi;
rz=0+iz0;

drx=-aa*sin(2*pi*sin((pi*t)/(2*T))^2)*((5*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(8*pi*sin((pi*t)/(2*T))^2))/T + (2*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T) - (2*pi^2*aa*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T;
dry=(2*pi^2*aa*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T - aa*cos(2*pi*sin((pi*t)/(2*T))^2)*((5*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(8*pi*sin((pi*t)/(2*T))^2))/T + (2*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T);
drz=0;

ddrx=aa*sin(2*pi*sin((pi*t)/(2*T))^2)*((8*pi^3*cos((pi*t)/(2*T))^2*sin(8*pi*sin((pi*t)/(2*T))^2))/T^2 - (8*pi^3*sin((pi*t)/(2*T))^2*sin(8*pi*sin((pi*t)/(2*T))^2))/T^2 + (5*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*sin((pi*sin((pi*t)/(2*T))^2)/6)^5)/(36*T^2) - (5*pi^3*cos((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(12*T^2) - (pi^3*cos((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T^2 + (5*pi^3*sin((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(12*T^2) + (pi^3*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T^2 + (128*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*cos(8*pi*sin((pi*t)/(2*T))^2))/T^2 - (5*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)^2*sin((pi*sin((pi*t)/(2*T))^2)/6)^3)/(9*T^2) + (4*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2)^2)/T^2 - (4*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*cos(2*pi*sin((pi*t)/(2*T))^2))/T^2) - (pi^3*aa*cos((pi*t)/(2*T))^2*cos(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2 + (pi^3*aa*sin((pi*t)/(2*T))^2*cos(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2 - (4*pi^2*aa*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos(2*pi*sin((pi*t)/(2*T))^2)*((5*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(8*pi*sin((pi*t)/(2*T))^2))/T + (2*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T))/T + (4*pi^4*aa*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*sin(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2;
ddry=aa*cos(2*pi*sin((pi*t)/(2*T))^2)*((8*pi^3*cos((pi*t)/(2*T))^2*sin(8*pi*sin((pi*t)/(2*T))^2))/T^2 - (8*pi^3*sin((pi*t)/(2*T))^2*sin(8*pi*sin((pi*t)/(2*T))^2))/T^2 + (5*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*sin((pi*sin((pi*t)/(2*T))^2)/6)^5)/(36*T^2) - (5*pi^3*cos((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(12*T^2) - (pi^3*cos((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T^2 + (5*pi^3*sin((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(12*T^2) + (pi^3*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T^2 + (128*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*cos(8*pi*sin((pi*t)/(2*T))^2))/T^2 - (5*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*cos((pi*sin((pi*t)/(2*T))^2)/6)^2*sin((pi*sin((pi*t)/(2*T))^2)/6)^3)/(9*T^2) + (4*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2)^2)/T^2 - (4*pi^4*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*cos(2*pi*sin((pi*t)/(2*T))^2))/T^2) + (pi^3*aa*cos((pi*t)/(2*T))^2*sin(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2 - (pi^3*aa*sin((pi*t)/(2*T))^2*sin(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2 + (4*pi^4*aa*cos((pi*t)/(2*T))^2*sin((pi*t)/(2*T))^2*cos(2*pi*sin((pi*t)/(2*T))^2)*(sin((pi*sin((pi*t)/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t)^2)) + 2*cos(8*pi*sin((pi*t)/(2*T))^2)))/T^2 + (4*pi^2*aa*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(2*pi*sin((pi*t)/(2*T))^2)*((5*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*cos((pi*sin((pi*t)/(2*T))^2)/6)*sin((pi*sin((pi*t)/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*sin(8*pi*sin((pi*t)/(2*T))^2))/T + (2*pi^2*cos((pi*t)/(2*T))*sin((pi*t)/(2*T))*exp(cos(2*pi*sin((pi*t)/(2*T))^2))*sin(2*pi*sin((pi*t)/(2*T))^2))/T))/T;
ddrz=0;

deltD=delt_D(t);
delts=delt_s(t);

dr=[drx;dry;drz];    %对应的是文章中的r_dot
ddr=[ddrx;ddry;ddrz];

[fx,fy,fz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q);

[J,DJ]=ZYNjdj(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q,dq); 
df=J*dq;

f=[fx(7);fy(7);fz(7)];
r=[rx;ry;rz];

nerrEC=norm(dr-df);
nerrZnn=norm(r-f);

znngamma=znn_gamma(nerrEC,nerrZnn);
znn_d=dr+K*(r-f);
znn_W=norm(q-q0)^2*eye(n,n);
znn_Q=[znn_W, J'; J, zeros(m,m)];                 %A(t)
znn_dotQ=[2*(q-q0)'*dq*eye(n,n), DJ'; DJ, zeros(m,m)];       %\dot A(t)


k=5;                                                %对应文章中的c等式中的gamma，这里取0表示不考虑重复运动的约束
znn_q=k*(q-q0);                                         %对应文章中的c
znn_u=[-znn_q; znn_d];
znn_dotu=[-k*dq; ddr+K*(dr-J*dq)];

dotyy=-(znn_dotQ+deltD)*yy-znngamma*AFMlinear(znn_Q*yy-znn_u)+znn_dotu+delts;
% dotyy=-(znn_dotQ+deltD)*yy-gamma0*exp(t)*AFMlinear(znn_Q*yy-znn_u)+znn_dotu+delts;
t
doty=[dq;dotyy;];

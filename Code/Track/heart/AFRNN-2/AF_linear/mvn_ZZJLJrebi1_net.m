function doty=mvn_ZZJLJrebi1_net(t,y)

global alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T m gamma0 n k K aa;

persistent ix0 iy0 iz0 ix1 iy1 iz1 t_last;

% figure(1);
% plotfis(a2);



   q=y(1:n);                         %读取初始状
   dq=y(n+1:n+n);                        %读取Theta dot，对应文章中的Theta dot
   yy=y(n+1:n+n+m);


if t==0
    [x,y,z]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q0);  %计算初始位置
    ix0=x(7);iy0=y(7);iz0=z(7);
%     disp(ix0);
%     disp(iy0);
%     disp(iz0);
end

% t_time=t-Time*T;                    %每个周期内的相对于周期起点的运行时间，对该变量初始化计算出处于本周期中的哪一个时刻
% 
% if t_time>T
%    Time=Time+1;
%    t_time=t_time-T;
% end
rx=aa*(2*cos(t)-cos(2*t))+ix0-0.05;
ry=aa*(2*sin(t)-sin(2*t))+iy0;
rz=iz0;
drx=aa*(-2*sin(t)+2*sin(2*t));
dry=aa*(2*cos(t)-2*cos(2*t));
drz=0;
ddrx=-aa*(-2*cos(t)+4*cos(2*t));
ddry=-aa*(-2*sin(t)+4*sin(2*t)); 
ddrz=0;
dr=[drx;dry;drz];

[fx,fy,fz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q);
f=[fx(7);fy(7);fz(7)];
r=[rx;ry;rz];
% disp(rx);
% disp(ry);
% disp(rz);
% disp(fx(7));
% disp(fy(7));
% disp(fz(7));

[J,DJ]=ZYNjdj(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q,dq); 
df=J*dq;
nerrEC=norm(dr-df);
nerrZnn=norm(r-f);
fis_a = readfis('fis2');
fis_out=evalfis([nerrZnn,nerrEC],fis_a);
v=norm(fis_out,2);
gamma = gamma0 +v;



% gamma = gamma0;

dr=[drx;dry;drz];           %对应的是文章中的r_dot
ddr=[ddrx;ddry;ddrz];
znn_d=dr+K*(r-f);

% znn_W=eye(n,n);
% znn_Q=[znn_W, J'; J, zeros(m,m)];
% znn_dotQ=[zeros(n,n), DJ'; DJ, zeros(m,m)];

%znn_W=eye(n,n);
znn_W=norm(q-q0)^2*eye(n,n);
znn_Q=[znn_W, J'; J, zeros(m,m)];                 %A(t)
znn_dotQ=[2*(q-q0)'*dq*eye(n,n), DJ'; DJ, zeros(m,m)];       %\dot A(t)


k=5;                                                %对应文章中的c等式中的gamma，这里取0表示不考虑重复运动的约束
% zz=k*(q-q0);
znn_q=k*(q-q0);                                         %对应文章中的c
znn_u=[-znn_q; znn_d];
znn_dotu=[-k*dq; ddr+K*(dr-J*dq)];


%%%%%%%%%%%%%%%%%扰动
% if t>0.3*T && t<0.8*T
if t>0
delt_D=0.2*[sin(t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(3*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             cos(2*t/T*pi) cos(2*t/T*pi) -cos(4*t/T*pi) sin(t/T*pi) -sin(2*t/T*pi) -cos(3*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(3*t/T*pi) -sin(4*t/T*pi) cos(5*t/T*pi) -cos(3*t/T*pi) cos(4*t/T*pi) -sin((5*t/T*pi)) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(2*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             -sin(2*t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(2*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi)];
delt_s=0.2*[2*sin(t/T*pi);4*cos(2*t/T*pi);-5*sin(3*t/T*pi); cos(3*t/T*pi); -sin(3*t/T*pi); 3*cos(t/T*pi);-sin(t/T*pi);-sin(2*t/T*pi);cos(t/T*pi)];
else 
    delt_D=0;
    delt_s=0;
end

dotyy=-(znn_dotQ+delt_D)*yy-gamma*AFMlinear(znn_Q*yy-znn_u)+znn_dotu+delt_s;
t
doty=[dq;dotyy];

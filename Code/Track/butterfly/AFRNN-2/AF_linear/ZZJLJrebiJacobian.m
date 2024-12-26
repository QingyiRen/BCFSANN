function Jacob0=ZZJLJrebiJacobian(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q)
	mys=sin(q);
    s1=mys(1);s2=mys(2);s3=mys(3);s4=mys(4);s5=mys(5);s6=mys(6);
    myc=cos(q);
    c1=myc(1);c2=myc(2);c3=myc(3);c4=myc(4);c5=myc(5);c6=myc(6);
    
    alpha_4 = alpha_i(4);alpha_5 = alpha_i(5);
    d_1 = d_i(1);d_2 = d_i(2);d_3 = d_i(3);d_4 = d_i(4);
    d_5 = d_i(5);d_6 = d_i(6);
%     sp = sin(phi);
%     cp = cos(phi);
    sa4 = sin(alpha_4);ca4 = cos(alpha_4);
    sa5 = sin(alpha_5);ca5 = cos(alpha_5);
    
%以下的AW1和AW2实际上是由于机械臂的安装或建模问题最终引起了画图的一点小小误差，通过基点的一个旋转矩阵偏移一个小角度让笔尖画出了的轨迹在平面内
    sp1=sin(0.02);cp1=cos(0.02);
    sp2=sin(0.05);cp2=cos(0.05);
    
    AW1 = [1, 0, 0, 0;...
    0, cp1, sp1, 0;...
    0, -sp1, cp1, 0;...
    0, 0, 0, 1];
    
    AW2 = [cp2, 0, -sp2, 0;...
    0, 1, 0, 0;...
    sp2, 0, cp2, 0;...
    0, 0, 0, 1];

    AW1=AW1*AW2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%以下为Jaco机械臂旋转矩阵
    A1 = [c1, 0, s1, 0;...
    s1, 0, -c1, 0;...
    0, 1, 0, d_1;...
    0, 0, 0, 1];
    DA1 = [-s1, 0, c1, 0;...
    c1, 0, s1, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA1 = [-c1, 0, -s1, 0;...
    -s1, 0, c1, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    
    A2 = [c2, s2, 0, D(2)*c2;...
    s2, -c2, 0, D(2)*s2;...
    0, 0, -1, 0;...
    0, 0, 0, 1];
    DA2 = [-s2, c2, 0, D(2)*-s2;...
    c2, s2, 0, D(2)*c2;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA2 = [-c2, -s2, 0, D(2)*-c2;...
    -s2, c2, 0, D(2)*-s2;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    
    A3 = [c3, 0, s3, 0;...
    s3, 0, -c3, 0;...
    0, 1, 0, d_3;...
    0, 0, 0, 1];
    DA3 = [-s3, 0, c3, 0;...
    c3, 0, s3, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA3 = [-c3, 0, -s3, 0;...
    -s3, 0, c3, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0]; 
    
    A4 = [c4, -ca4*s4, sa4*s4, 0;...
    s4, ca4*c4, -sa4*c4, 0;...
    0, sa4, ca4, d_4;...
    0, 0, 0, 1];
    DA4 = [-s4, -ca4*c4, sa4*c4, 0;...
    c4, ca4*-s4, sa4*s4, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA4 = [-c4, ca4*s4, sa4*-s4, 0;...
    -s4, ca4*-c4, sa4*c4, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    
    A5 = [c5, -ca5*s5, sa5*s5, 0;...
    s5, ca5*c5, -sa5*c5, 0;...
    0, sa5, ca5, d_5;...
    0, 0, 0, 1];
    DA5 = [-s5, -ca5*c5, sa5*c5, 0;...
    c5, ca5*-s5, sa5*s5, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA5 = [-c5, ca5*s5, sa5*-s5, 0;...
    -s5, ca5*-c5, sa5*c5, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];

    A6 = [c6, s6, 0, 0;...
    s6, -c6, 0, 0;...
    0, 0, -1, d_6;...
    0, 0, 0, 1];
    DA6 = [-s6, c6, 0, 0;...
    c6, s6, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];
    DDA6 = [-c6, -s6, 0, 0;...
    -s6, c6, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0];

%PUMA560旋转矩阵
%     A1=[c1,0,s1,0;...
% 		s1,0,-c1,0;...
% 		0,1,0,0;...
% 		0,0,0,1];
% 	DA1=[-s1,0,c1,0;...
% 		 c1,0,s1,0;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA1=[-c1,0,-s1,0;...
% 		  -s1,0,c1,0;...
% 		  0,0,0,0;...
% 		  0,0,0,0];
% 	A2=[c2,-s2,0,a2*c2;...
% 		s2,c2,0,a2*s2;...
% 		0,0,1,0;...
% 		0,0,0,1];
% 	DA2=[-s2,-c2,0,-a2*s2;...
% 		 c2,-s2,0,a2*c2;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA2=[-c2,s2,0,-a2*c2;...
% 		  -s2,-c2,0,-a2*s2;...
% 		  0,0,0,0;...
% 		  0,0,0,0];
% 	A3=[c3,0,-s3,c3*a3;...
% 		s3,0,c3,s3*a3;...
% 		0,-1,0,d3;...
% 		0,0,0,1];
% 	DA3=[-s3,0,-c3,-s3*a3;...
% 		 c3,0,-s3,c3*a3;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA3=[-c3,0,s3,-c3*a3;...
% 		  -s3,0,-c3,-s3*a3;...
% 		  0,0,0,0;...
% 		  0,0,0,0];
% 	A4=[c4,0,s4,0;...
% 		s4,0,-c4,0;...
% 		0,1,0,d4;...
% 		0,0,0,1];
% 	DA4=[-s4,0,c4,0;...
% 		 c4,0,s4,0;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA4=[-c4,0,-s4,0;...
% 		  -s4,0,c4,0;...
% 		  0,0,0,0;...
% 		  0,0,0,0];
% 	A5=[c5,0,-s5,0;...
% 		s5,0,c5,0;...
% 		0,-1,0,0;...
% 		0,0,0,1];
% 	DA5=[-s5,0,-c5,0;...
% 		 c5,0,-s5,0;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA5=[-c5,0,s5,0;...
% 		  -s5,0,-c5,0;...
% 		  0,0,0,0;...
% 		  0,0,0,0];
% 	A6=[c6,-s6,0,0;...
% 		s6,c6,0,0;...
% 		0,0,1,d6;...
% 		0,0,0,1];
% 	DA6=[-s6,-c6,0,0;...
% 		 c6,-s6,0,0;...
% 		 0,0,0,0;...
% 		 0,0,0,0];
% 	DDA6=[-c6,s6,0,0;...
% 		  -s6,-c6,0,0;...
% 		  0,0,0,0;...
% 		  0,0,0,0];

    Pen=[0.0;-pen_long;0.0;1.0];

	A6v=A6*Pen;
    
	A56v=A5*A6v;
	A456v=A4*A56v;
	A3456v=A3*A456v;
	A23456v=A2*A3456v;
	Je1=DA1*A23456v;
	D2A3456v=DA2*A3456v;
    Je2=A1*D2A3456v;
    D3A456v=DA3*A456v;
    A2D3A456v=A2*D3A456v;
    Je3=A1*A2D3A456v;
    D4A56v=DA4*A56v;
    A3D4A56v=A3*D4A56v;
    A23D4A56v=A2*A3D4A56v;
    Je4=A1*A23D4A56v;
    D5A6v=DA5*A6v;
    A4D5A6v=A4*D5A6v;
    A34D5A6v=A3*A4D5A6v;
    A234D5A6v=A2*A34D5A6v;
    Je5=A1*A234D5A6v;
    D6v=DA6*Pen;
    A5D6v=A5*D6v;
    A45D6v=A4*A5D6v;
    A345D6v=A3*A45D6v;
    A2345D6v=A2*A345D6v;
    Je6=A1*A2345D6v;
    Je1=AW1*Je1;
    Je2=AW1*Je2;
    Je3=AW1*Je3;
    Je4=AW1*Je4;
    Je5=AW1*Je5;
    Je6=AW1*Je6;
    
    Jacob0(:,1)=Je1(1:3,1);
    Jacob0(:,2)=Je2(1:3,1);
    Jacob0(:,3)=Je3(1:3,1);
    Jacob0(:,4)=Je4(1:3,1);
    Jacob0(:,5)=Je5(1:3,1);
    Jacob0(:,6)=Je6(1:3,1);
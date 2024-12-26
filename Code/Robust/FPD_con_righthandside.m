function  y=FPD_con_righthandside(t,x,gamma)
if nargin==2, gamma=1;end
diffW=[-sin(t) 2*cos(2*t) -3*sin(3*t)
        2*cos(2*t) -sin(t) 3*cos(3*t)
        -3*sin(3*t) 3*cos(3*t) 0]; 
diffu=[2*sin(2*t)
       -2*cos(2*t) 
       4*cos(4*t)];
matrixD=[ 5*sin(1*t) 4*cos(1*t) -3*sin(1*t)
   cos(2*t)   2*cos(2*t) -3*cos(2*t)
    sin(3*t)  -1*sin(4*t) cos(5*t) ];
vectorS=[ 2*sin(t)
    4*cos(2*t)
   -5*sin(3*t)];
y=-(diffW+matrixD)*x-gamma*AFlinear((LIVEmatrixW(t,x)*x-LIVEvectorU(t,x))+diffu+vectorS;


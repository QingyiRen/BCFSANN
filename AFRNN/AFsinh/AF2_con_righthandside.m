function y=AF2_con_righthandside(t,x,gamma)
if nargin==2, gamma=1;end
diffW=[-sin(t) 2*cos(2*t) -3*sin(3*t)
        2*cos(2*t) -sin(t) 3*cos(3*t)
        -3*sin(3*t) 3*cos(3*t) 0]; 
diffu=[2*sin(2*t)
       -2*cos(2*t) 
       4*cos(4*t)];
err=LIVEmatrixW(t,x)*x-LIVEvectorU(t,x);
EC=diffW*x-diffu;
nerrEC=norm(EC);
nerrZnn=norm(err);
fis_a = readfis('fis2_1');
fis_out=evalfis([nerrZnn,nerrEC],fis_a);
v=norm(fis_out,2);
y=-diffW*x-(gamma+v)*AFsinh(LIVEmatrixW(t,x)*x-LIVEvectorU(t,x))+diffu;
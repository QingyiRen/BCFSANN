function y=AF_con_righthandside(t,x,gamma)
if nargin==2, gamma=1;end
diffW=[-sin(t) 2*cos(2*t) -3*sin(3*t)
        2*cos(2*t) -sin(t) 3*cos(3*t)
        -3*sin(3*t) 3*cos(3*t) 0]; 
diffu=[2*sin(2*t)
       -2*cos(2*t) 
       4*cos(4*t)];
err=LIVEmatrixW(t,x)*x-LIVEvectorU(t,x);
nerrZnn=norm(err);
fis_a = readfis('fis');
fis_out=evalfis(nerrZnn,fis_a);
v=norm(fis_out,2);
y=-diffW*x-(gamma+v)*AFtunable(LIVEmatrixW(t,x)*x-LIVEvectorU(t,x))+diffu;

 


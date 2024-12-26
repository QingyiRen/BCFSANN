function a2=fis2()
z=10e-4;
a = newfis('fis2');

a = addvar(a,'input','e',[0,10]*z);
a = addmf(a,'input',1,'NB','gaussmf',[1.5,0]*z);
a = addmf(a,'input',1,'NM','gaussmf',[1.5,1.667]*z);
a = addmf(a,'input',1,'NS','gaussmf',[1.5,3.333]*z);
a = addmf(a,'input',1,'Z','gaussmf',[1.5,5]*z);
a = addmf(a,'input',1,'PS','gaussmf',[1.5,6.666]*z);
a = addmf(a,'input',1,'PM','gaussmf',[1.5,8.334]*z);
a = addmf(a,'input',1,'PB','gaussmf',[1.5,10]*z);

zz=0.001;
a = addvar(a,'input','de',[0,10]*zz);
a = addmf(a,'input',2,'NB','gaussmf',[1.5,0]*zz);
a = addmf(a,'input',2,'NM','gaussmf',[1.5,1.667]*zz);
a = addmf(a,'input',2,'NS','gaussmf',[1.5,3.333]*zz);
a = addmf(a,'input',2,'Z','gaussmf', [1.5,5]*zz);
a = addmf(a,'input',2,'PS','gaussmf',[1.5,6.666]*zz);
a = addmf(a,'input',2,'PM','gaussmf',[1.5,8.334]*zz);
a = addmf(a,'input',2,'PB','gaussmf',[1.5,10]*zz);


zzz=1000;

a = addvar(a,'output','v',[0,10]*zzz);
a = addmf(a,'output',1,'NB','gaussmf',[1.0,0]*zzz);
a = addmf(a,'output',1,'NM','gaussmf',[1.0,1.667]*zzz);
a = addmf(a,'output',1,'NS','gaussmf',[1.0,3.333]*zzz);
a = addmf(a,'output',1,'Z','gaussmf', [1.0,5]*zzz);
a = addmf(a,'output',1,'PS','gaussmf',[1.0,6.666]*zzz);
a = addmf(a,'output',1,'PM','gaussmf',[1.0,8.334]*zzz);
a = addmf(a,'output',1,'PB','gaussmf',[1.0,10]*zzz);
 
%建立模糊规则
rulelist=[
        1 1 7 1 1;
        1 2 7 1 1;
        1 3 6 1 1;
        1 4 6 1 1;
        1 5 5 1 1;
        1 6 4 1 1;
        1 7 4 1 1;
        2 1 7 1 1;
        2 2 7 1 1;
        2 3 6 1 1;
        2 4 5 1 1;
        2 5 5 1 1;
        2 6 4 1 1;
        2 7 4 1 1;
        3 1 6 1 1;
        3 2 5 1 1;
        3 2 6 1 1;
        3 3 6 1 1;
        3 4 5 1 1;
        3 5 4 1 1;
        3 6 3 1 1;
        3 7 2 1 1;
        4 1 6 1 1;
        4 2 5 1 1;
        4 3 5 1 1;
        4 4 4 1 1;
        4 5 3 1 1;
        4 6 2 1 1;
        4 7 2 1 1;
        5 1 5 1 1;
        5 2 5 1 1;
        5 3 4 1 1;
        5 4 3 1 1;
        5 5 3 1 1;
        5 6 2 1 1;
        5 7 2 1 1;
        6 1 4 1 1;
        6 2 4 1 1;
        6 3 3 1 1;
        6 5 2 1 1;
        6 6 2 1 1;
        6 7 1 1 1;
        7 1 4 1 1;
        7 2 3 1 1;
        7 3 3 1 1;
        7 4 2 1 1;
        7 5 2 1 1;
        7 6 1 1 1;
        7 7 1 1 1;];
a = addrule(a,rulelist);
a1 = setfis(a,'DefuzzMethod','centroid');
writefis(a1,'fis2');
a2 = readfis('fis2');


x1=0:10e-6:0.001;
y1=0:10e-5:0.01;
for ii=1:length(x1)
    for jj=1:length(y1)
     fis_out=evalfis([x1(ii),y1(ii)],a2);
     Z(ii,jj)=norm(fis_out,2);
    end
end


% [X,Y]=meshgrid(x1,y1);
% mesh(X,Y,Z);
figure(1);
plotmf(a,'input',1);
figure(2);
plotmf(a,'input',2);
figure(3);
plotmf(a,'output',1);



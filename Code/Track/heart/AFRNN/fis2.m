function a2=fis2()
 zz=10e-1;

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
 
%建立模糊规则
rulelist=[1 1 1 1;
         2 2 1 1;
         3 3 1 1];
a = addrule(a,rulelist);

a1 = setfis(a,'DefuzzMethod','centroid');
writefis(a1,'fis');
a2 = readfis('fis');

figure(1);
plotfis(a2);
figure(2);
plotmf(a,'input',1);
figure(3);
plotmf(a,'output',1);
  
showrule(a);
ruleview('fis');



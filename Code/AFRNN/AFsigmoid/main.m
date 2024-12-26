tspan=[0 8]; 
gamma=1;
for iter=1:6
y0=1*(0.5*rand(3,1)-ones(3,1));
AFRNNerror(y0,tspan,gamma);
% legend('FP-CDNN');
hold on;
AFRNN2error(y0,tspan,gamma);
% legend('AFRNN-2');
end
legend('AFRNN','FPDRNN');


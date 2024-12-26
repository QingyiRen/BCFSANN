tspan=[0 4]; 
%If gamma is 1, then the end-time is 10s for creating figures
%If gamma is 10,then the end-time is 1s  for creating figures
%If gamma is 100,then the end-time is 0.1s for creating figures
options=odeset('Mass', @LIVEmatrixW, 'MStateDep', 'none');
gamma=1;
for iter=1:6
    y0=1*(0.5*rand(3,1)-ones(3,1));
    [t,y]=ode45(@AF2_con_righthandside, tspan, y0, options, gamma);
         total=length(t);
    for i=1:total
    ts=t(i);
    W=LIVEmatrixW(ts,y);
    u=LIVEvectorU(ts,y);
    A=W*(y(i,:))';
    lambdastar(i,:)=u(1,1);
    lambda(i,:)=A(1,1);
    end 
    plot(t,lambda(1:length(t)),'r-','linewidth',1.5);
    xlabel('t(s)');
    legend('x_{1}(t) ');                                          
    hold on;
end
    plot(t,lambdastar(1:length(t)),'b--','linewidth',1.5);
    
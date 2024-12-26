format long;
load mvn_SYSdata;
load mvn_INITdata;
% global aa q0 T mu_p mu_v qP qM qDp qDm qDDp qDDm m Pinfty Minfty gamma myInf n;

interval=0.005;                         %最小步长
jj=0;
epsilon=0;
for ii=1:length(t)
    if(t(ii,1)>=epsilon)                %比累计最小步长大？这一步是用来干嘛的？
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        yn(jj,:)=y(ii,:);
        %%%%%%%%%%%%%%%%%
        doty=mvn_ZZJLJrebi1_net(tn(jj,1),yn(jj,:)');     %重新代入计算该时刻值神经网络输出
        ddqAll(jj,:)=doty((n+1):(n+n),1)';              %取出u_dot的前n个参数
        %%%%%%%%%%%%%%%%%
        epsilon=jj*interval;        
    elseif(ii==length(t))
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        yn(jj,:)=y(ii,:);
        %%%%%%%%%%%%%%%%%
        doty=mvn_ZZJLJrebi1_net(tn(jj,1),yn(jj,:));
        ddqAll(jj,:)=doty((n+1):(n+n),1)';
        %%%%%%%%%%%%%%%%%
        epsilon=jj*interval;
    end
end
clear t y;
t=tn;
y=yn;
size(t)
size(y)
clear tn yn;
qAll=y(:,1:n);                                      %取出计算的Theta
size(qAll)
uAll=y(:,(n+1):(n+n+m));                            %取出计算的u
size(uAll)
dqAll=uAll(:,1:n);                                  %取出u的前n个
size(dqAll)
size(ddqAll)
save mvn_SNDdata t qAll dqAll uAll ddqAll;

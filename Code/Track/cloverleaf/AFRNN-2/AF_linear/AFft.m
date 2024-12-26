function y=AFft(E)
if nargin==1, k1=1;k2=1;k3=1;r1=0.5;
end
y=sign(E).*(1*(abs(E)).^(0.5)+1*(abs(E))+1*(abs(E)).^(1/0.5));
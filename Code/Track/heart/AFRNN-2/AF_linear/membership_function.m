x = 0:0.01:10;

y1= trimf(x, [1 5 9]);
figure(1)
plot(x,y1,'LineWidth',2);

y2= gbellmf(x, [2 4 5]);
figure(2)
plot(x,y2,'LineWidth',2);

y3 = gaussmf(x,[1,5]);
figure(3)
plot(x,y3,'LineWidth',2);

y4 = sigmf(x,[1.5,5]);
figure(4)
plot(x,y4,'LineWidth',2);
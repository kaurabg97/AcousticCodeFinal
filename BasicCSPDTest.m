x = linspace(1,100,1000);
y = linspace(2,101,1000);
Xw = sin(x);
Yw = cos(x+pi/3);
figure;
plot(x,Xw)
hold on
plot(x,Yw)

%%
figure;
[croSpecB1,FrB1] = cpsd(Xw,Yw,[],[],500000,100000);
semilogx((FrB1),abs(croSpecB1),'LineWidth',1.2,'Color','r')
phsDifB = angle(croSpecB1);

figure;
semilogx((FrB1),phsDifB,'o','Color','k')
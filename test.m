
x = linspace(0,10,100)
%x1 = exp(x)
x1 = sin(x)
figure(1)
subplot(2,1,1)
plot(x,x1)

[Rxx lags] = xcorr(x1,x1) %v lags means Tau's
subplot(2,1,2)
plot(lags, Rxx)

%% Ccross corelation

x = linspace(0,10,100)
x1 = sin(x)
x2 = sin(x+pi)
figure(1)
title('Phase diff Pi')
subplot(2,1,1)
plot(x,x1)

[Rxx lags] = xcorr(x1,x2) %v lags means Tau's
loc = max(Rxx)
loc2 = find(Rxx == loc)
samplelag = (loc2)
subplot(2,1,2)
title('Phase diff Pi')
plot(lags, Rxx)
 X = fft(x1)
 stem(abs(X))

 %% CSPD calculation


x = linspace(0,10,100)
x1 = sin(x)
x2 = sin(x+pi/6)
[CSPD,Fr] = cpsd(x1,x2,[],[],length(x));
phsDifB1 = angle(CSPD)
figure(1)
subplot(211)
semilogx(Fr,abs(CSPD),'LineWidth',1.2,'Color','r')
subplot(212)
semilogx(Fr,phsDifB1/pi,'o','Color','r')




clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a = 1; %Test 2,0.5,10
c = 1;
varv = 1;
varw = 1;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);
varx = 1;

%time passing
for n=1:100
    w = normrnd(0,varw);
    x = [x a*x(1,n)+w];
    v = normrnd(0,varv);
    y = [y c*x(1,(n+1))+v];
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,n+1)*c*varx/vary];
end

%hold all
subplot(1,2,1), plot(0:n,x,'k')
subplot(1,2,1), hold on, plot(0:n,y,'b')
subplot(1,2,1), hold on, plot(0:n,xhat,'r')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

subplot(1,2,2), plot(0:n,y/c-x,'b')
subplot(1,2,2), hold on, plot(0:n,xhat-x,'r')
xlabel('n=time')
legend2 = legend('y/c-x','xhat-x');

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

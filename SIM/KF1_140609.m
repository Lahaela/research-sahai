clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a = 2; %Test 2,0.5,10
c = 100;
varv = 10; %channel noise
varw = 1; %system noise
n = 25;

%initialize
x = zeros(1,n);
y = zeros(1,n);
xhat = zeros(1,n);
xhatmem = zeros(1,n);

x(1,1) = normrnd(0,1); %x[0], matlab is 1-indexed
y(1,1) = c*x(1,1)+normrnd(0,varv);
xhat(1,1) = y(1,1)*c*1/(c^2+varv);
xhatmem(1,1) = 0;
varx = [1];
vary = [c^2+varv];

%time passing
for t=1:(n-1)
    x(1,t+1) = a*x(1,t)+normrnd(0,varw);
    y(1,t+1) = c*x(1,(t+1))+normrnd(0,varv);
    varx = [varx a^2*varx + varw];
    vary = [vary c^2*varx + varv];
    xhat(1,t+1) = y(1,t+1)*c*varx(1,t+1)/vary(1,t+1);
    star = y(1,t+1)*a^4*c^3*varx(1,t)^2*(vary(1,t+1)-c^2*varx(1,t+1))/(a^2*vary(1,t+1)^2*c^2*varx(1,t)+a^2*vary(1,t+1)^2*varv-c^4*varx(1,t)^2*(c^2*a^2*varx(1,t)+c^2*varw+varv));
    xhatmem(1,t+1) = -star+xhat(1,t+1);
    t
end

%hold all
subplot(1,2,1), stem(0:t,x,'k','filled')
subplot(1,2,1), hold on, stem(0:t,y,'b')
subplot(1,2,1), hold on, stem(0:t,xhat,'r')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

subplot(1,2,2), stem(0:t,y/c-x,'b')
subplot(1,2,2), hold on, stem(0:t,xhat-x,'r','filled')
subplot(1,2,2), hold on, stem(0:t,xhatmem-x, 'g','filled')
xlabel('n=time')
legend2 = legend('y/c-x','xhat-x');

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])
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
    star = (a^2*c*varx(1,t)*vary(1,t)*(vary(1,t)-c^2*varx(1,t)))/(a^2*c^2*varx(1,t)*(c^2*varx(1,t)*(varx(1,t)+varv)-2*c*varx(1,t)*vary(1,t)+vary(1,t)^2)-2*a*c^2*varv*varx(1,t)*vary(1,t)+vary(1,t)^2*(c^2*varw+varv));
    star = star*(y(1,t+1)-c^2*a*varx(1,t)/vary(1,t)*y(1,t));
    xhatmem(1,t+1) = star+(a*c*varx(1,t)/vary(1,t)*y(1,t));
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

%%

clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a = 2; %Test 2,0.5,10
c = 1;
varv = 10;
varw = 1;
n = 25;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);
xhatmem = y*c/(c^2+varv);
varx = 1;
vary = c^2+varv;

%time passing
for t=1:(n-1)
    w = normrnd(0,varw);
    x = [x a*x(1,t)+w];
    v = normrnd(0,varv);
    y = [y c*x(1,(t+1))+v];
    star = (a^2*c*varx*vary*(vary-c^2*varx))/(a^2*c^2*varx*(c^2*varx*(varx+varv)-2*c*varx*vary+vary^2)-2*a*c^2*varv*varx*vary+vary^2*(c^2*varw+varv));
    star = star*(y(1,t+1)-c^2*a*varx/vary*y(1,t));
    xhatmem = [xhatmem star+(a*c*varx/vary*y(1,t))];
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
end

%hold all
subplot(1,2,1), stem(0:t,x,'k')
subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, stem(0:t,xhat,'r')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

subplot(1,2,2), stem(0:t,abs(y/c-x),'b','filled')
subplot(1,2,2), hold on, stem(0:t,abs(xhat-x),'r','filled')
subplot(1,2,2), hold on, stem(0:t,abs(xhatmem-x),'g','filled')
xlabel('n=time')
legend2 = legend('y/c-x','xhat-x');

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a1 = 1; %Test 2,0.5,10
a2 = 5;
c1 = 1;
c2 = 1;
n = 10;

%initialize
x1 = normrnd(0,1,[1,1]); %x[0]
x2 = normrnd(0,1,[1,1]);
y = c1*x1+c2*x2;

varx1 = 1;
varx2 = 1;
vary = c1^2+c2^2;

xhat1 = y*c1*varx1/vary;
xhat2 = y*c2*varx2/vary;
sigman1 = (1-c1*varx1/vary)^2*varx1+(c2*varx1/vary)^2*varx2;

%time passing
for t=1:(n-1)
    x1 = [x1 a1*x1(1,t)];
    x2 = [x2 a2*x2(1,t)];
    y = [y c1*x1(1,t+1)+c2*x2(1,t+1)];
    varx1 = a1^2*varx1;
    varx2 = a2^2*varx2;
    vary = c1^2*varx1 + c2^2*varx2;
    xhat1 = [xhat1 y(1,t+1)*c1*varx1/vary];
    xhat2 = [xhat2 y(1,t+1)*c2*varx2/vary];
    sigman1 = [sigman1 (1-c1*varx1/vary)^2*varx1+(c2*varx1/vary)^2*varx2];
end

%hold all
subplot(1,2,1), stem(0:t,x1,'k','filled')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, stem(0:t,xhat1,'r','filled')
subplot(1,2,1), hold on, stem(0:t, x2, 'k')
subplot(1,2,1), hold on, stem(0:t, xhat2, 'r')
xlabel('n=time')
legend1 = legend('x1','xhat1','x2','xhat2','Location','Best');

%subplot(1,2,2), stem(0:t,abs(y*c1/(c1+c2)-x1),'b','filled')
subplot(1,2,2), hold on, stem(0:t,xhat1-x1,'r','filled')
%subplot(1,2,2), hold on, stem(0:t, sqrt(sigman1), 'k')
subplot(1,2,2), hold on, stem(0:t,xhat2-x2,'g')
%subplot(1,2,2), hold on, stem(0:t, abs(xhat1-x1+xhat2-x2), 'k', 'filled')
xlabel('n=time')
legend2 = legend('xhat1-x1','xhat2-x2');
subplot(1,2,2), title('Error')

suptitle(['a1= ' num2str(a1) '; a2 = ' num2str(a2)])
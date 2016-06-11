%Kalman Filter Dynamics with Memory

clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a1 = 2; %If a<1, error converges to 0
a2 = 3;
c1 = 1; %If c != 1, simulation fails (xhatmem worse than xhat)
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
xhatmem1 = y*c1*varx1/vary;
xhatmem2 = y*c2*varx2/vary;

sigman1 = varx1 -2*c1*varx1^2/vary + (c1*varx1/vary)^2*vary; %cov(xn-xhatn)
sigman2 = varx2 -2*c2*varx2^2/vary + (c2*varx2/vary)^2*vary;
%sigman1 = (1-c1*varx1/vary)^2*varx1+(c2*varx1/vary)^2*varx2;
sigmat = varx1*varx2/vary*(c1*c2-c1-c2);
sigmap1 = c1*c2*varx2*varx1/vary;
sigmap2 = c1*c2*varx2*varx1/vary;
sigmac1 = c1^2*varx1^2/vary;
sigmac2 = c2^2*varx2^2/vary;

%time passing
for t=1:(n-1)
    x1 = [x1 a1*x1(1,t)];
    x2 = [x2 a2*x2(1,t)];
    y = [y c1*x1(1,t+1)+c2*x2(1,t+1)];
    
    %sigmaprime = c1*c2*varx2*varx1/vary;
    %sigmat = a1*a2*varx1*varx2/vary*(c1*c2-c1-c2);
    
    sn1 = a1^2*sigman1; %generating sn for this iteration
    sn2 = a2^2*sigman2;
    
    knd = c1^2*sn1+c2^2*sn2+2*c1*c2*sigmat;
    ytilde = y(1,t+1)-c1*a1*xhatmem1(1,t)-c2*a2*xhatmem2(1,t);
    
    kn1 = (sn1*c1-a1*a2*c2*sigmap1(1,t))/knd;
    kn2 = (sn2*c2-a2*a1*c1*sigmap2(1,t))/knd;
    xhatmem1 = [xhatmem1 a1*xhatmem1(1,t)+kn1*ytilde];
    xhatmem2 = [xhatmem2 a2*xhatmem2(1,t)+kn2*ytilde];
    
    %sigman1 = (1-kn1*c1)*sn1; %generate new sigma_n for next iteration
    %sigman2 = (1-kn2*c2)*sn2;
    sigman1 = varx1 -2*c1*varx1^2/vary + (c1*varx1/vary)^2*vary; %cov(xn-xhatn)
    sigman2 = varx2 -2*c2*varx2^2/vary + (c2*varx2/vary)^2*vary;
    sigmat = -kn1*c2*(1-kn2*c2)*sn2 - kn2*c1*(1-kn1*c1)*sn1 + a1*a2*(2*kn1*kn2*c2*c1+1-kn1*c1-kn2*c2)*sigmat; %generating new sigmat
    
    sigmap1 = [sigmap1 a1*a2*sigmap1+kn2*(c1*varx1-c1*a1^2*sigmac1(1,t)-c2*a2*a1*sigmap1)];
    sigmap2 = [sigmap2 a1*a2*sigmap2+kn1*(c2*varx2-c2*a2^2*sigmac2(1,t)-c1*a1*a2*sigmap2)];
    sigmac1 = [sigmac1 a1^2*sigmac1+kn1*c1*varx1-kn1*a1^2*sigmac1-kn1*c2*a1*a2*sigmap1(1,t)];
    sigmac2 = [sigmac2 a2^2*sigmac2+kn2*c2*varx2-kn2*a2^2*sigmac2-kn2*c1*a1*a2*sigmap2(1,t)];
    
    varx1 = a1^2*varx1;
    varx2 = a2^2*varx2;
    vary = c1^2*varx1 + c2^2*varx2;
    
    xhat1 = [xhat1 y(1,t+1)*c1*varx1/vary];
    xhat2 = [xhat2 y(1,t+1)*c2*varx2/vary];
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
subplot(1,2,2), hold on, stem(0:t,xhatmem1-x1, 'k', 'filled')
%subplot(1,2,2), hold on, stem(0:t, abs(xhat1-x1+xhat2-x2), 'k', 'filled')
xlabel('n=time')
legend2 = legend('xhat1-x1','xhat2-x2','xhatmem1-x1');
subplot(1,2,2), title('Error')

suptitle(['a1= ' num2str(a1) '; a2 = ' num2str(a2)])

sum(abs(xhat1-x1))/n
sum(abs(xhatmem1-x1))/n

%%
%Noiseless dynamics with memory

clc;
clear all;

a1 = 5;
a2 = 2;
c1 = 1;
c2 = 1;
n = 25;

%initialize
x1 = normrnd(0,1,[1,1]); %x[0]
x2 = normrnd(0,1,[1,1]);
y = c1*x1+c2*x2;

varx1 = 1;

xhat1 = y*c1/(c1^2+c2^2);
xhat2 = y*c2/(c1^2+c2^2);

varw1 = 5;
varw2 = 5;

%time passing
for t=1:(n-1)
    x1 = [x1 a1*x1(1,t)+normrnd(0,varw1)];%added noise
    x2 = [x2 a2*x2(1,t)+normrnd(0,varw2)];
    y = [y c1*x1(1,t+1)+c2*x2(1,t+1)];
    
    xhat2 = [xhat2 a2*(y(1,t+1)-a1*y(1,t))/c2/(a2-a1)];
    %xhat1 = [xhat1 (y(1,t)-((y(1,t+1)-a1*y(1,t))/(a2-a1)))/c1];
    %xhat1 = [xhat1 a1*(y(1,t)-c2*xhat2(1,t))/c1]; %because the first xhat sucks
    xhat1 = [xhat1 a1*(a2*y(1,t)-y(1,t+1))/c1];
    
    varx1 = a1^2*varx1 + varw1;
end

subplot(1,2,1), stem(0:t,x1,'k','filled')
subplot(1,2,1), hold on, stem(0:t, x2, 'k')
subplot(1,2,1), hold on, stem(0:t, xhat1, 'r','filled')
subplot(1,2,1), hold on, stem(0:t, xhat2, 'r')
xlabel('n=time')
legend1 = legend('x1','x2','xhat1','xhat2','Location','Best');

subplot(1,2,2), stem(0:t,abs(x1-xhat1),'r','filled')
subplot(1,2,2), hold on, stem(0:t, abs(x2-xhat2), 'g', 'filled')
xlabel('n=time')
%legend2 = legend('xhat1-x1','xhat2-x2','xhatmem1-x1');
subplot(1,2,2), title('Error')

suptitle(['a1= ' num2str(a1) '; a2 = ' num2str(a2)])

sum(abs(x1-xhat1))/n
sum(abs(x2-xhat2))/n
sum((xhat1-x1).^2)/n
sum((xhat2-x2).^2)/n

a1^2*(1-a2+a1)^2*varx1 + (1+a1)^2*varw1 + (a1*c2/c1)^2*varw2
(1-a2/(a2-a1))^2*varw2 + (a2*c1/c2/(a2-a1))^2*varw1

%%

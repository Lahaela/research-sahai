%Kalman Filter with Additive Noise

clc;
clear all;

a = 1; %Test 2,0.5,10
c = 1;
varv = 1; %channel noise
varw = 1;
n = 500;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);
xhatmem = y*c/(c^2+varv);
varx = 1;
vary = c^2+varv;

sigman = varx -2*c^2*varx^2/vary + c^2*varx^2/vary^2*(c^2*varx+varv);

%sn = a^2+varw - 2*a^2*c^2/(c^2+varv) + a^2*c^4/(c^2+varv)^2 + a^2*c^2/(c^2+varv)^2*varv;


%time passing
for t=1:(n-1)
    w = normrnd(0,varw);
    x = [x a*x(1,t)+w];
    v = normrnd(0,varv);
    y = [y c*x(1,(t+1))+v];
    
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn*c/(c^2*sn+varv);
    xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t+1)-c*a*xhatmem(1,t))];
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn*c)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
end

xtilde = x-xhat;
xmemtilde = x-xhatmem;

subplot(3,2,1), hold on, plot(0:t,x,'k')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(3,2,1), hold on, plot(0:t,xhat,'r')
subplot(3,2,1), hold on, plot(0:t,xhatmem, 'g')
xlabel('n=time')
legend1 = legend('x','xhat','xhatmem','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(3,2,2), hold on, plot(0:t,xtilde,'r')
subplot(3,2,2), hold on, plot(0:t,xmemtilde,'b')
xlabel('n=time')
legend2 = legend('x-xhat','x-xhatmem','Location','Best');
subplot(3,2,2), title('Error')

subplot(3,2,3), hist(xmemtilde);
subplot(3,2,3), title('Kalman Filter Estimation Error')

xmemtilde_vec = linspace(min(xmemtilde),max(xmemtilde));
data = zeros(length(xmemtilde_vec),1);
for i = 1:length(xmemtilde_vec)
    data(i,1) = sum(xmemtilde < xmemtilde_vec(i));
end
data = data./length(xmemtilde); 

subplot(3,2,4), plot(xmemtilde_vec,data,'LineWidth',2)
subplot(3,2,4), xlabel('m = magnitude of error'), ylabel('P(Xmemtilde < m)')
subplot(3,2,4), title('CDF of KF Estimation Error')

xmemtilde = (x-xhatmem).^2; %This is to make the bottom plots error squared

subplot(3,2,5), hist(xmemtilde);
subplot(3,2,5), title('Kalman Filter Estimation Error Squared')

xmemtilde_vec = linspace(min(xmemtilde),max(xmemtilde));
data = zeros(length(xmemtilde_vec),1);
for i = 1:length(xmemtilde_vec)
    data(i,1) = sum(xmemtilde < xmemtilde_vec(i));
end
data = data./length(xmemtilde); 

subplot(3,2,6), plot(xmemtilde_vec,data,'LineWidth',2)
subplot(3,2,6), xlabel('m = magnitude of error squared'), ylabel('P((X-Xhat)^2 < m)')
subplot(3,2,6), title('CDF of KF Estimation Error Squared')

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig1.pdf');

sum(abs(xhat-x))/n
sum(abs(xhatmem-x))/n

%%
%Kalman Filter Ignoring Multiplicative Noise

clc;
clear all;

a = 1; %Test 2,0.5,10
c = 1;
varv = 1; %multiplicative noise, 1/p
varw = 1;%system noise
n = 1000;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = (c+normrnd(0,varv))*x;
xhat = y*c*1/c^2;
xhatmem = y*c/c^2;
varx = 1;
vary = c^2;

sigman = varx -2*c^2*varx^2/vary + c^2*varx^2/vary^2*(c^2*varx);

%sn = a^2+varw - 2*a^2*c^2/(c^2+varv) + a^2*c^4/(c^2+varv)^2 + a^2*c^2/(c^2+varv)^2*varv;


%time passing
for t=1:(n-1)
    w = normrnd(0,varw);
    x = [x a*x(1,t)+w];
    v = normrnd(0,varv);
    y = [y (c+v)*x(1,(t+1))];
    
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn*c/(c^2*sn);
    xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t+1)-c*a*xhatmem(1,t))];
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn*c)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    varx = a^2*varx + varw;
    vary = c^2*varx;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
end

xtilde = x-xhat;
xmemtilde = x-xhatmem;

subplot(3,2,1), hold on, plot(0:t,x,'k')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(3,2,1), hold on, plot(0:t,xhat,'r')
subplot(3,2,1), hold on, plot(0:t,xhatmem, 'g')
xlabel('n=time')
legend1 = legend('x','xhat','xhatmem','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(3,2,2), hold on, plot(0:t,xtilde,'r')
subplot(3,2,2), hold on, plot(0:t,xmemtilde,'b')
xlabel('n=time')
legend2 = legend('x-xhat','x-xhatmem','Location','Best');
subplot(3,2,2), title('Error')

subplot(3,2,3), hist(xmemtilde);
subplot(3,2,3), title('Kalman Filter Estimation Error')

xmemtilde_vec = linspace(min(xmemtilde),max(xmemtilde));
data = zeros(length(xmemtilde_vec),1);
for i = 1:length(xmemtilde_vec)
    data(i,1) = sum(xmemtilde < xmemtilde_vec(i));
end
data = data./length(xmemtilde); 

subplot(3,2,4), plot(xmemtilde_vec,data,'LineWidth',2)
subplot(3,2,4), xlabel('m = magnitude of error'), ylabel('P(Xmemtilde < m)')
subplot(3,2,4), title('CDF of KF Estimation Error')

xmemtilde = (x-xhatmem).^2; %This is to make the bottom plots error squared

subplot(3,2,5), hist(xmemtilde);
subplot(3,2,5), title('Kalman Filter Estimation Error Squared')

xmemtilde_vec = linspace(min(xmemtilde),max(xmemtilde));
data = zeros(length(xmemtilde_vec),1);
for i = 1:length(xmemtilde_vec)
    data(i,1) = sum(xmemtilde < xmemtilde_vec(i));
end
data = data./length(xmemtilde); 

subplot(3,2,6), plot(xmemtilde_vec,data,'LineWidth',2)
subplot(3,2,6), xlabel('m = magnitude of error squared'), ylabel('P((X-Xhat)^2 < m)')
subplot(3,2,6), title('CDF of KF Estimation Error Squared')


suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','fig2_v0.01.pdf');

sum(abs(xhat-x))/n
sum(abs(xhatmem-x))/n

%%
%Kalman Filter Additive+Multiplicative Noise

clc;
clear all;

a = 1; %must be <1
%c = 1;
varv = 1; %additive noise
varw = 1; %system noise
p = 0.01; %multiplicative noise, test 0.01, 1, 10, 100
n = 100;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = x+normrnd(0,varv);

xhatt = y/(1+varv);
xhatr = [0]; %no initial xhatr
st = xhatt(1,1)-xhatr(1,1); % no initial xhatr so st = xhatt
z = (1+normrnd(0,1/p))*st;
xhatf = xhatr(1,1)+p/(p+1)*z(1,1);

varx = 1;
vary = 1+varv;

sigman = varx -2*varx^2/vary + varx^2/vary^2*(varx+varv);

%time passing
for t=1:(n-1)
    x = [x a*x(1,t)+normrnd(0,varw)];
    y = [y x(1,(t+1))+normrnd(0,varv)];
    
    %transmitter side
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn/(sn+varv);
    xhatt = [xhatt a*xhatt(1,t)+kn*(y(1,t+1)-a*xhatt(1,t))]; %t|t
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    xhatr = [xhatr a*(xhatr(1,t)+p/(p+1)*(z(1,t)))]; %t|t-1, -zhat(1,t)=0 because gamma =1
    
    st = [st xhatt(1,t+1)-xhatr(1,t+1)];
    
    z = [z (1+normrnd(0,1/p))*st(1,t+1)];
    xhatf = [xhatf xhatr(1,t+1)+p/(p+1)*z(1,t+1)];
    
    varx = a^2*varx + varw;
    vary = varx + varv;
    
end

xttilde = x-xhatt;
xrtilde = x-xhatr;
xftilde = x-xhatf;

subplot(3,2,1), hold on, plot(0:t,x,'k')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(3,2,1), hold on, plot(0:t,xhatt, 'g')
subplot(3,2,1), hold on, plot(0:t,xhatf,'r')
xlabel('n=time')
legend1 = legend('x','xhatt','xhatf','Location','Best');
subplot(3,2,1), title('Signal')

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(3,2,2), hold on, plot(0:t,xttilde,'b')
subplot(3,2,2), hold on, plot(0:t,xftilde,'r')
xlabel('n=time')
legend2 = legend('x-xhatt','x-xhatf','Location','Best');
subplot(3,2,2), title('Error')

subplot(3,2,3), hist(xftilde,100);
subplot(3,2,3), ylabel(['n = ' num2str(n)]), title('Kalman Filter Final Estimation Error')

xttilde_vec = linspace(min(xttilde),max(xttilde));
data = zeros(length(xttilde_vec),1);
for i = 1:length(xttilde_vec)
    data(i,1) = sum(xttilde < xttilde_vec(i));
end
data = data./length(xttilde); 

xftilde_vec = linspace(min(xftilde),max(xftilde));
dataf = zeros(length(xftilde_vec),1);
for j = 1:length(xftilde_vec)
    dataf(j,1) = sum(xftilde < xftilde_vec(j));
end
dataf = dataf./length(xftilde);

subplot(3,2,4), plot(xttilde_vec,data,'Color','b','LineWidth',2)
subplot(3,2,4), hold on, plot(xftilde_vec,dataf,'Color','r','LineWidth',2)
subplot(3,2,4), xlabel('m = magnitude of error'), ylabel('P(X-tilde < m)')
subplot(3,2,4), title('CDF of KF Estimation Error')
legend3 = legend('xttilde','xftilde','Location','Best');

%Bottom plots are error squared
xttilde = (x-xhatt).^2;
xftilde = (x-xhatf).^2;

subplot(3,2,5), hist(xftilde,100);
subplot(3,2,5), ylabel(['n = ' num2str(n)]), title('Kalman Filter Final Estimation Error Squared')

xttilde_vec = linspace(min(xttilde),max(xttilde));
data = zeros(length(xttilde_vec),1);
for i = 1:length(xttilde_vec)
    data(i,1) = sum(xttilde < xttilde_vec(i));
end
data = data./length(xttilde); 

xftilde_vec = linspace(min(xftilde),max(xftilde));
dataf = zeros(length(xftilde_vec),1);
for j = 1:length(xftilde_vec)
    dataf(j,1) = sum(xftilde < xftilde_vec(j));
end
dataf = dataf./length(xftilde);

subplot(3,2,6), plot(xttilde_vec,data,'Color','b','LineWidth',2)
subplot(3,2,6), hold on, plot(xftilde_vec,dataf,'Color','r','LineWidth',2)
subplot(3,2,6), xlabel('m = magnitude of error squared'), ylabel('P(X-tilde^2 < m)')
subplot(3,2,6), title('CDF of KF Estimation Error Squared')
legend4 = legend('xttilde','xftilde','Location','Best');

suptitle(['a = ' num2str(a) '; varv = ' num2str(varv) '; varw = ' num2str(varw) '; p = ' num2str(p)])

% legend3 = legend('xttilde','xftilde','Location','Best');
set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig3_p0.01_n100.pdf');

sum(abs(xttilde))/n
sum(abs(xrtilde))/n
sum(abs(xftilde))/n
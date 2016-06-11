clc;
clear all;
set(0,'DefaultAxesFontSize', 14);

a = 1; %Test 2,0.5,10
c = 1;
varv = 1;
varw = 1;
n = 40;

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
    
    knum = a^2*c*varx-c^3*a^2*varx^2/vary+c*varw;
    kden = c^2*varw+varv+(c^2*a*varx/vary)^2*varv+c^2*a^2*(1-c^2*varx/vary)^2*varx;
    
    xhatmem = [xhatmem a*c*varx/vary*y(1,t)+knum/kden*(y(1,t+1)-c^2*a*varx/vary*y(1,t))];
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
    t
end

%hold all
subplot(1,2,1), stem(0:t,x,'k')
subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, stem(0:t,xhat,'r')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

%subplot(1,2,2), stem(0:t,abs(y/c-x),'b','filled')
subplot(1,2,2), hold on, stem(0:t,abs(xhat-x),'r','filled')
subplot(1,2,2), hold on, stem(0:t,abs(xhatmem-x),'g','filled')
xlabel('n=time')
legend2 = legend('y/c-x','xhat-x','xhatmem-x');

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

%%
%Kalman Filter (Prediction of Present)
clc
clear all

a = 1; %Test 2,0.5,10
c = 1;
varv = 0.09; %channel noise
varw = 0.04;
n = 10;

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

subplot(1,2,1), stem(0:t,x,'k','filled')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, stem(0:t,xhat,'r')
subplot(1,2,1), hold on, stem(0:t, xhatmem, 'g','filled')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(1,2,2), hold on, plot(0:t,abs(xhat-x),'r')
subplot(1,2,2), hold on, plot(0:t,abs(xhatmem-x),'b')
xlabel('n=time')
legend2 = legend('xhat-x','xhatmem-x');
subplot(1,2,2), title('Error')

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

%print -dpdf output2.pdf

sum(abs(xhat-x))/n
sum(abs(xhatmem-x))/n

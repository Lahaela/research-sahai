clc
clear all

%Noiseless Dynamics Two Observation Vector

A = [2 0; 0 3];
C = [1 2];
n = 25;

X = [normrnd(0,1); normrnd(0,1)];
Y = C*X;

Xhat = [Y(1,1)*C(1,1)/(C(1,1)^2+C(1,2)^2); Y(1,1)*C(1,2)/(C(1,1)^2+C(1,2)^2)];

varw1 = 1;
varw2 = 1;
varW = [1; 1];

for t=1:(n-1)
    W = [normrnd(0,varw1); normrnd(0,varw2)];
    
    X = [X A*X(1:2,t)+W];
    Y = [Y C*X(1:2,t+1)];
    
    Xhat = [Xhat A*([C; C*A]^(-1)*[Y(1,t); Y(1,t+1)])];
end

stem(0:t,abs(X(1,:)-Xhat(1,:)),'r','filled')
hold on, stem(0:t, abs(X(2,:)-Xhat(2,:)), 'g')
xlabel('n=time')
legend2 = legend('xhat1-x1','xhat2-x2');
title('Error')

%suptitle(['a1= ' num2str(a1) '; a2 = ' num2str(a2)])

sum(abs(X(1,:)-Xhat(1,:)))/n
sum(abs(X(2,:)-Xhat(2,:)))/n
sum((X(1,:)-Xhat(1,:)).^2)/n
sum((X(2,:)-Xhat(2,:)).^2)/n

(1+A(1,1)/(A(2,2)-A(1,1)))^2*varw1 + (A(1,1)*C(1,2)/C(1,1)/(A(2,2)-A(1,1)))^2*varw2
(1-A(2,2)/(A(2,2)-A(1,1)))^2*varw2 + (A(2,2)*C(1,1)/C(1,2)/(A(2,2)-A(1,1)))^2*varw1

%%

clc
clear all

%Noiseless Dynamics Two Observation Vector with Y(0), Y(n)

A = [2 0; 0 3];
C = [1 1];
n = 25;

X = [normrnd(0,1); normrnd(0,2)];
Y = C*X;

Xhat = [Y(1,1)*C(1,1)/(C(1,1)^2+C(1,2)^2); Y(1,1)*C(1,2)/(C(1,1)^2+C(1,2)^2)];

varw1 = 1;
varw2 = 1;
varW = [1; 1];

for t=1:(n-1)
    W = [normrnd(0,varw1); normrnd(0,varw2)];
    
    X = [X A*X(1:2,t)+W];
    Y = [Y C*X(1:2,t+1)];
    
    Xhat = [Xhat A^(t+1)*([C; C*A^(t+1)]^(-1)*[Y(1,1); Y(1,t+1)])];
end

hold on
stem(0:t,abs(X(1,:)-Xhat(1,:)),'r','filled')
stem(0:t, abs(X(2,:)-Xhat(2,:)), 'g', 'filled')
xlabel('n=time')
%legend2 = legend('xhat1-x1','xhat2-x2','xhatmem1-x1');
title('Error')

%suptitle(['a1= ' num2str(a1) '; a2 = ' num2str(a2)])

sum(abs(X(1,:)-Xhat(1,:)))/n
sum(abs(X(2,:)-Xhat(2,:)))/n
sum((X(1,:)-Xhat(1,:)).^2)/n
sum((X(2,:)-Xhat(2,:)).^2)/n

%%
%Control Attempt 1

clc
clear all

a = 1; %Test 2,0.5,10
c = 2;
varv = 0.000001; %channel noise
varw = 0;
n = 100;

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
    x = [x a*(x(1,t)-xhatmem(1,t))];
    v = normrnd(0,varv);
    y = [y c*x(1,(t+1))+v];
    
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn*c/(c^2*sn+varv);
    xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t+1)-c*a*xhatmem(1,t))];
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn*c)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    varx = a^2*(1-c^2*varx/vary)^2*varx + varw + (a*c*varx/vary)^2*varv;
    %varx = a^2*varx + varw + a^2*c^2*varx^2/vary; %add in var(u)?
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
end

subplot(1,2,1), plot(0:t,x,'k')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, plot(0:t,xhat,'r')
subplot(1,2,1), hold on, plot(0:t, xhatmem, 'g')
xlabel('n=time')
legend1 = legend('x','xhat','xhatmem','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(1,2,2), hold on, plot(0:t,abs(xhat-x),'k')
subplot(1,2,2), hold on, plot(0:t,abs(xhatmem-x),'b')
%subplot(1,2,2), hold on, plot(0:t, 1./sqrt(1:500),'r')
xlabel('n=time')
legend2 = legend('xhat-x','xhatmem-x');
subplot(1,2,2), title('Error')

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

sum(abs(xhat-x))/n
sum(abs(xhatmem-x))/n
sum(abs(x))/n

%'Error'
%a*((1-c^2*varx/vary)^2*varx + (c*varx/vary)^2*varv)
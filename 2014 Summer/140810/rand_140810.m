%Kalman Filter Schenato Try 2
clc;
clear all;

n = 200;
M = 10000;

a = 0.99; %must be <1
%c = 1;
varv = 1; %additive noise
varw = 1; %system noise
% p = 100; %multiplicative noise, test 0.01, 1, 10, 100
varx = 1;
vary = 1+varv;

P = [100 10 1 0.01];

for p=P
Xttilde= [];
Xrtilde = [];
Xftilde = [];

for m=1:M
%initialize
x = normrnd(0,sqrt(varx),[1,1]); %x[0]
y = x+normrnd(0,sqrt(varv));

xhatt = y/(1+varv);
xhatr = [0]; %no initial xhatr
st = xhatt(1,1)-xhatr(1,1); % no initial xhatr so st = xhatt
z = (1+normrnd(0,sqrt(1/p)))*st;
xhatf = xhatr(1,1)+p/(p+1)*z(1,1);

sigman = varx -2*varx^2/vary + varx^2/vary^2*(varx+varv);

%time passing
for t=1:(n-1)
    x = [x a*x(1,t)+normrnd(0,sqrt(varw))];
    y = [y x(1,(t+1))+normrnd(0,sqrt(varv))];
    
    %transmitter side
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn/(sn+varv);
    xhatt = [xhatt a*xhatt(1,t)+kn*(y(1,t+1)-a*xhatt(1,t))]; %t|t
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    xhatr = [xhatr a*xhatf(1,t)]; %t|t-1, -zhat(1,t)=0 because gamma =1 (xhatr(1,t)+p/(p+1)*(z(1,t)))
    
    st = [st xhatt(1,t+1)-xhatr(1,t+1)];
    
    z = [z (1+normrnd(0,sqrt(1/p)))*st(1,t+1)];
    xhatf = [xhatf xhatr(1,t+1)+p/(p+1)*z(1,t+1)];
    
end

Xttilde = [Xttilde (x-xhatt)'];
Xrtilde = [Xrtilde (x-xhatr)'];
Xftilde = [Xftilde (x-xhatf)'];
end

i=find(p==P);
subplot(2,2,i), plot(0:t,var(Xttilde'),'b')
subplot(2,2,i), hold on, plot(0:t, var(Xftilde'),'r')
subplot(2,2,i), hold on, plot(0:t, var(Xrtilde'),'g')
subplot(2,2,i), title(['P = ' num2str(p)]), xlabel('n = time')

end

subplot(2,2,1), legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')
subplot(2,2,4), legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')


%Squared Error
% Xttilde = Xttilde.^2;
% Xrtilde = Xrtilde.^2;
% Xftilde = Xftilde.^2;

suptitle(['Variance of Error; A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','schenato_varyP.pdf');
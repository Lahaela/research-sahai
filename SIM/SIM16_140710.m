% Rajasekaran vs. Schenato: R

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 0.5; %Test 2,0.5,10
b = 1;
%c = 1;
varu = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_u = 2;
varv = 0.001; %additive observation noise
varw = 0.001;%system noise

n = 500;
M = 100;
Xmtildep = []; %prediction
Xmtildee = []; %estimate

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,varx,[1,1]); %x[0]
    %x = [10];
    u = normrnd(mu_u,varu);
    y = u*x+normrnd(0,varv);

    xhatm = y*mu_u*varx/vary;
    xhatp = [0]; %no initial prediction

    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,varw)];
        u = normrnd(mu_u,varu);
        y = [y u*x(1,(t+1))+normrnd(0,varv)];
        
        varx = a^2*varx+b^2*varw;
        
        xhatp = [xhatp a*xhatm(1,t)];
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u*sn_+varu*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n

    end
    
    Xmtildep = [Xmtildep; x-xhatp];
    Xmtildee = [Xmtildee; x-xhatm];
end

subplot(4,2,1), plot(0:t,mean(abs(Xmtildep)), 'r')
subplot(4,2,1), hold on, plot(0:t,mean(abs(Xmtildee)),'b')
subplot(4,2,1), title('Average Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X_p}$$','$$\tilde{X_e}$$','Location','Best');

subplot(4,2,2), plot(0:t,var(Xmtildep),'r')
subplot(4,2,2), hold on, plot(0:t,var(Xmtildee),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X_p}$$','$$\tilde{X_e}$$','Location','Best');

Xmtildep = Xmtildep.^2;
Xmtildee = Xmtildee.^2;

subplot(4,2,3), hist(Xmtildee(:,1),20)
subplot(4,2,3), title('Squared Error of $$\hat{X}_e(0)$$')

[vec, data] = cdfld(Xmtildee(:,1)); [vec2, data2] = cdfld(Xmtildep(:,1));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')

subplot(4,2,5), hist(Xmtildee(:,n/2),20)
subplot(4,2,5), title(['Squared Error of $$\hat{X}_e$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xmtildee(:,n/2)); [vec2, data2] = cdfld(Xmtildep(:,n/2));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
legend('$$\tilde{X_e}$$','$$\tilde{X_p}$$','Location','Best');

subplot(4,2,7), hist(Xmtildee(:,n),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}_e$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtildee(:,n)); [vec2, data2] = cdfld(Xmtildep(:,n));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','raj_a05_p001.pdf');

%%
%Rajasekaran vs Schenato: S

clc;
clear all;

n = 500;
M = 1000;

a = 0.5; %must be <1
%c = 1;
varv = 1; %additive noise
varw = 1; %system noise
p = 1; %multiplicative noise, test 0.01, 1, 10, 100
varx = 1;
vary = 1+varv;

Xttilde= [];
Xrtilde = [];
Xftilde = [];

for m=1:M
%initialize
x = normrnd(0,varx,[1,1]); %x[0]
y = x+normrnd(0,varv);

xhatt = y/(1+varv);
xhatr = [0]; %no initial xhatr
st = xhatt(1,1)-xhatr(1,1); % no initial xhatr so st = xhatt
z = (1+normrnd(0,1/p))*st;
xhatf = xhatr(1,1)+p/(p+1)*z(1,1);

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
    
    xhatr = [xhatr a*xhatf(1,t)]; %t|t-1, -zhat(1,t)=0 because gamma =1 (xhatr(1,t)+p/(p+1)*(z(1,t)))
    
    st = [st xhatt(1,t+1)-xhatr(1,t+1)];
    
    z = [z (1+normrnd(0,1/p))*st(1,t+1)];
    xhatf = [xhatf xhatr(1,t+1)+p/(p+1)*z(1,t+1)];
    
end

Xttilde = [Xttilde; (x-xhatt)];
Xrtilde = [Xrtilde; (x-xhatr)];
Xftilde = [Xftilde; (x-xhatf)];
end

subplot(4,2,1), plot(0:t,mean(abs(Xftilde)),'b',0:t,mean(abs(Xrtilde)),'r')
subplot(4,2,1), title('Average Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

subplot(4,2,2), plot(0:t, var(Xftilde),'b',0:t, var(Xrtilde),'r')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

%Squared Error
Xrtilde = Xrtilde.^2;
Xftilde = Xftilde.^2;

subplot(4,2,3), hist(Xftilde(:,1),50)
subplot(4,2,3), title('Squared Error of $$\hat{Xf(0)}$$')

[vec2, data2] = cdfld(Xftilde(:,1)); [vec3, data3] = cdfld(Xrtilde(:,1));
subplot(4,2,4), plot(vec2, data2, 'b', 'LineWidth',2), hold on, plot(vec3, data3, 'r', 'LineWidth',2)
legend({'$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')
subplot(4,2,4), title('CDF of $$\hat{X}(0)$$ Squared Error'), xlabel('e = Magnitude'), ylabel('$$P((X(0)-\hat{X}(0))^2 < e)$$')

subplot(4,2,5), hist(Xftilde(:,n/2),50)
subplot(4,2,5), title(['Squared Error of $$\hat{X}f$$(' num2str(n/2) ')'])

[vec2, data2] = cdfld(Xftilde(:,n/2)); [vec3, data3] = cdfld(Xrtilde(:,n/2));
subplot(4,2,6), plot(vec2, data2, 'b', 'LineWidth',2), hold on, plot(vec3, data3, 'r', 'LineWidth',2)
subplot(4,2,6), title(['CDF of $$\hat{X}$$(' num2str(n/2) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n/2) '$$)-\hat{X}($$' num2str(n/2) '$$))^2 <$$ e)'])

subplot(4,2,7), hist(Xftilde(:,n),50)
subplot(4,2,7), title(['Squared Error of $$\hat{X}f$$(' num2str(n) ')'])

[vec2, data2] = cdfld(Xftilde(:,n)); [vec3, data3] = cdfld(Xrtilde(:,n));
subplot(4,2,8), plot(vec2, data2, 'b', 'LineWidth',2), hold on, plot(vec3, data3, 'r', 'LineWidth',2)
subplot(4,2,8), title(['CDF of $$\hat{X}$$(' num2str(n) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n) ')-$$\hat{X}$$(' num2str(n) '$$))^2 <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; P = ' num2str(p) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','schenato_a05_p1_m1000.pdf');

%%
% Rajasekaran vs. Schenato: R

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all

a = 0.5; %Test 2,0.5,10
b = 1;
%c = 1;
varu = 100; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_u = 1;
varv = 1; %additive observation noise
varw = 1;%system noise

n = 50;
M = 1;
Xmtildep = []; %prediction
Xmtildee = []; %estimate
Xmp = [];
Xme = [];
X = [];

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,varx,[1,1]); %x[0]
    u = normrnd(mu_u,varu);
    y = u*x+normrnd(0,varv);

    xhatm = y*mu_u*varx/vary;
    xhatp = [0]; %no initial prediction

    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,varw)];
        u = normrnd(mu_u,varu);
        y = [y u*x(1,(t+1))+normrnd(0,varv)];
        
        varx = a^2*varx+b^2*varw;
        
        xhatp = [xhatp a*xhatm(1,t)];
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u*sn_+varu*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n

    end
    
%     Xmp = [Xmp; xhatp];
%     Xme = [Xme; xhatm];
%     X = [X; x];
    
    Xmtildep = [Xmtildep; x-xhatp];
    Xmtildee = [Xmtildee; x-xhatm];
end

%subplot(4,2,1), plot(0:t,mean(abs(Xmtildep)), 'r')
%subplot(4,2,1), hold on, plot(0:t,mean(abs(Xmtildee)),'b')
plot(0:t, x,'k')
hold on, plot(0:t, xhatp,'r')
hold on, plot(0:t,xhatm,'b')
title('X'), xlabel('n = time'), ylabel('Magnitude')
legend('X','$$X_p$$','$$X_e$$','Location','Best');

% subplot(2,2,2), plot(0:t,var(x'),'k')
% subplot(2,2,2), hold on, plot(0:t,var(xhatp'),'r',0:t,var(xhatm'),'b')
% subplot(2,2,2), title('Variance of X'), xlabel('n = time'), ylabel('Magnitude')
% legend('X','$$X_p$$','$$X_e$$','Location','Best');

% Xmtildep = Xmtildep.^2;
% Xmtildee = Xmtildee.^2;
% 
% subplot(4,2,3), hist(Xmtildee(:,1),20)
% subplot(4,2,3), title('Squared Error of $$\hat{X}_e(0)$$')
% 
% [vec, data] = cdfld(Xmtildee(:,1)); [vec2, data2] = cdfld(Xmtildep(:,1));
% subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
% subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')
% 
% subplot(4,2,5), hist(Xmtildee(:,n/2),20)
% subplot(4,2,5), title(['Squared Error of $$\hat{X}_e$$(' num2str(n/2) ')'])
% 
% [vec, data] = cdfld(Xmtildee(:,n/2)); [vec2, data2] = cdfld(Xmtildep(:,n/2));
% subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
% legend('$$\tilde{X_e}$$','$$\tilde{X_p}$$','Location','Best');
% 
% subplot(4,2,7), hist(Xmtildee(:,n),20)
% subplot(4,2,7), title(['Squared Error of $$\hat{X}_e$$(' num2str(n) ')'])
% 
% [vec, data] = cdfld(Xmtildee(:,n)); [vec2, data2] = cdfld(Xmtildep(:,n));
% subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','output0.pdf');

 mean(abs(x-xhatm))
  mean(abs(x-xhatp))
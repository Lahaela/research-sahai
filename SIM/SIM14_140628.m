% Kalman Filter Schenato Error Var Curve
% *** UNFINISHED***

clc;
clear all;

n = 100;
M = 100;

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

Xttilde = [Xttilde (x-xhatt)'];
Xrtilde = [Xrtilde (x-xhatr)'];
Xftilde = [Xftilde (x-xhatf)'];
end

subplot(4,2,1), plot(0:t,mean(Xttilde'),'b')
subplot(4,2,1), hold on, plot(0:t,mean(Xftilde'),'r'), plot(0:t,mean(Xrtilde'),'g')
subplot(4,2,1), title('Average Error'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

subplot(4,2,2), plot(0:t,var(Xttilde'),'b')
subplot(4,2,2), hold on, plot(0:t, sigman,'k','LineWidth',2)
subplot(4,2,2), hold on, plot(0:t, var(Xftilde'),'r')
subplot(4,2,2), hold on, plot(0:t, var(Xrtilde'),'g')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

%Squared Error
Xttilde = Xttilde.^2;
Xrtilde = Xrtilde.^2;
Xftilde = Xftilde.^2;

subplot(4,2,3), hist(Xftilde(1,:),M/5)
subplot(4,2,3), title('Squared Error of $$\hat{Xf(0)}$$')

[vec, data] = cdfld(Xttilde(1,:)); [vec2, data2] = cdfld(Xftilde(1,:)); [vec3, data3] = cdfld(Xrtilde(1,:));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')
subplot(4,2,4), title('CDF of $$\hat{X}(0)$$ Squared Error'), xlabel('e = Magnitude'), ylabel('$$P((X(0)-\hat{X}(0))^2 < e)$$')

subplot(4,2,5), hist(Xftilde(n/2,:),M/5)
subplot(4,2,5), title(['Squared Error of $$\hat{X}f$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xttilde(n/2,:)); [vec2, data2] = cdfld(Xftilde(n/2,:)); [vec3, data3] = cdfld(Xrtilde(n/2,:));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
subplot(4,2,6), title(['CDF of $$\hat{X}$$(' num2str(n/2) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n/2) '$$)-\hat{X}($$' num2str(n/2) '$$))^2 <$$ e)'])

subplot(4,2,7), hist(Xftilde(n,:),M/5)
subplot(4,2,7), title(['Squared Error of $$\hat{X}f$$(' num2str(n) ')'])

[vec, data] = cdfld(Xttilde(n,:)); [vec2, data2] = cdfld(Xftilde(n,:)); [vec3, data3] = cdfld(Xrtilde(n,:));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
subplot(4,2,8), title(['CDF of $$\hat{X}$$(' num2str(n) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n) ')-$$\hat{X}$$(' num2str(n) '$$))^2 <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; P = ' num2str(p) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','fig3_p100_m1000.pdf');

%%
%Gireeja NonCoherence Empirical Variance

clc
clear all

set(0,'defaulttextinterpreter','latex')

A = [0.1 0.5 0.9 1.2];
colors = {'r', 'b', 'g', 'k'};
%c = 1;
varx = 1;
varv = 1; %channel noise
varw = 1;
n = 1000;
mu_c = 1;
varc = 1;
M = 1000;

legentries = cell(1,length(A));

for i=1:length(A)
%     Xtilde = [];
%     Xtildet = [];
    X = [];
%     Xhat = [];
    for m=1:M
    %initialize
    c = normrnd(mu_c,sqrt(varc));
    x = normrnd(0,sqrt(varx),[1,1]); %x[0]
    y = c*x+normrnd(0,sqrt(varv));

    vary = varc*varx+varv;
    
%     xhat = y*mu_c*varx/vary;
%     xtildet = (1-mu_c*c*varx/vary)^2*varx+(mu_c*varx/vary)^2*varv;

    %time passing
    for t=1:(n-1)
        c = normrnd(mu_c,sqrt(varc));
        x = [x A(i)*(x(1,t)-mu_c/(mu_c^2+varc)*y(1,t))+normrnd(0,sqrt(varw))]; %-xhat(1,t)
        y = [y c*x(1,(t+1))+normrnd(0,sqrt(varv))];

        %varx = A(i)^2*varx + varw;
        
%         xhat = [xhat y(1,t+1)*mu_c*varx/vary];
%         xtildet = [xtildet (1-mu_c*c*varx/vary)^2*varx+(mu_c*varx/vary)^2*varv];
    end
%     Xtilde = [Xtilde (x-xhat)'];
%     Xtildet = [Xtildet; xtildet];
    X = [X; x];
%     Xhat = [Xhat; xhat];
    
    end
    
    Varxt = [1]; %outside of m loop because theoretical, update varx if necessary
    for t = 1:(n-1)
        varx = A(i)^2*(1-mu_c^2/(mu_c^2+varc))*varx + varw + (A(i)*mu_c/(mu_c^2+varc))^2*varv;
        vary = varc*varx+varv;
        Varxt = [Varxt varx];
    end
    
    legentries{1,i} = ['A=' num2str(A(1,i))];
    
    subplot(2,2,1), hold on, plot(0:t,mean(X),'Color',colors{i});
    subplot(2,2,1), title('Mean X')
    
    subplot(2,2,2), hold on, plot(0:t,var(X),'Color',colors{i});
    subplot(2,2,2), hold on, plot(0:t, Varxt,'m','LineWidth',1.5);
    subplot(2,2,2), title('Variance of X')
    
    [vec, data] = cdfld(X(n,:)); 
    subplot(2,2,3), hold on, plot(vec,data,'Color',colors{i},'LineWidth',2);
    subplot(2,2,3), title('CDF of X at time n')
    
%     Xtilde = Xtilde.^2;
    X2 = X.^2;
    [vec, data] = cdfld(X2(n,:));
    subplot(2,2,4), hold on, plot(vec,data,'Color',colors{i}, 'LineWidth',2);
    subplot(2,2,4), title('CDF of $$X^2$$ at time n')
end

subplot(2,2,1), legend(legentries,'Location','Best');
%subplot(2,2,2), legend(legentries,'Location','Best');
subplot(2,2,4), legend(legentries,'Location','Best');
subplot(2,2,3), legend(legentries,'Location','Best');

suptitle(['Gireeja NonCoherence Paper Varying A M = ' num2str(M)])

%test = mean(Xtildet);

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','output0_2.pdf');

%%
%Gireeja NonCoherence Varying A Curve

clc
clear all

set(0,'defaulttextinterpreter','latex')

A = 0:0.1:1.1;
colors = {'r', 'b', 'g', 'k', 'm'};
%c = 1;
b = 1;
varv = 1; %channel noise
varw = 1;
n = 200;
mu_c = 1;
varc = 1;

Z = [];

for i=1:length(A)

    %initialize
    varx = 1;
    vary = (mu_c^2+varc)*varx+varv;
    
    alphac = mu_c*varx/((mu_c^2+varc^2)*varx+varv);
    
    xtildet = [((1-2*alphac*mu_c+alphac^2*(mu_c^2+varc))*varx+alphac^2*varv)];
    
    %time passing
    for t=1:(n-1)

        sn_ = A(i)^2*xtildet(t)+b^2*varw;
        kf_ = mu_c*sn_/(mu_c^2*sn_+varc*varx+varv);
        xtildet = [xtildet (1-kf_*mu_c)*sn_]; %generating new sigma_n
        
        varx = A(i)^2*varx+b^2*varw;
%         varx = A(i)^2*(1-mu_c^2/(mu_c^2+varc))*varx + varw + (A(i)*mu_c/(mu_c^2+varc))^2*varv;
        vary = (mu_c^2+varc)*varx+varv;
    end

Z = [Z xtildet(n)];

end

plot(A, Z,'m','LineWidth',2); 
title(['Gireeja NonCoherence Paper at N = ' num2str(n)])
xlabel('A'), ylabel('Theoretical Variance of Estimation Error')
% legend('Theoretical','Empirical','Location','Best')

%suptitle(['Gireeja NonCoherence Paper A to Estimation Error Curve M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
print('-dpdf','-r100','output1.pdf');
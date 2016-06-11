%Gireeja NonCoherence Control

clc
clear all

A = [0.1 0.5 0.75];
colors = {'r', 'b', 'g', 'k', 'm'};
%c = 1;
varx = 1;
varv = 1; %channel noise
varw = 1;
n = 100;
mu_c = 1;
varc = 0.1;

for i=1:length(A)
    %initialize
    c = normrnd(mu_c,varc);
    x = normrnd(0,varx,[1,1]); %x[0]
    y = c*x+normrnd(0,varv);

    vary = varc*varx+varv;
    
    xhat = y*mu_c*varx/vary;

    %time passing
    for t=1:(n-1)
        c = normrnd(mu_c,varc);
        x = [x A(i)*(x(1,t))+normrnd(0,varw)]; %-xhat(1,t)
        y = [y c*x(1,(t+1))+normrnd(0,varv)];

        varx = A(i)^2*varx + varw;
        %varx = A(i)^2*(1-mu_c*c*varx/vary)^2*varx + varw + (A(i)*mu_c*varx/vary)^2*varv;
        vary = varc*varx+varv;
        xhat = [xhat y(1,t+1)*mu_c*varx/vary];
    end

    xtilde = x-xhat;
    
    xtilde_vec = linspace(min(xtilde),max(xtilde));
    data = zeros(length(xtilde_vec),1);
    for a = 1:length(xtilde_vec)
        data(a,1) = sum(xtilde < xtilde_vec(a));
    end
    data = data./length(xtilde); 
    subplot(1,2,1), hold on, plot(xtilde_vec,data,'Color',colors{i},'LineWidth',2);
    subplot(1,2,1), title('CDF of Estimation Error')
    h = legend('A=0.1','A=0.5','A=1','A=2','Location','Best');
    
    xtilde = (x-xhat).^2;
    xtilde_vec = linspace(min(xtilde),max(xtilde));
    data = zeros(length(xtilde_vec),1);
    for a = 1:length(xtilde_vec)
        data(a,1) = sum(xtilde < xtilde_vec(a));
    end
    data = data./length(xtilde); 
    subplot(1,2,2), hold on, plot(xtilde_vec,data,'Color',colors{i}, 'LineWidth',2);
    subplot(1,2,2), title('CDF of Estimation Error Squared')
end

suptitle(['Gireeja NonCoherence Paper Varying A VarC = ' num2str(varc)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,6],'PaperPosition',[0 0 12 6]);
% print('-dpdf','-r100','fig5_c1.pdf');

%%
%Kalman Filter with Additive Noise Try 2

clc;
clear all;

set(0,'defaulttextinterpreter','latex')

a = 1; %Test 2,0.5,10
c = 1;
varv = 10; %channel noise
varw = 1;
n = 500; %timesteps
M = 1000; %number of trials
Xtilde_ = [];
Xmtilde_ = [];

Xhat = [];
Xmem = [];
X = [];
Y = [];

for m=1:M
%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,sqrt(varv));
xhat = y*c*1/(c^2+varv);
xhatmem = y*c/(c^2+varv);
varx = 1;
vary = c^2*varx+varv;

%sigman = varx -2*c^2*varx^2/vary + c^2*varx^2/vary^2*(c^2*varx+varv);
sigman = (1-varx/vary)^2*varx + (varx/vary)^2*varv; %assumes c = 1
Sigman = [sigman];

%time passing
for t=1:(n-1)
    w = normrnd(0,sqrt(varw));
    x = [x a*x(1,t)+w];
    v = normrnd(0,sqrt(varv));
    y = [y c*x(1,(t+1))+v];
    
    sn = a^2*sigman+varw; %generating sn for this iteration
    kn = sn*c/(c^2*sn+varv);
    xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t+1)-c*a*xhatmem(1,t))];
    %sigman = (1-kn*c)^2 + kn^2*varv;
    sigman = (1-kn*c)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    
    Sigman = [Sigman sigman];
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
end

xtilde = x-xhat;
xmemtilde = x-xhatmem;

Xhat = [Xhat; xhat];
Xmem = [Xmem; xhatmem];
X = [X; x];
Y = [Y; y];

Xtilde_ = [Xtilde_ xtilde'];
Xmtilde_ = [Xmtilde_ xmemtilde'];

end

subplot(4,2,1), plot(0:t,mean(abs(Xmtilde_')),'b')
subplot(4,2,1), hold on, plot(0:t,mean(abs(Xtilde_')),'r')
subplot(4,2,1), title('Average Error'), xlabel('n = time'), ylabel('Magnitude')
legend('Xmtilde','Xtilde','Location','Best')

subplot(4,2,2), plot(0:t,var(Xmtilde_'),'b')
subplot(4,2,2), hold on, plot(0:t, var(Xtilde_'),'r',0:t,Sigman,'m')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}_{kf}$$','$$\tilde{X}$$','$$\tilde{X}calc$$','Location','Best')

Xtilde_ = Xtilde_.^2;
Xmtilde_ = Xmtilde_.^2;

subplot(4,2,3), hist(Xmtilde_(1,:),20)
subplot(4,2,3), title('Squared Error of $$\hat{X_m(0)}$$')

[vec, data] = cdfld(Xmtilde_(1,:)); [vec2, data2] = cdfld(Xtilde_(1,:));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2, data2, 'r', 'LineWidth',2)
legend('Xmtilde','Xtilde','Location','Best')
subplot(4,2,4), title('CDF of $$\hat{X}_m(0)$$ Squared Error'), xlabel('e = Magnitude'), ylabel('$$P((X(0)-\hat{X}_m(0))^2 < e)$$')

subplot(4,2,5), hist(Xmtilde_(n/2,:),20)
subplot(4,2,5), title(['Squared Error of $$\hat{X}_m$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xmtilde_(n/2,:)); [vec2, data2] = cdfld(Xtilde_(n/2,:));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2, data2, 'r', 'LineWidth',2)
subplot(4,2,6), title(['CDF of $$\hat{X}_m$$(' num2str(n/2) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n/2) '$$)-\hat{X}_m($$' num2str(n/2) '$$))^2 <$$ e)'])

subplot(4,2,7), hist(Xmtilde_(n,:),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}_m$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(n,:)); [vec2, data2] = cdfld(Xtilde_(n,:));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2, data2, 'r', 'LineWidth',2)
subplot(4,2,8), title(['CDF of $$\hat{X}_m$$(' num2str(n) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n) ')-$$\hat{X}_m$$(' num2str(n) '$$))^2 <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','kfadd_m1000.pdf');

%%
%Kalman Filter Ignoring Multiplicative Noise Try 2

clc;
clear all;

a = 1; %Test 2,0.5,10
% c = 1;
varv = 1/100; %multiplicative noise, 1/p, 0.01, 1, 10, 100
varw = 1;%system noise
n = 500;
M = 1000;
Xmtilde_ = [];

for m=1:M
    %initialize
    x = normrnd(0,1,[1,1]); %x[0]
    y = (1+normrnd(0,sqrt(varv)))*x;

    xhatmem = y;

    sigman = 0;

    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        x = [x a*x(1,t)+w];
        v = normrnd(0,sqrt(varv));
        y = [y (1+v)*x(1,(t+1))];

        sn = a^2*sigman+varw; %generating sn for this iteration
        kn = 1;
        xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t+1)-a*xhatmem(1,t))];
        %sigman = (1-kn*c)^2 + kn^2*varv;
        sigman = (1-kn)*sn; %generate new sigma_n for next iteration
        %sn = a^2*sigman+varw;

    end

    xmemtilde = x-xhatmem;
    
    Xmtilde_ = [Xmtilde_ xmemtilde'];
end

subplot(4,2,1), plot(0:t,mean(abs(Xmtilde_')),'b')
subplot(4,2,1), title('Average Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')

subplot(4,2,2), plot(0:t,var(Xmtilde_'),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')

Xmtilde_ = Xmtilde_.^2;

subplot(4,2,3), hist(Xmtilde_(1,:),20)
subplot(4,2,3), title('Squared Error of $$\hat{X_m(0)}$$')

[vec, data] = cdfld(Xmtilde_(1,:));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), title('CDF of $$\hat{X}_m(0)$$ Squared Error'), xlabel('e = Magnitude'), ylabel('$$P((X(0)-\hat{X}_m(0))^2 < e)$$')

subplot(4,2,5), hist(Xmtilde_(n/2,:),20)
subplot(4,2,5), title(['Squared Error of $$\hat{X}_m$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xmtilde_(n/2,:));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), title(['CDF of $$\hat{X}_m$$(' num2str(n/2) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n/2) '$$)-\hat{X}_m($$' num2str(n/2) '$$))^2 <$$ e)'])

subplot(4,2,7), hist(Xmtilde_(n,:),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}_m$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(n,:));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\hat{X}_m$$(' num2str(n) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n) ')-$$\hat{X}_m$$(' num2str(n) '$$))^2 <$$ e)'])

suptitle(['A = ' num2str(a) '; Power = ' num2str(1/varv) '; W = ' num2str(varw) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','kfavm_p100.pdf');

%%
%Kalman Filter Schenato Try 2
clc;
clear all;

n = 500;
M = 1000;

a = 0.99; %must be <1
%c = 1;
varv = 1; %additive noise
varw = 1; %system noise
p = 100; %multiplicative noise, test 0.01, 1, 10, 100
varx = 1;
vary = 1+varv;

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

subplot(4,2,1), plot(0:t,mean(abs(Xttilde')),'b')
subplot(4,2,1), hold on, plot(0:t,mean(abs(Xftilde')),'r'), plot(0:t,mean(abs(Xrtilde')),'g')
subplot(4,2,1), title('Average Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

subplot(4,2,2), plot(0:t,var(Xttilde'),'b')
subplot(4,2,2), hold on, plot(0:t, var(Xftilde'),'r')
subplot(4,2,2), hold on, plot(0:t, var(Xrtilde'),'g')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')

%Squared Error
Xttilde = Xttilde.^2;
Xrtilde = Xrtilde.^2;
Xftilde = Xftilde.^2;

subplot(4,2,3), hist(Xftilde(1,:),50)
subplot(4,2,3), title('Squared Error of $$\hat{Xf(0)}$$')

[vec, data] = cdfld(Xttilde(1,:)); [vec2, data2] = cdfld(Xftilde(1,:)); [vec3, data3] = cdfld(Xrtilde(1,:));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
legend({'$$\tilde{X}^t_{t|t}$$','$$\tilde{X}^r_{t|t}$$','$$\tilde{X}^r_{t|t-1}$$'},'Location','Best','Interpreter','Latex')
subplot(4,2,4), title('CDF of $$\hat{X}(0)$$ Squared Error'), xlabel('e = Magnitude'), ylabel('$$P((X(0)-\hat{X}(0))^2 < e)$$')

subplot(4,2,5), hist(Xftilde(n/2,:),50)
subplot(4,2,5), title(['Squared Error of $$\hat{X}f$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xttilde(n/2,:)); [vec2, data2] = cdfld(Xftilde(n/2,:)); [vec3, data3] = cdfld(Xrtilde(n/2,:));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
subplot(4,2,6), title(['CDF of $$\hat{X}$$(' num2str(n/2) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n/2) '$$)-\hat{X}($$' num2str(n/2) '$$))^2 <$$ e)'])

subplot(4,2,7), hist(Xftilde(n,:),50)
subplot(4,2,7), title(['Squared Error of $$\hat{X}f$$(' num2str(n) ')'])

[vec, data] = cdfld(Xttilde(n,:)); [vec2, data2] = cdfld(Xftilde(n,:)); [vec3, data3] = cdfld(Xrtilde(n,:));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2, data2, 'r', 'LineWidth',2), plot(vec3, data3, 'g', 'LineWidth',2)
subplot(4,2,8), title(['CDF of $$\hat{X}$$(' num2str(n) ') Squared Error']), xlabel('e = Magnitude'), ylabel(['P((X(' num2str(n) ')-$$\hat{X}$$(' num2str(n) '$$))^2 <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; P = ' num2str(p) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','schenato_p100.pdf');

%%
%Gireeja NonCoherence Try 2

%set(0,'defaultlegendlocation','Best')

clc
clear all

A = [0.1 0.5 0.9 1.2];
colors = {'r', 'b', 'g', 'k', 'm'};
%c = 1;
varx = 1;
varv = 1; %channel noise
varw = 1;
n = 100;
mu_c = 1;
varc = 1;
M = 500;

legentries = cell(1,length(A));

for i=1:length(A)
    Xtilde = [];
    for m=1:M
    %initialize
    c = normrnd(mu_c,sqrt(varc));
    x = normrnd(0,sqrt(varx),[1,1]); %x[0]
    y = c*x+normrnd(0,sqrt(varv));

    vary = varc*varx+varv;
    
    xhat = y*mu_c*varx/vary;

    %time passing
    for t=1:(n-1)
        c = normrnd(mu_c,sqrt(varc));
        x = [x A(i)*(x(1,t)-mu_c/(mu_c^2+varc)*y(1,t))+normrnd(0,sqrt(varw))]; %-xhat(1,t)
        y = [y c*x(1,(t+1))+normrnd(0,sqrt(varv))];

        %varx = A(i)^2*varx + varw;
        varx = A(i)^2*(1-mu_c^2/(mu_c^2+varc))*varx + varw + (A(i)*mu_c/(mu_c^2+varc))^2*varv;
        vary = varc*varx+varv;
        xhat = [xhat y(1,t+1)*mu_c*varx/vary];
    end
    
    Xtilde = [Xtilde (x-xhat)'];
    
    end
    
    legentries{1,i} = ['A=' num2str(A(1,i))];
    
    subplot(2,2,1), hold on, plot(0:t,mean(Xtilde'),'Color',colors{i});
    subplot(2,2,1), title('Mean Estimation Error')
    
    subplot(2,2,2), hold on, plot(0:t,var(Xtilde'),'Color',colors{i});
    subplot(2,2,2), title('Variance of Estimation Error')
    
    [vec, data] = cdfld(Xtilde(n,:)); 
    subplot(2,2,3), hold on, plot(vec,data,'Color',colors{i},'LineWidth',2);
    subplot(2,2,3), title('CDF of Estimation Error')
    
    Xtilde = Xtilde.^2;
    [vec, data] = cdfld(Xtilde(n,:));
    subplot(2,2,4), hold on, plot(vec,data,'Color',colors{i}, 'LineWidth',2);
    subplot(2,2,4), title('CDF of Estimation Error Squared')
end

subplot(2,2,1), legend(legentries,'Location','Best');
subplot(2,2,2), legend(legentries,'Location','Best');
subplot(2,2,4), legend(legentries,'Location','Best');
subplot(2,2,3), legend(legentries,'Location','Best');

suptitle(['Gireeja NonCoherence Paper Varying A M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig5_test_6.pdf');
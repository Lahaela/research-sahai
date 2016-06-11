% Kalman Filter Kalman Gain

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 0.99; %1 doesn't work?
b = 1;
%c = 1;
varc = 1;
mu_c = 1;
varv = 1;
varw = 1;

n = 1000;
M = 1000;

Xmtilde_ = []; % new kalman filter

KFGain = zeros(n,M);

for m=1:M
    varx = 1;
    vary = varc*varx+varv;
    %initialize
    x = normrnd(0,sqrt(varx),[1,1]); %x[0]
    c = normrnd(mu_c,sqrt(varc));
    y = c*x+normrnd(0,sqrt(varv));

    xhatm = y*mu_c*varx/vary;

    sigman_ = (1-c*mu_c*varx/vary)^2*varx+(mu_c*varx/vary)^2*varv;


    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,sqrt(varw))];
        c = normrnd(mu_c,sqrt(varc));
        y = [y c*x(1,(t+1))+normrnd(0,sqrt(varv))];
        
        varx = a^2*varx+b^2*varw;
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_c*sn_/(mu_c^2*sn_+varc*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_c*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_c)*sn_; %generating new sigma_n
        
        KFGain(1+t,m) = kf_;

    end
    
    Xmtilde_ = [Xmtilde_; x-xhatm];
end

subplot(1,2,1), hold on, plot(0:t,mean(KFGain),'m'), ylim([0, 0.13])
subplot(1,2,1), title('Kalman Filter Gain'), xlabel('n = time'), ylabel('Magnitude')

subplot(1,2,2), hold on, plot(0:t,mean(Xmtilde_.^2),'b','LineWidth',2)
subplot(1,2,2), title('MSE'), xlabel('n = time'), ylabel('Magnitude')

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','kfmul_a1_m1000_1.pdf');

%%
%UDelay Sweeping Controls
clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 1.05;
mu_c = 1;
varc = 1;
%c = 1;
varw = 1; %system noise
varv = 1; %observation noise

n = 30;
M = 1000;

k = 3; %delay

for d=0.52:0.02:0.6
%d = a^k*mu_c/(mu_c^2+varc);

    X = [];
    %Y = [];
    %U = [];

    for m=1:M
        mu_x = 0; %doesn't matter
        var_x = 1;

        %initialize
        xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
        c = normrnd(mu_c, sqrt(varc));
        yn = c*xideal + normrnd(0, sqrt(varv));
        un = -d*yn;

        %time passing
        for t=1:(n-1)
            mu_x = a*mu_x;

            w = normrnd(0,sqrt(varw));
            v = normrnd(0, sqrt(varv));
            c = normrnd(mu_c, sqrt(varc));

            if mod(t,k)==0
                xideal = [xideal a*xideal(1,t)+w+un(1,t-k+1)]; %+1 bc of Matlab 1-indexing
            else
                xideal = [xideal a*xideal(1,t)+w];
            end

            yn = [yn c*xideal(t+1)+v];
            un = [un -d*yn(t+1)];

        end

        X = [X; xideal];
        %Y = [Y; yn];
        %U = [U; un];
    end


    Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
    hold all
    plot(0:t, mean(Xsq))
end

title('MSE ($$X^2$$)'), xlabel('n = time'), ylabel('Magnitude')
% legend('0.52','0.54','0.56','0.58','0.6','Location','Best')


suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M)])

display(a^k*mu_c/(mu_c^2+varc))

set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
print('-dpdf','-r100','k3_fig4_2.pdf');

%%
%Packet Drops Theoretical

clc
clear all
close all
%set(0,'defaulttextinterpreter','latex')

Bound = [];

for L = 1:10

pe = 0.5;
alphap = 0.5; %n = alpha * L

n = L*alphap; %I control
m = 2*L; %I control
k = m-n;

es = (k)/(n);
boundexp = m*log(0.5*es+0.5) - log(es)*(k); %>= k

boundexp_2 = 2*L*log(0.5*(2*L/(alphap*L-1))) - ((2-alphap)*L+1)*log(((2-alphap)*L+1)/(alphap*L-1)); %>= (k+1)

% binbound = [];
% for i=k+1:m
%     binbound = [binbound nchoosek(m,i)*pe^(i)*(1-pe)^(m-i)];
% end

Bound = [Bound exp(boundexp)]; %exp(boundexp) or sum(binbound)
end

plot(1:10, Bound,'m','LineWidth',2)
%plot([2,4,6,8,10], Bound,'m','LineWidth',2)
title('Pe = 0.5, Alphap = 0.5 Bin Direct Computation')
xlabel('L'), ylabel('P(Drops > k) = P(communication breakdown)')

% set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
% print('-dpdf','-r100','bound_3.pdf');

%%
%Packet Drops Empirical

clc 
clear all
close all

L = 2;
pe = 0.5;
alphap = 0.5; %I control
trials = 100000;

n = L*alphap; %I control
m = 2*L; %I control
k = m-n;
breakdown = k/m;

y = zeros(m,trials);
x = rand(m,trials);
y(x<pe) = 1;

z = zeros(1,trials);
z(sum(y)/m > breakdown) = 1;

%bound = sum(y)/m;

mean(z)
    
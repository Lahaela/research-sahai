%Learning CRand

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 1; %
varx = 1;
varv = 1;
mu_c = 1;
varc = 1;

n = 100;
M = 1000;

X = [];
CRandHat = [];
for m = 1:M
    x = normrnd(0,sqrt(varx),[1,1]);
    v = normrnd(0,sqrt(varv));
    crand = normrnd(mu_c, sqrt(varc));
    y = crand*x+v;
    
    crandhat = [mu_c];
    for t = 1:(n-1) %remember t = past, t+1 = present
        if mod(t,2)==0
            u = 1;
        else
            u = 1;
        end
        
        x = [x a*x(t)+u];
        y = [y crand*x(t+1)+normrnd(0,sqrt(varv))];
        crandhat = [crandhat y(t+1)-a*y(t)];
        
        
    end
    crandhat = crandhat - crand; %computing the difference between reality
    X = [X; x];
    CRandHat = [CRandHat; crandhat];
end

subplot(1,2,1), plot(0:t, mean(abs(CRandHat)),'LineWidth',2)
subplot(1,2,1), title(['Abs Error in Learning Crand M = ' num2str(M)])

CRandHatMSE = mean(CRandHat.^2);

subplot(1,2,2), plot(0:t, CRandHatMSE, 'LineWidth',2)
subplot(1,2,2), title(['MSE in Learning CRand M = ' num2str(M)])
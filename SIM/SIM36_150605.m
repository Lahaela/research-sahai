%Summer 2015 Model

set(0,'defaulttextinterpreter','latex')

% clc;
% clear all;
% close all;
figure()
%1 packet = 1 bit
a = 1.22; %growth of the system
b = 1; %rate of info, bits / time
D = 3; %delay
alph = 2/3; %rate of channel
code = (D+1)*b; %bits in codeword
R = alph * code;
drop_p = 2^(-(code - R + 1)); %drop probability
thresh = 

varx = 1;
varw = 1;
varv = 1;
mu_c = 1;
varc = 2^(-2*R);

n = 50*D; %time of simulation
M = 1000; %number of trials
XD = []; %aggregate info from all trials

ud = a^D*mu_c/(mu_c^2+varc); %control for delay

%Empirical
for m = 1:M
    %initialize
    xd = normrnd(0,sqrt(varx),[1,1]); %x[0] ~ N(0, 1)
    c = normrnd(mu_c,sqrt(varc)); %c(0) ~ N(1, sigma_c)
    yd = c*xd+normrnd(0,sqrt(varv));
    
    for t = 1:(n-1)
        w = normrnd(0,sqrt(varw));
        v = normrnd(0,sqrt(varv));
        c = normrnd(mu_c,sqrt(varc));
        
        control = ud*yd(t-D);
        if abs(control) > thresh
            control = thresh*control/abs(control); %preserve sign
        end
        if mod(t, D+1)==0 && rand(1,1) > drop_p
            xd = [xd a*(xd(t)-ud*yd(t-D))+w];
        else
            xd = [xd a*xd(t)+w];
        end
        xd = [xd a*(xd(t)-ud*yd(t-D))+w];
        yd = [yd c*xd(t+1)+v];
    end
    
    XD = [XD; xd];
end
ZD = XD.^2;

%Theoretical Calculation
XDT = [1];
for t = 1:(n-1)
    %
end

%Plots
subplot(1,3,1), plot(0:t, mean(abs(XD)),'r', 'LineWidth',2)
subplot(1,3,1), title('Mean(Abs(X))'), xlabel('n = time')
%subplot(1,3,1), legend('XSystem [dd]','Location','Best')

subplot(1,3,2), hold on, plot(0:t,mean(ZD),'r', 'LineWidth',2)
%subplot(1,3,2), plot(0:t, XDT, 'm', 'LineWidth',2)
subplot(1,3,2), title('$$Mean(X^2)$$'), xlabel('n = time')
subplot(1,3,2), legend('Empirical','Theoretical','Location','Best');

[vec, data] = cdfld(ZD(:,n));
subplot(1,3,3), plot(log(vec),log(1-data'),'r', 'LineWidth',2)
subplot(1,3,3), title(['CCDF of ' '$$X^2($$' num2str(n) ')']), xlabel('Log X^2'), ylabel('Log CCDF')
%subplot(1,3,3), legend('XSystem [dd]','Location','Best')

suptitle(['A = ' num2str(a) ' ; Delay = ' num2str(D) ' ; Drop = ' num2str(drop_p) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

polyfit(log(vec(1,8500:9500)),log((1-data(8500:9500,1))'),1)
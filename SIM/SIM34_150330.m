%Spring 2014 Model
% Something is really really broken

set(0,'defaulttextinterpreter','latex')

% clc;
% clear all;
% close all;

figure()
a = 1.03; %growth of system
b = 1; %rate of information, bits per time
D = 3; %delay
alph = 2/3; %percent of code that's message
code = (D+1)*b; %bits in codeword
R = alph * code; %bits of information/message
drop_p = 2.^(-(code - R + 1)); %Reed Solomon codes

varx = 1; %this formulation doesn't need re-initialization
varv = 1; %additive observation noise
varw = 1;%system noise
mu_c = 1; %multiplicative/quantitative noise
varc = 2^(-2*R); %finding from Cover & Thomas--Rate Distortion

n = 50*D; %time per simulation
M = 1000; %number of trials

u = mu_c/(mu_c^2+varc); %control for ideal system
ud = 0; %conctrol for system with Delay

XIdeal = []; XD = [];

for m=1:M
    %initialize
    x = normrnd(0,sqrt(varx),[1,1]); %x[0] ~ N(0, 1)
    c = normrnd(mu_c,sqrt(varc)); %c(0) ~ N(1, sigma_c)
    y = c*x+normrnd(0,sqrt(varv));
    
    %delay + ideal start the same
    xd = x; yd = y;
    varxd = [varx];
    
    for t= 1:(n-1)
        w = normrnd(0,sqrt(varw)); v = normrnd(0,sqrt(varv));
        c = normrnd(mu_c,sqrt(varc));
        
        %No delay, no packet drops
        x = [x a*(x(t)-u*y(t))+w];
        y = [y c*x(t+1)+v];
        
        %delay
        if mod(t, D+1)==0 && rand(1,1) > drop_p
            ud = a^D*u;
        else
            ud = 0;
        end
        
        %Can't access negative indices
        try
            xd = [xd a*(xd(t)-ud*yd(t-D))+w];
            %varxd = [varxd a^2*(varxd(t)+(-2*ud*mu_c*a^D+ud^2*(mu_c^2+varc))*varxd(t-D)+ud^2*varv)+varw];
        catch
            xd = [xd a*xd(t)+w];
            %varxd = [varxd a^2*varxd(t)+varw];
        end
        
        yd = [yd c*xd(t+1)+v];
        
    end
    XIdeal = [XIdeal; x]; XD = [XD; xd];
    %Theoretical calculations
    XIT = [1]; XDT = [1]; %Ideal and Delay
    
    %don't need to re-initialize varx
    for t = 1:(n-1)
        varx = a^2*((1 - 2*u*mu_c + u^2*(mu_c^2+varc))*varx + u^2*varv) + varw;
        XIT = [XIT varx];
        
        g = 1 - drop_p;
        %XDT = [XDT 1];
        if mod(t,D+1) == 1 && t>D
            varxd = a^2*(XDT(t)-2*ud*g*mu_c*a^D*XDT(t-D)+ud^2*g*(mu_c^2+varc)*XDT(t-D)+ud^2*g*varv)+varw;
            XDT = [XDT varxd];
            %XDT = [XDT a^2*(XDT(t)+(-2*g*ud*mu_c*a^D+g*ud^2*(mu_c^2+varc))*XDT(t-D)+g*ud^2*varv)+varw];
        else
            XDT = [XDT a^2*XDT(t)+varw];
        end
    end
end

ZIdeal = XIdeal.^2; ZD = XD.^2;

%subplot(1,3,1), plot(0:t, mean(abs(XIdeal)),'b')
subplot(1,3,1), hold on, plot(0:t, mean(abs(XD)),'r')
subplot(1,3,1), title('Mean(Abs(X))'), xlabel('n = time')
subplot(1,3,1), legend('XSystem [dd]','Location','Best') % 'XIdeal',

%subplot(1,3,2), plot(0:t,mean(ZIdeal),'b')
%subplot(1,3,2), hold on, plot(0:t,XIT,'g')
subplot(1,3,2), hold on, plot(0:t,mean(ZD),'r',0:t,XDT,'m') %
subplot(1,3,2), title('$$Mean(X^2)$$'), xlabel('n = time')
subplot(1,3,2), legend('XSystem [dd]','X[dd]Theory','Location','Best'); %'XIdeal','XTheoretical',

%[vec, data] = cdfld(ZIdeal(:,n)); 
[vec2, data2] = cdfld(ZD(:,n));
%subplot(1,3,3), plot(log(vec),log(1-data'),'b')
subplot(1,3,3), hold on, plot(log(vec2),log(1-data2'),'r')
subplot(1,3,3), title(['CCDF of ' '$$X^2($$' num2str(n) ')']), xlabel('Log X^2'), ylabel('Log CCDF')
subplot(1,3,3), legend('XSystem [dd]','Location','Best') % 'XIdeal',

suptitle(['A = ' num2str(a) ' ; Delay = ' num2str(D) ' ; Drop = ' num2str(drop_p) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

polyfit(log(vec2(1,8500:9500)),log((1-data2(8500:9500,1))'),1)
%mean(ZD(:,[2:D:end])) %start with 1 for peak, start with 2 for low pt
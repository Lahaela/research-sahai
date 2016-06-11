%Blocklength

clc
clear all
hold all;

set(0,'defaulttextinterpreter','latex');
grid minor

% k = 2;
% a = 1.05;
mu_c = 1;
varc = 1;

%alph = 0.75; % n = alpha*L, blocklength = 2L
Alph = [0.25, 0.5, 0.75, 1];
%Alph = [0.5];
values = [];

for alph=Alph
    j = 1:0.01:10; %every jth timestep
    k = j-1; %actual delay
    % dropp = 1-(1./k); % dropp = (k-1)./k;
    % dropp = 1./k;
    %dropp = 2.^(-k);
    dropp = 2.^((alph-2)*k);
    % dropp = 0.5;

    % gammap = 0:0.01:1;
    % dropp = 1-gammap;

    z = (mu_c^2+varc)./(dropp.*mu_c^2+varc);
    a = z.^(1./(2*j));
    values = [values; a];
    % k = log(z)./(2*log(a));

    plot(k,a,'LineWidth',2)
    
end
xlabel('Delay K'), ylabel('Threshold for A')
title('Theoretical Effect of Delay on A, Drop = $$2^{-(2-alph)k}$$')
legend('0.25','0.5','0.75','1','Location','Best')

set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
print('-dpdf','-r100','150219_2.pdf');

%%
%Finding optimum alph-k pairs

%%
%Delay Gireeja Control Drops New Prob Model

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.17; %Test 1.1
b = 1; %rate of information through channel
D = 10; %delay
alph = 0.5;
msg = D*b;
R = alph * msg;

%varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varc = 2^(-2*R);
varv = 1; %additive observation noise
varw = 1;%system noise

%gammap = 1-2.^(-D);
%gammap = 1-2.^(R-msg);
drop_p = 2.^(R-msg);

%n = 200;
n = 5*D;
M = 100;

u = a*mu_c/(mu_c^2+varc);

XIdeal = [];
XSys = [];
for m=1:M
    varx = 1;
    vary = (varc+mu_c^2)*varx+varv;
    %initialize
    xideal = normrnd(0,sqrt(varx),[1,1]); %x[0]
    c = normrnd(mu_c,sqrt(varc));
    yideal = c*xideal+normrnd(0,sqrt(varv));

    xsys = xideal;
    ysys = yideal;
    varxsys = varx;
    varysys = vary;
    
    g = 1; %flag for packet drop

    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        c = normrnd(mu_c,sqrt(varc));
        v = normrnd(0,sqrt(varv));

        %baseline with no delay/packet drop
        xideal = [xideal a*xideal(1,t)+w-u*yideal(1,t)];
        yideal = [yideal c*xideal(1,(t+1))+v];

        %new varx
        varx = ((a^2-2*a*u*mu_c+u^2*(mu_c^2+varc))*varx + varw +u^2*varv);
        vary = (varc+mu_c^2)*varx+varv;

        %delay
        if mod(t,D) == 1 %packet gets sent at beginning of cycle
            if rand(1,1) < drop_p
                g = 0;
            else
                g = 1;
            end
            
            xsys = [xsys a*xsys(t)-(g*u)*ysys(t)+w];
            varxsys = ((a^2-2*a*(g*u)*mu_c+g*u^2*(mu_c^2+varc))*varxsys +varw +g*u^2*varv);
        else
            xsys = [xsys a*xsys(t)+w];
            varxsys = a^2*varxsys+varw;
        end
        ysys = [ysys c*xsys(t+1)+v];
        varysys = (varc+mu_c^2)*varxsys+varv;

    end

    XIdeal = [XIdeal; xideal];
    XSys = [XSys; xsys];
    
    %theoretical calculations of variance
    XIT = [1]; %X Ideal Theoretical
    XSysT = [1]; %X System Theoretical
    for t=1:(n-1)
        XIT = [XIT ((a^2-2*a*u*mu_c+u^2*(mu_c^2+varc))*XIT(t) +varw +u^2*varv)];
        
        if mod(t,D) == 1
            XSysT = [XSysT ((a^2-a*(1-drop_p)*u*mu_c)*XSysT(t) +varw +(1-drop_p)*u^2*varv)];
        else
            XSysT = [XSysT a^2*XSysT(t)+varw];
        end
    end
end

ZIdeal = XIdeal.^2;
ZSys = XSys.^2;

subplot(1,3,1), plot(0:t, mean(abs(XIdeal)),'b')
subplot(1,3,1), hold on, plot(0:t, mean(abs(XSys)),'r')
subplot(1,3,1), title('Mean(Abs(X))'), xlabel('n = time')
subplot(1,3,1), legend('XIdeal','XSystem [dd]','Location','Best')

subplot(1,3,2), plot(0:t,mean(ZIdeal),'b')
subplot(1,3,2), hold on, plot(0:t,XIT,'g')
subplot(1,3,2), hold on, plot(0:t,mean(ZSys),'r',0:t,XSysT,'m')
subplot(1,3,2), title('$$Mean(X^2)$$'), xlabel('n = time')
subplot(1,3,2), legend('XIdeal','XTheoretical','XSystem [dd]','X[dd]Theory','Location','Best');

[vec, data] = cdfld(ZIdeal(:,n)); [vec2, data2] = cdfld(ZSys(:,n));
subplot(1,3,3), plot(log(vec),log(1-data'),'b')
subplot(1,3,3), hold on, plot(log(vec2),log(1-data2'),'r')
subplot(1,3,3), title(['CCDF of ' '$$X^2($$' num2str(n) ')']), xlabel('Log X^2'), ylabel('Log CCDF')
subplot(1,3,3), legend('XIdeal','XSystem [dd]','Location','Best')


suptitle(['A = ' num2str(a) ' ; Delay = ' num2str(D) ' ; Drop = ' num2str(drop_p) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

polyfit(log(vec2(1,8500:9500)),log((1-data2(8500:9500,1))'),1)
mean(ZSys(:,[2:D:end])) %start with 1 for peak, start with 2 for low pt

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','sp15_a117_d50.pdf');
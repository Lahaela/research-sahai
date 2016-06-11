% Gireeja Control No Drops ReWrite

set(0,'defaulttextinterpreter','latex')

a = 3;

varv = 1;
varw = 1;
mu_c = 1;
varc = 1;

n = 100;
M = 10000;

X = [];

d = a*mu_c/(mu_c^2+varc);

for m=1:M
    x = normrnd(0,1,[1,1]);
    c = normrnd(mu_c,sqrt(varc));
    y = c*x + normrnd(0,sqrt(varv));
    
    for t=1:n-1
        w = normrnd(0,sqrt(varw));
        c = normrnd(mu_c,sqrt(varc));
        v = normrnd(0,sqrt(varv));
        
        x = [x (a*x(t)+w-d*y(t))];
        y = [y (c*x(t+1) +v)];
        
    end
    X = [X; x];
end

Z = X.^2;

figure(1), plot(0:t,mean(abs(X))), title('Mean(Abs(X))')
figure(2), plot(0:t,mean(Z)), title(['Mean(' '$$X^2$$' ')'])

[vec, data] = cdfld(Z(:,n));
figure(3), plot(log(vec),log(1-data)), xlim([0 8]), grid minor, title(['LogLog CCDF of' '$$X^2$$' '(100)'])

%%
%Delay Gireeja Control Drops

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.05; %Test 1.1
b = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
k = 5; %this is j in the theoretical curve!
%gammap = 1-(1/(k-1)); %probability a packet will go through
gammap = 0.5;

n = 100;
M = 1000;

d = a*mu_c/(mu_c^2+varc);

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
    
    g = 1;

    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        c = normrnd(mu_c,sqrt(varc));
        v = normrnd(0,sqrt(varv));

        xideal = [xideal a*xideal(1,t)+b*w-d*yideal(1,t)];
        yideal = [yideal c*xideal(1,(t+1))+v];

        %new varx
        varx = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx +b^2*varw +d^2*varv);
        vary = (varc+mu_c^2)*varx+varv;

        %delay
        if mod(t,k) == 0 %packet gets sent at beginning of cycle
            if rand(1,1) < gammap
                g = 1;
            else
                g = 0;
            end
            
            xsys = [xsys a*xsys(t)-(g*d)*ysys(t)+b*w];
            varxsys = ((a^2-2*a*(g*d)*mu_c+g*d^2*(mu_c^2+varc))*varxsys +b^2*varw +g*d^2*varv);
        else
            xsys = [xsys a*xsys(t)+b*w];
            varxsys = a^2*varxsys+b^2*varw;
        end
        ysys = [ysys c*xsys(t+1)+v];
        varysys = (varc+mu_c^2)*varxsys+varv;

    end

    XIdeal = [XIdeal; xideal];
    XSys = [XSys; xsys];
    
    %theoretical calculations
    XIT = [1]; %X Ideal Theoretical
    XSysT = [1]; %X System Theoretical
    for t=1:(n-1)
        XIT = [XIT ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*XIT(t) +b^2*varw +d^2*varv)];
        
        if mod(t,k) == 0
            XSysT = [XSysT ((a^2-a*gammap*d*mu_c)*XSysT(t) +b^2*varw +gammap*d^2*varv)];
        else
            XSysT = [XSysT a^2*XSysT(t)+b^2*varw];
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
subplot(1,3,3), title(['CCDF of ' '$$X^2($$' num2str(n) ')']), xlabel(['Log ' '$$X^2$$']), ylabel('Log CCDF')
subplot(1,3,3), legend('XIdeal','XSystem [dd]','Location','Best')


suptitle(['A = ' num2str(a) ' ; K = ' num2str(k) ' ; Drop = ' num2str(1-gammap) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

polyfit(log(vec2(1,8000:9000)),log((1-data2(8000:9000,1))'),1)

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','ddcontrol_base_p1.pdf');

%%
%Theoretical Curves A, K, P

%clc
%clear all

set(0,'defaulttextinterpreter','latex');
hold all;
figure()

% k = 2;
% a = 1.05;
mu_c = 1;
varc = 1;

j = 1:50; %every jth timestep
k = j-1; %actual delay
% dropp = 1-(1./k); % dropp = (k-1)./k;
%dropp = 1./k;
dropp = 2.^(-k);
%varc = 2.^(-2*k);
%dropp = 0.5;

% gammap = 0:0.01:1;
% dropp = 1-gammap;

z = (mu_c^2+varc)./(dropp.*mu_c^2+varc);
a = z.^(1./(2*j));
% k = log(z)./(2*log(a));

plot(k,a,'LineWidth',2)
xlabel('Delay K'), ylabel('Threshold for A')
title('Theoretical Effect of Delay on A, Drop = $$2^{(-k)}$$')

% set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
% print('-dpdf','-r100','dd_theoretical_xkya_4.pdf');
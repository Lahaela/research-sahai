%Delay Gireeja Control No Drops Ver 2

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.1; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
% k = 10;

n = 100;
M = 1000;

K = [1 2 3 4];

d = a*mu_c/(mu_c^2+varc);

for k=K
    XIdeal = [];% kalman filter
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

        XSysT = varxsys;
        XIT = varx;

        %time passing
        for t=1:(n-1)
            w = normrnd(0,sqrt(varw));
            c = normrnd(mu_c,sqrt(varc));
            v = normrnd(0,sqrt(varv));

            xideal = [xideal a*xideal(1,t)+b*w-d*yideal(1,t)];
            yideal = [yideal c*xideal(1,(t+1))+v];

            %new varx
            varx = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx +varw +d^2*varv);
            vary = (varc+mu_c^2)*varx+varv;

            %delay
            if mod(t,k) == 0
                xsys = [xsys a*xsys(t)-d*ysys(t)+b*w];
                varxsys = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varxsys +varw +d^2*varv);
            else
                xsys = [xsys a*xsys(t)+b*w];
                varxsys = a^2*varxsys+b^2*varw;
            end
            ysys = [ysys c*xsys(t+1)+v];
            varysys = (varc+mu_c^2)*varxsys+varv;

            XSysT = [XSysT varxsys];
            XIT = [XIT varx];
        end

        XIdeal = [XIdeal; xideal];
        XSys = [XSys; xsys];
    end
    
    i = find(k==K);
    subplot(2,2,i), plot(0:t,var(XIdeal),'b')
    subplot(2,2,i), hold on, plot(0:t,XIT,'g')
    subplot(2,2,i), hold on, plot(0:t,var(XSys),'r',0:t,XSysT,'m')
    subplot(2,2,i), title(['K = ' num2str(k)]), xlabel('n = time'), ylabel('Magnitude')

end

subplot(2,2,1), legend('$$X$$','XCalc','$$X_d$$','$$X_dCalc$$','Location','Best');
subplot(2,2,4), legend('$$X$$','XCalc','$$X_d$$','$$X_dCalc$$','Location','Best');

suptitle(['Variance of State A = ' num2str(a) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delay_a11_varyK.pdf');

%%
%Delay Gireeja Control No Drops Ignore this

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 0.5; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
% varw = 0.5;%system noise
k = 10;

n = 100;
M = 1000;
XIdeal = [];% kalman filter
XSys = [];

W = [0 0.5 1 10];

for varw=W
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

        XSysT = varxsys;
        XIT = varx;

        %time passing
        for t=1:(n-1)
            w = normrnd(0,sqrt(varw));
            c = normrnd(mu_c,sqrt(varc));
            v = normrnd(0,sqrt(varv));

            xideal = [xideal a*xideal(1,t)+b*w-a*mu_c/(mu_c^2+varc)*yideal(1,t)];
            yideal = [yideal c*xideal(1,(t+1))+v];

            %new varx
            alphaprime = mu_c/(mu_c^2+varc);
            varx = a^2*(1 -2*alphaprime*mu_c +alphaprime^2*(varc+mu_c^2))*varx +b^2*varw + (a*alphaprime)^2*varv;
            vary = (varc+mu_c^2)*varx+varv;

            %delay
            if mod(t,k) == 0
                alphasys = mu_c/(mu_c^2+varc);
                xsys = [xsys a*(xsys(t)-alphasys*ysys(t))+b*w];
                varxsys = a^2*(1 -2*alphasys*mu_c +alphasys^2*(varc+mu_c^2))*varxsys +b^2*varw + (a*alphasys)^2*varv;
            else
                xsys = [xsys a*xsys(t)+b*w];
                varxsys = a^2*varxsys+b^2*varw;
            end
            ysys = [ysys c*xsys(t+1)+v];
            varysys = (varc+mu_c^2)*varxsys+varv;

            XSysT = [XSysT varxsys];
            XIT = [XIT varx];
        end

        XIdeal = [XIdeal; xideal];
        XSys = [XSys; xsys];
    end

    i = find(varw==W);
    subplot(2,2,i), plot(0:t,var(XIdeal),'b')
    subplot(2,2,i), hold on, plot(0:t,XIT,'g')
    subplot(2,2,i), hold on, plot(0:t,var(XSys),'r',0:t,XSysT,'m')
    subplot(2,2,i), title('Variance of State'), xlabel('n = time'), ylabel('Magnitude')
    
end

subplot(2,2,1), legend('$$X$$','XCalc','$$X_d$$','$$X_dCalc$$','Location','Best');

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M) '; K = ' num2str(k)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delay_a05_wlim_k10.pdf');

%%
%Theoretical E[ |X(n)| ]

% a = 1; %Test 1.1
b = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1;
varw = 1;

n = 100;

A = [0.5 1 1.4 1.5];
for a=A

    varx = 1;
    absexpx = [sqrt(2/pi)];

    d = a*mu_c/(mu_c^2+varc);

    for t = 1:n-1
        varx = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx +varw +d^2*varv);
        vary = (varc+mu_c^2)*varx+varv;

        absexpx = [absexpx sqrt(varx*2/pi)]; 
    end
    i = find(a==A);
    subplot(2,2,i), plot(0:t, absexpx,'m')
    subplot(2,2,i), title(['A = ' num2str(a)]), xlabel('n = time')
end

suptitle(['E[|X(n)|]; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','expect_abs_x_varyA.pdf');

%%
%Empirical | X(n) |

%Delay Gireeja Control No Drops Ver 2

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
% close all;

% a = 1.1; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
% k = 10;

n = 100;
M = 1000;

A = [1.4 1.5 2 3];

for a=A
    d = a*mu_c/(mu_c^2+varc);
    XIdeal = [];
    for m=1:M
        varx = 1;
        vary = (varc+mu_c^2)*varx+varv;
        %initialize
        xideal = normrnd(0,sqrt(varx),[1,1]); %x[0]
        c = normrnd(mu_c,sqrt(varc));
        yideal = c*xideal+normrnd(0,sqrt(varv));

        %time passing
        for t=1:(n-1)
            w = normrnd(0,sqrt(varw));
            c = normrnd(mu_c,sqrt(varc));
            v = normrnd(0,sqrt(varv));

            xideal = [xideal a*xideal(1,t)+b*w-d*yideal(1,t)];
            yideal = [yideal c*xideal(1,(t+1))+v];

            %new varx
            varx = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx +varw +d^2*varv);
            vary = (varc+mu_c^2)*varx+varv;

        end

        XIdeal = [XIdeal; xideal];
    end
    z = XIdeal.^2;
    
    [vec, data] = cdfld(z(:,n));
    
    i = find(a==A);
    figure(1), subplot(2,2,i), plot(log(vec),log(1-data),'b')
%     subplot(2,2,i), hold on, plot(0:t,absexpx,'m')
    figure(1), subplot(2,2,i), title(['A = ' num2str(a)]), xlabel('Log X'), ylabel('Log CCDF'), xlim([0 100])
    polyfit(log(vec),log(1-data'),1)

end

% subplot(2,2,1), legend('$$|X(n)|$$','$$E[|X(n)|]$$','Location','Best');
% subplot(2,2,4), legend('$$|X(n)|$$','$$E[|X(n)|]$$','Location','Best');

suptitle(['CCDF of ' '$$X^2$$' '(100) ; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','LogLogCCDF_VaryA.pdf');

%%
%Delay Gireeja Control No Drops Ver 2

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.04; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
k = 3;

n = 100;
M = 10000;

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

    XSysT = varxsys;
    XIT = varx;

    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        c = normrnd(mu_c,sqrt(varc));
        v = normrnd(0,sqrt(varv));

        xideal = [xideal a*xideal(1,t)+b*w-d*yideal(1,t)];
        yideal = [yideal c*xideal(1,(t+1))+v];

        %new varx
        varx = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx +varw +d^2*varv);
        vary = (varc+mu_c^2)*varx+varv;

        %delay
        if mod(t,k) == 0
            xsys = [xsys a*xsys(t)-d*ysys(t)+b*w];
            varxsys = ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varxsys +varw +d^2*varv);
        else
            xsys = [xsys a*xsys(t)+b*w];
            varxsys = a^2*varxsys+b^2*varw;
        end
        ysys = [ysys c*xsys(t+1)+v];
        varysys = (varc+mu_c^2)*varxsys+varv;

        XSysT = [XSysT varxsys];
        XIT = [XIT varx];
    end

    XIdeal = [XIdeal; xideal];
    XSys = [XSys; xsys];
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

suptitle(['A = ' num2str(a) ' ; K = ' num2str(k) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

polyfit(log(vec2(1,8000:9000)),log((1-data2(8000:9000,1))'),1)

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','delay_control_comparison.pdf');

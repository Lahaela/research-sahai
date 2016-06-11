% Test of Stuff

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.5; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
k = 10;
gammap = 0.5; %probability packet will be dropped

n = 100;
M = 5000;
XIdeal = [];% kalman filter
XSys = [];

for m=1:M
    varx = 1;
    vary = (varc+mu_c^2)*varx+varv;
    %initialize
    xideal = normrnd(0,sqrt(varx),[1,1]); %x[0]
    c = normrnd(mu_c,sqrt(varc));
    yideal = c*xideal+normrnd(0,sqrt(varv));
    
    XIT = varx;
    
    %time passing
    for t=1:(n-1)
        xideal = [xideal a*(xideal(1,t)-mu_c*varx/vary*yideal(t))+b*normrnd(0,sqrt(varw))];
        c = normrnd(mu_c,sqrt(varc));
        yideal = [yideal c*xideal(1,(t+1))+normrnd(0,sqrt(varv))];
        
        %new varx
%         varx = a^2*(1-2*mu_c^2*varx/vary+(mu_c*varx/vary)^2*varc)*varx+b^2*varw+(a*mu_c*varx/vary)^2*varv;
%         varx = a^2*varx+b^2*varw;
%         varx = a^2*(1 -2*(1+mu_c)^2*varx/vary +((1+mu_c)*varx/vary)^2*(1+2*mu_c+varc))*varx +b^2*varw +(a*(1+mu_c)*varx/vary)^2*varv;
        alphaprime = mu_c*varx/vary;
        varx = a^2*(1 -2*alphaprime*mu_c +alphaprime^2*(varc+mu_c^2))*varx +b^2*varw + (a*alphaprime)^2*varv;
        vary = (varc+mu_c^2)*varx+varv;

        XIT = [XIT varx];
    end
    
    XIdeal = [XIdeal; xideal];
end

subplot(1,2,1), plot(0:t,mean(abs(XIdeal)), 'b')
subplot(1,2,1), title('Average Abs(State)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

subplot(1,2,2), plot(0:t,var(XIdeal),'b')
subplot(1,2,2), hold on, plot(0:t,XIT,'g')
subplot(1,2,2), title('Variance of State'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M) '; K = ' num2str(k) '; GammaP = ' num2str(gammap)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','test_leah_a15.pdf');
%%
%Delay Leah Control No Drops

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 0.5; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise
k = 10;
gammap = 0.5; %probability packet will be dropped

n = 100;
M = 1000;
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
        
        xideal = [xideal a*xideal(1,t)+b*w-a*mu_c*varx/vary*yideal(1,t)];
        yideal = [yideal c*xideal(1,(t+1))+v];
        
        %new varx
        alphaprime = mu_c*varx/vary;
        varx = a^2*(1 -2*alphaprime*mu_c +alphaprime^2*(varc+mu_c^2))*varx +b^2*varw + (a*alphaprime)^2*varv;
        vary = (varc+mu_c^2)*varx+varv;

        %delay
        if mod(t,k) == 0
            alphasys = mu_c*varxsys/varysys;
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

subplot(1,2,1), plot(0:t,mean(abs(XIdeal)), 'b')
subplot(1,2,1), hold on, plot(0:t,mean(abs(XSys)),'r')
subplot(1,2,1), title('Average Abs(State)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

subplot(1,2,2), plot(0:t,var(XIdeal),'b')
subplot(1,2,2), hold on, plot(0:t,var(XSys),'r',0:t,XSysT,'m',0:t,XIT,'g')
subplot(1,2,2), title('Variance of State'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

% XIdeal = XIdeal.^2;
% Xtildereal = Xtildereal.^2;
% 
% subplot(4,2,3), hist(Xtildereal(:,1),20)
% subplot(4,2,3), title('Squared Error of $$\hat{X}(0)$$')
% 
% [vec, data] = cdfld(Xmtilde_(:,1)); [vec2, data2] = cdfld(Xtildereal(:,1));
% subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
% subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')
% 
% subplot(4,2,5), hist(Xtildereal(:,n/2),20)
% subplot(4,2,5), title(['Squared Error of $$\hat{X}$$(' num2str(n/2) ')'])
% 
% [vec, data] = cdfld(Xmtilde_(:,n/2)); [vec2, data2] = cdfld(Xtildereal(:,n/2));
% subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
% legend('$$\tilde{X}_{kf}$$','$$\tilde{X}$$','Location','Best');
% 
% subplot(4,2,7), hist(Xtildereal(:,n),20)
% subplot(4,2,7), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])
% 
% [vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xtildereal(:,n));
% subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M) '; K = ' num2str(k) '; GammaP = ' num2str(gammap)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delaydrop_u100.pdf');

%%
%Test Gireeja NonCoherence Control
set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.5; %Test 1.1
b = 1;
%c = 1;
varc = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_c = 1;
varv = 0; %additive observation noise
varw = 0;%system noise

n = 100;
M = 1000;
XIdeal = [];% kalman filter

for m=1:M
    varx = 1;
    vary = (varc+mu_c^2)*varx+varv;
    %initialize
    xideal = normrnd(0,sqrt(varx),[1,1]); %x[0]
    c = normrnd(mu_c,sqrt(varc));
    yideal = c*xideal+normrnd(0,sqrt(varv));
    
    XIT = varx;
    
    d = a*mu_c/(mu_c^2+varc);
    %time passing
    for t=1:(n-1)
%         xideal = [xideal a*(xideal(1,t)-mu_c/(mu_c^2+varc)*yideal(t))+b*normrnd(0,sqrt(varw))];
        xideal = [xideal (a*xideal(t) +b*normrnd(0,sqrt(varw)) -d*yideal(t))];
        c = normrnd(mu_c,sqrt(varc));
        yideal = [yideal c*xideal(1,(t+1))+normrnd(0,sqrt(varv))];
        
        %new varx
%         varx = a^2*(1-2*mu_c^2*varx/vary+(mu_c*varx/vary)^2*varc)*varx+b^2*varw+(a*mu_c*varx/vary)^2*varv;
%         varx = a^2*varx+b^2*varw;
%         varx = a^2*(1 -2*(1+mu_c)^2*varx/vary +((1+mu_c)*varx/vary)^2*(1+2*mu_c+varc))*varx +b^2*varw +(a*(1+mu_c)*varx/vary)^2*varv;
%         alphaprime = mu_c/(mu_c^2+varc);
%         varx = a^2*(1 -2*alphaprime*mu_c +alphaprime^2*(varc+mu_c^2))*varx +b^2*varw + (a*alphaprime)^2*varv;
        varx = (a^2 -2*a*d*mu_c + d^2*(mu_c^2+varc))*varx +b^2*varw + d^2*varv;
        vary = (varc+mu_c^2)*varx+varv;

        XIT = [XIT varx];
    end
    
    XIdeal = [XIdeal; xideal];
end

subplot(1,2,1), plot(0:t,mean(abs(XIdeal)), 'b')
subplot(1,2,1), title('Average Abs(State)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

subplot(1,2,2), plot(0:t,var(XIdeal),'b')
subplot(1,2,2), hold on, plot(0:t,XIT,'g')
subplot(1,2,2), title('Variance of State'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M)])

figure(2),plot(var(XIdeal))

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','test_gireeja_a15.pdf');

%%
%Delay Gireeja Control No Drops

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
varw = 10;%system noise
k = 10;

n = 100;
M = 5000;
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

subplot(1,2,1), plot(0:t,mean(abs(XIdeal)), 'b')
subplot(1,2,1), hold on, plot(0:t,mean(abs(XSys)),'r')
subplot(1,2,1), title('Average Abs(State)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$X$$','$$X_d$$','Location','Best');

subplot(1,2,2), plot(0:t,var(XIdeal),'b')
subplot(1,2,2), hold on, plot(0:t,XIT,'g')
subplot(1,2,2), hold on, plot(0:t,var(XSys),'r',0:t,XSysT,'m')
subplot(1,2,2), title('Variance of State'), xlabel('n = time'), ylabel('Magnitude')
subplot(1,2,2), legend('$$X$$','XCalc','$$X_d$$','$$X_dCalc$$','Location','Best');

subplot(1,2,1), legend('$$X$$','$$X_d$$','Location','Best');

% XIdeal = XIdeal.^2;
% Xtildereal = Xtildereal.^2;
% 
% subplot(4,2,3), hist(Xtildereal(:,1),20)
% subplot(4,2,3), title('Squared Error of $$\hat{X}(0)$$')
% 
% [vec, data] = cdfld(Xmtilde_(:,1)); [vec2, data2] = cdfld(Xtildereal(:,1));
% subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
% subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')
% 
% subplot(4,2,5), hist(Xtildereal(:,n/2),20)
% subplot(4,2,5), title(['Squared Error of $$\hat{X}$$(' num2str(n/2) ')'])
% 
% [vec, data] = cdfld(Xmtilde_(:,n/2)); [vec2, data2] = cdfld(Xtildereal(:,n/2));
% subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
% legend('$$\tilde{X}_{kf}$$','$$\tilde{X}$$','Location','Best');
% 
% subplot(4,2,7), hist(Xtildereal(:,n),20)
% subplot(4,2,7), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])
% 
% [vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xtildereal(:,n));
% subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
% subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
% subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M) '; K = ' num2str(k)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delay_a05_wlim_k10.pdf');

%%
%Test Gireeja NonCoherence Ver 2

set(0,'defaulttextinterpreter','latex')

clc;
clear all;
close all;

a = 1.5;
varv = 1;
varw = 1;
mu_c = 1;
varc = 1;

n = 100;
M = 5000;

d = a*mu_c/(mu_c^2+varc);
XIdeal = [];

for m=1:M
    varx = [1];

    x = normrnd(0,1,[1,1]);
    y = normrnd(mu_c,sqrt(varc))*x +normrnd(0,sqrt(varv));

    for t = 1:(n-1)
        x = [x a*x(t)+normrnd(0,sqrt(varw))-d*y(t)];
        y = [y (normrnd(mu_c,sqrt(varc))*x(t+1) +normrnd(0,sqrt(varv)))];

        varx = [varx ((a^2-2*a*d*mu_c+d^2*(mu_c^2+varc))*varx(t) +varw +d^2*varv)];

    end
    XIdeal = [XIdeal; x];
end
% subplot(2,2,find(varw==W)), plot(0:t,var(XIdeal,1),'b',0:t,varx,'m')
% subplot(2,2,find(varw==W)), title(['VarW = ' num2str(varw)])
% 
% 
% subplot(2,2,4), legend('Empirical','Theoretical','Location','Best')
% subplot(2,2,1), legend('Empirical','Theoretical','Location','Best')

subplot(2,3,1), plot(0:t,mean(abs(XIdeal)),'b')
subplot(2,3,1), title('Average Abs(State)')

subplot(2,3,2), plot(0:t,var(XIdeal,1),'b')
subplot(2,3,2), hold on, plot(0:t,varx,'m')
subplot(2,3,2), title('Variance of State'), legend('Empirical','Theoretical','Location','Best')
% 
subplot(2,3,3), plot(0:t,var(XIdeal))
subplot(2,3,3), title('Variance of XIdeal')

XIdealV = XIdeal.^2;

[vec, data] = cdfld(XIdealV(:,1));
subplot(2,3,4), plot(log(vec),log(1-data),'LineWidth',2)
subplot(2,3,4), title('CCDF of XIdeal(0)')
polyfit(log(vec),(log(1-data))',1)

[vec2, data2] = cdfld(XIdealV(:,n/2));
subplot(2,3,5), plot(log(vec2),log(1-data2),'LineWidth',2)
subplot(2,3,5), title(['CCDF of XIdeal(' num2str(n/2) ')'])
polyfit(log(vec2),(log(1-data2))',1)

[vec3, data3] = cdfld(XIdealV(:,n));
subplot(2,3,6), plot(log(vec3),log(1-data3),'LineWidth',2)
subplot(2,3,6), title(['CCDF of XIdeal(' num2str(n) ')'])
polyfit(log(vec3),(log(1-data3))',1)

suptitle(['Variance of State; A = ' num2str(a) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','gireeja_140802_3.pdf');
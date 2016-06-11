% Delay Estimation No Packet Drops

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.1; %Test 1.1
b = 1;
%c = 1;
varu = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_u = 1;
varv = 0; %additive observation noise
varw = 0;%system noise
k = 10;

n = 50;
M = 500;
Xmtilde_ = [];% kalman filter
Xtildereal = [];

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,varx,[1,1]); %x[0]
    u = normrnd(mu_u,varu);
    y = u*x+normrnd(0,varv);

    xhatm = y*mu_u*varx/vary;
    xhatreal = [0]; %mu_x

    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,varw)];
        u = normrnd(mu_u,varu);
        y = [y u*x(1,(t+1))+normrnd(0,varv)];
        
        varx = a^2*varx+b^2*varw;
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u^2*sn_+varu*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n
        
        if mod(t+1,k) == 0
            xhatreal = [xhatreal xhatm(t+1)];
        else
            xhatreal = [xhatreal a*xhatreal(t)];
        end
    end
    
    Xmtilde_ = [Xmtilde_; x-xhatm];
    Xtildereal = [Xtildereal; x-xhatreal];
end

subplot(4,2,1), plot(0:t,mean(abs(Xtildereal)), 'r')
subplot(4,2,1), hold on, plot(0:t,mean(abs(Xmtilde_)),'b')
subplot(4,2,1), title('Average Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X_{kf}}$$','Location','Best');

subplot(4,2,2), plot(0:t,var(Xtildereal),'r')
subplot(4,2,2), hold on, plot(0:t,var(Xmtilde_),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X_{kf}}$$','Location','Best');

Xmtilde_ = Xmtilde_.^2;
Xtildereal = Xtildereal.^2;

subplot(4,2,3), hist(Xtildereal(:,1),20)
subplot(4,2,3), title('Squared Error of $$\hat{X}(0)$$')

[vec, data] = cdfld(Xmtilde_(:,1)); [vec2, data2] = cdfld(Xtildereal(:,1));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')

subplot(4,2,5), hist(Xtildereal(:,n/2),20)
subplot(4,2,5), title(['Squared Error of $$\hat{X}$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xmtilde_(:,n/2)); [vec2, data2] = cdfld(Xtildereal(:,n/2));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
legend('$$\tilde{X_{kf}}$$','$$\tilde{X}$$','Location','Best');

subplot(4,2,7), hist(Xtildereal(:,n),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xtildereal(:,n));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig.pdf');

%%
%Delay Estimation + No Info

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.1; %Test 1.1
b = 1;
%c = 1;
varu = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_u = 1;
varv = 0; %additive observation noise
varw = 0;%system noise
k = 10;

n = 50;
M = 500;
Xmtilde_ = [];% kalman filter
Xtildereal = [];
Xtildenoinfo = [];

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,varx,[1,1]); %x[0]
    u = normrnd(mu_u,varu);
    y = u*x+normrnd(0,varv);

    xhatm = y*mu_u*varx/vary;
    xhatreal = [0]; %mu_x
    xhatnoinfo = xhatm; %you get the first value)

    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,varw)];
        u = normrnd(mu_u,varu);
        y = [y u*x(1,(t+1))+normrnd(0,varv)];
        
        varx = a^2*varx+b^2*varw;
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u^2*sn_+varu*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n
        
        if mod(t,k-1) == 0
            xhatreal = [xhatreal xhatm(t+1)];
        else
            xhatreal = [xhatreal a*xhatreal(t)];
        end
        
        xhatnoinfo = [xhatnoinfo a*xhatnoinfo(t)];
    end
    
    Xmtilde_ = [Xmtilde_; x-xhatm];
    Xtildereal = [Xtildereal; x-xhatreal];
    Xtildenoinfo = [Xtildenoinfo; x-xhatnoinfo];
end

subplot(1,2,1), plot(0:t,mean(abs(Xtildereal)), 'r')
subplot(1,2,1), hold on, plot(0:t,mean(abs(Xmtilde_)),'b',0:t,mean(abs(Xtildenoinfo)),'m')
subplot(1,2,1), title('Average Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','$$\tilde{X}_{no}$$','Location','Best');

subplot(1,2,2), plot(0:t,var(Xtildereal),'r')
subplot(1,2,2), hold on, plot(0:t,var(Xmtilde_),'b',0:t,var(Xtildenoinfo),'m')
subplot(1,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','$$\tilde{X}_{no}$$','Location','Best');

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig3_no_v10.pdf');

%%
%Delay Estimation Packet Drops

%***This may have been lost during a Force Quit***
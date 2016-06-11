%Delay Estimation Drop Packets
%***NEED TO CHECK***

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
gammap = 0.5; %probability packet will be dropped
g=1;

n = 50;
M = 1000;
Xmtilde_ = [];% kalman filter
Xtildereal = [];

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,sqrt(varx),[1,1]); %x[0]
    u = normrnd(mu_u,sqrt(varu));
    y = u*x+normrnd(0,sqrt(varv));

    xhatm = y*mu_u*varx/vary;
    xhatdrop = xhatm; %begins with same estimate
    xhatreal = [0]; %mu_x

    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;
    sigmandrop = sigman_; %starts out with same estimate

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,sqrt(varw))];
        u = normrnd(mu_u,sqrt(varu));
        y = [y u*x(1,(t+1))+normrnd(0,sqrt(varv))];
        
        varx = a^2*varx+b^2*varw;
        
        % ideal case
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u^2*sn_+varu*varx+varv);        
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n
        
        %packet drops
        sndrop = a^2*sigmandrop+b^2*varw;
        if mod(t,k) == 0
            if rand(1,1) > gammap
                g = 1;
            else
                g = 0;
            end
        end
        kfdrop = g*mu_u*sndrop/(mu_u^2*sndrop+varu*varx+varv);
        
        xhatdrop = [xhatdrop a*xhatdrop(1,t)+kfdrop*(y(1,t+1)-mu_u*a*xhatdrop(1,t))];
        sigmandrop = (1-kfdrop*mu_u)*sndrop; %generating new sigma_n
        
%         ['g = ' num2str(g) '; kfdrop = ' num2str(kfdrop) '; sigmandrop = ' num2str(sigmandrop) '; sndrop = ' num2str(sndrop)]
        
        %delay
        if mod(t,k) == 0
            xhatreal = [xhatreal xhatdrop(t+1)];
        else
            xhatreal = [xhatreal a*xhatreal(t)];
        end
    end
    
    Xmtilde_ = [Xmtilde_; x-xhatm];
    Xtildereal = [Xtildereal; x-xhatreal];
end

subplot(4,2,1), plot(0:t,mean(abs(Xtildereal)), 'r')
subplot(4,2,1), hold on, plot(0:t,mean(abs(Xmtilde_)),'b')
subplot(4,2,1), title('Average Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','Location','Best');

subplot(4,2,2), plot(0:t,var(Xtildereal),'r')
subplot(4,2,2), hold on, plot(0:t,var(Xmtilde_),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','Location','Best');

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
legend('$$\tilde{X}_{kf}$$','$$\tilde{X}$$','Location','Best');

subplot(4,2,7), hist(Xtildereal(:,n),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xtildereal(:,n));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M) '; K = ' num2str(k) '; GammaP = ' num2str(gammap)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delaydrop_u100.pdf');

%%
%Delay Control

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
gammap = 0.5; %probability packet will be dropped

n = 50;
M = 5000;
Xmtilde_ = [];% kalman filter
Xtildereal = [];

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    xideal = normrnd(0,sqrt(varx),[1,1]); %x[0]
    xsys = xideal;
    u = normrnd(mu_u,sqrt(varu));
    y = u*xideal+normrnd(0,sqrt(varv));
    ysys = yideal;

    xhatm = y*mu_u*varx/vary;
    xhatdrop = xhatm; %begins with same estimate
    xhatreal = [0]; %mu_x
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        u = normrnd(mu_u,sqrt(varu));
        v = normrnd(0,sqrt(varv));
        xideal = [xideal a*xideal(1,t)+b*w];
        y = [y u*xideal(1,(t+1))+v];
        
        varx = a^2*varx+b^2*varw;
        vary = varu*varx+varv;
        
        %packet drops
        sndrop = a^2*sigmandrop+b^2*varw;
        if rand(1,1) > gammap
            kfdrop = mu_u*sndrop/(mu_u^2*sndrop+varu*varx+varv);
        else
            kfdrop = 0;
        end
        
        xhatdrop = [xhatdrop a*xhatdrop(1,t)+kfdrop*(y(1,t+1)-mu_u*a*xhatdrop(1,t))];
        sigmandrop = (1-kfdrop*mu_u)*sndrop; %generating new sigma_n
        
        %delay
        if mod(t,k) == 0
            xhatreal = [xhatreal xhatdrop(t+1)];
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
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','Location','Best');

subplot(4,2,2), plot(0:t,var(Xtildereal),'r')
subplot(4,2,2), hold on, plot(0:t,var(Xmtilde_),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X}$$','$$\tilde{X}_{kf}$$','Location','Best');

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
legend('$$\tilde{X}_{kf}$$','$$\tilde{X}$$','Location','Best');

subplot(4,2,7), hist(Xtildereal(:,n),20)
subplot(4,2,7), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xtildereal(:,n));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M) '; K = ' num2str(k) '; GammaP = ' num2str(gammap)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delaydrop_u100.pdf');
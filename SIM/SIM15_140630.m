% Kalman Filter Multiplicative Noise Rajasekaran

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.0; %Test 2,0.5,10
b = 1;
%c = 1;
varu = 1; %multiplicative noise, 1/p, 0.01, 0.1, 1, 100
mu_u = 1;
varv = 1; %additive observation noise
varw = 1;%system noise

n = 1000;
M = 1000;
Xmtilde = []; %ignoring mul noise
Xmtilde_ = []; % new kalman filter

%Sigman_ = zeros(n,M); %test of whether u is necessary, numbers converge

for m=1:M
    varx = 1;
    vary = varu*varx+varv;
    %initialize
    x = normrnd(0,varx,[1,1]); %x[0]
    u = normrnd(mu_u,varu);
    y = u*x+normrnd(0,varv);

    xhata = y*varx/vary;
    xhatm = y*mu_u*varx/vary;

    sigman = varx -2*varx^2/vary + varx^2/vary^2*(varx+varv);
    sigman_ = (1-u*mu_u*varx/vary)^2*varx+(mu_u*varx/vary)^2*varv;
    
    %Sigman_(1,m) = sigman_;

    %time passing
    for t=1:(n-1)
        x = [x a*x(1,t)+b*normrnd(0,varw)];
        u = normrnd(mu_u,varu);
        y = [y u*x(1,(t+1))+normrnd(0,varv)];

        sn = a^2*sigman+varw; %generating sn for this iteration
        kn = sn/(sn+varv);
        xhata = [xhata a*xhata(1,t)+kn*(y(1,t+1)-a*xhata(1,t))];
        %sigman = (1-kn*c)^2 + kn^2*varv;
        sigman = (1-kn)*sn; %generate new sigma_n for next iteration
        %sn = a^2*sigman+varw;
        
        varx = a^2*varx+b^2*varw;
        
        sn_ = a^2*sigman_+b^2*varw;
        kf_ = mu_u*sn_/(mu_u*sn_+varu*varx+varv);
        xhatm = [xhatm a*xhatm(1,t)+kf_*(y(1,t+1)-mu_u*a*xhatm(1,t))];
        sigman_ = (1-kf_*mu_u)*sn_; %generating new sigma_n
        
        %Sigman_(1+t,m) = sigman_;

    end
    
    Xmtilde = [Xmtilde; x-xhata];
    Xmtilde_ = [Xmtilde_; x-xhatm];
end

subplot(4,2,1), plot(0:t,mean(Xmtilde), 'r')
subplot(4,2,1), hold on, plot(0:t,mean(Xmtilde_.^2),'b')
subplot(4,2,1), title('Average Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X_a}$$','$$\tilde{X_m}$$','Location','Best');

subplot(4,2,2), plot(0:t,var(Xmtilde),'r')
subplot(4,2,2), hold on, plot(0:t,var(Xmtilde_),'b')
subplot(4,2,2), title('Variance of Error'), xlabel('n = time'), ylabel('Magnitude')
legend('$$\tilde{X_a}$$','$$\tilde{X_m}$$','Location','Best');

Xmtilde = Xmtilde.^2;
Xmtilde_ = Xmtilde_.^2;

subplot(4,2,3), hist(Xmtilde_(:,1),M/5)
subplot(4,2,3), title('Squared Error of $$\hat{X}_m(0)$$')

[vec, data] = cdfld(Xmtilde_(:,1)); [vec2, data2] = cdfld(Xmtilde(:,1));
subplot(4,2,4), plot(vec, data, 'b', 'LineWidth',2)
subplot(4,2,4), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,4), title('CDF of $$\tilde{X}^2(0)$$'), xlabel('e = Magnitude'), ylabel('$$P(\tilde{X}^2(0) < e)$$')

subplot(4,2,5), hist(Xmtilde_(:,n/2),M/5)
subplot(4,2,5), title(['Squared Error of $$\hat{X}_m$$(' num2str(n/2) ')'])

[vec, data] = cdfld(Xmtilde_(:,n/2)); [vec2, data2] = cdfld(Xmtilde(:,n/2));
subplot(4,2,6), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,6), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,6), title(['CDF of $$\tilde{X}^2$$(' num2str(n/2) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n/2) '$$) <$$ e)'])
legend('$$\tilde{X_m}$$','$$\tilde{X_a}$$','Location','Best');

subplot(4,2,7), hist(Xmtilde_(:,n),M/5)
subplot(4,2,7), title(['Squared Error of $$\hat{X}_m$$(' num2str(n) ')'])

[vec, data] = cdfld(Xmtilde_(:,n)); [vec2, data2] = cdfld(Xmtilde(:,n));
subplot(4,2,8), plot(vec, data, 'b','LineWidth',2)
subplot(4,2,8), hold on, plot(vec2,data2,'r','LineWidth',2)
subplot(4,2,8), title(['CDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varu) '; W = ' num2str(varw) '; U = ' num2str(varu) '; M = ' num2str(M)])

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','kfmul_a1_m1000_1.pdf');
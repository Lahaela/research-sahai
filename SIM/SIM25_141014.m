%%
%Tangent lines to x^2

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

c = 1;

xprime = [0:100:1500];
y = xprime.^2;
yprime_alpha = 2*c*xprime;
yprime_beta = -yprime_alpha.*xprime + y;
yprime = yprime_alpha.*xprime + yprime_beta;

%plot(x,y,'b','LineWidth',2)
%plot(x,yprime,'m','LineWidth',2)

%%
%Genie gives tangent line estimation

set(0,'defaulttextinterpreter','latex')

a = 1.1; %Test 1.1
c = 1;
varv = 1; %additive observation noise
varw = 1;%system noise

n = 50;
M = 10000;
Xhattilde = [];% estimate
X = [];
Y = [];

for m=1:M
    mu_x = 0;
    var_x = 1;
    %secondmom_x = mu_x^2 + var_x;
    %mu_y = 2*c*mu_x; %!change
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    xideal = abs(xideal);
    yn = c*xideal.^2 + normrnd(0,sqrt(varv)); %note *xideal to be a matrix !change
    tmp = abs(y-yn(1));
    [idx idx] = min(tmp);
    
    if yprime_alpha(idx)==0
        xhat = [mu_x];
    else
        xhat = [(yn(1)-yprime_beta(idx))/yprime_alpha(idx)]; %cov(x,y)/var(y) * (y-E[y]) + E[x]; !change
        %xhat = [xprime(idx)];
    end
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        v = normrnd(0,sqrt(varv));
        xideal = [xideal a*xideal(1,t)+abs(w)];
        yn = [yn c*(xideal(1,(t+1)))^2+v];
        tmp = abs(y-yn(t+1));
        [idx idx] = min(tmp);
        
        mu_x = a*mu_x;
        %secondmom_x = a^2*secondmom_x + varw;
        %var_x = secondmom_x - mu_x^2;
        %mu_y = 2*c*mu_x; %!change
        
        if yprime_alpha(idx)==0
            xhat = [xhat mu_x];
        else
            xhat = [xhat (yn(t+1)-yprime_beta(idx))/yprime_alpha(idx)];
            %xhat = [xhat xprime(idx)];
        end
    end
    
    Xhattilde = [Xhattilde; xideal-xhat];
    X = [X; xideal];
    Y = [Y; yn];
end

subplot(2,2,1), plot(0:t,mean(abs(Xhattilde)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xhattilde.^2),'b','LineWidth',2)
subplot(2,2,2), title('Mean Squared Error'), xlabel('n = time'), ylabel('Magnitude')

Xhattilde = Xhattilde.^2;

subplot(2,2,3), hist(Xhattilde(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xhattilde(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(c) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','nonlinear_tangent_3.pdf');

%%
%Square root non-linear solution

set(0,'defaulttextinterpreter','latex')

a = 2; %Test 1.1
c = 1;
varv = 10; %additive observation noise
varw = 10;%system noise

n = 100;
M = 5000;
Xhattilde = [];% estimate
X = [];
Y = [];

for m=1:M
    mu_x = 5;
    var_x = 1;
    %secondmom_x = mu_x^2 + var_x;
    %mu_y = 2*c*mu_x; %!change
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    %xideal = abs(xideal);
    yn = c*xideal.^2 + normrnd(0,sqrt(varv)); %note *xideal to be a matrix !change
    
    xhat = (yn/c).^(0.5);
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        v = normrnd(0,sqrt(varv));
        xideal = [xideal a*xideal(1,t)+w];
        yn = [yn c*(xideal(1,(t+1)))^2+v];
        
        mu_x = a*mu_x;

        xhat = [xhat (yn(t+1)/c)^(0.5)];
    end
    
    Xhattilde = [Xhattilde; xideal-xhat];
    X = [X; xideal];
    Y = [Y; yn];
end

subplot(2,2,1), plot(0:t,mean(abs(Xhattilde)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(Error)'), xlabel('n = time'), ylabel('Magnitude')

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xhattilde.^2),'b','LineWidth',2)
subplot(2,2,2), title('Mean Squared Error'), xlabel('n = time'), ylabel('Magnitude')

Xhattilde = Xhattilde.^2;

subplot(2,2,3), hist(Xhattilde(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xhattilde(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(c) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','delaydrop_u100.pdf');

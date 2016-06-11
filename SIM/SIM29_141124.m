%Controller noise uniform

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 1.65;
b1 = 0;
b2 = 10;

n = 50;
M = 10000;

X = [];
Y = [];
U = [];

for m=1:M
    mu_x = 10;
    var_x = 0;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    yn = xideal;
    un = -2/(b1+b2)*yn;
    
    %time passing
    for t=1:(n-1)
        mu_x = a*mu_x;

        b = unifrnd(b1,b2);

        xideal = [xideal a*(xideal(1,t)+b*un(t))];
        
        yn = [yn xideal(t+1)];
        un = [un -2/(b1+b2)*yn(t+1)];
        
    end
    
    X = [X; xideal];
    Y = [Y; yn];
    U = [U; un];
end

subplot(2,2,1), semilogy(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xsq),'b','LineWidth',2)
subplot(2,2,2), title('MSE ($$X^2$$)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

baseline = logspace(10^-10,1);
absX = abs(X);

[vec, data] = cdfld(absX(:,n));
subplot(2,2,4), hold on, loglog(vec, 1-data, 'b','LineWidth',2)
%subplot(2,2,4), hold on, loglog(baseline, 1./baseline, 'm', 'LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$Abs X$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])
%subplot(2,2,4), xlim([

suptitle(['A = ' num2str(a) '; B1 = ' num2str(b1) '; B2 = ' num2str(b2) '; M = ' num2str(M)])

polyfit(log(vec(1,8000:9000)),log((1-data(8000:9000,1))'),1)

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig2.pdf');

%%
%Controller noise normal

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 2;
mu_b = 5.5;
varb = 1/12;

n = 50;
M = 1000;

X = [];
Y = [];
U = [];

for m=1:M
    mu_x = 0;
    var_x = 1;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    yn = xideal;
    un = -1/mu_b*yn;
    
    %time passing
    for t=1:(n-1)
        mu_x = a*mu_x;

        b = normrnd(mu_b, sqrt(varb));

        xideal = [xideal a*(xideal(1,t)+b*un(t))];
        
        yn = [yn xideal(t+1)];
        un = [un -1/mu_b*yn(t+1)];
        
    end
    
    X = [X; xideal];
    Y = [Y; yn];
    U = [U; un];
end

subplot(2,2,1), semilogy(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xsq),'b','LineWidth',2)
subplot(2,2,2), title('MSE ($$X^2$$)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

baseline = logspace(10^-10,1);
absX = abs(X);

[vec, data] = cdfld(absX(:,n));
subplot(2,2,4), hold on, loglog(vec, 1-data, 'b','LineWidth',2)
%subplot(2,2,4), hold on, loglog(baseline, 1./baseline, 'm', 'LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$Abs X$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])
%subplot(2,2,4), xlim([

suptitle(['A = ' num2str(a) '; B1 = ' num2str(b1) '; B2 = ' num2str(b2) '; M = ' num2str(M)])

polyfit(log(vec(1,8000:9000)),log((1-data(8000:9000,1))'),1)

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','udelay_6.pdf');
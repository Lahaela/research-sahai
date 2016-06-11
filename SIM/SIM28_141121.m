%Udelay + control

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 1.04;
mu_c = 1;
varc = 1;
%c = 1;
varw = 0; %system noise
varv = 0; %observation noise

n = 200;
M = 10000;

k = 10; %delay
d = a^k*mu_c/(mu_c^2+varc);
%d = a^k/c;

X = [];
Y = [];
U = [];

for m=1:M
    mu_x = 0; %doesn't matter
    var_x = 1;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    c = normrnd(mu_c, sqrt(varc));
    yn = c*xideal + normrnd(0, sqrt(varv));
    %un = a^k*(mu_x + (c*secondmom_x/(c^2*secondmom_x + varv))*(yn-c*mu_x));
    un = -d*yn;
    
    %time passing
    for t=1:(n-1)
        mu_x = a*mu_x;
        
        w = normrnd(0,sqrt(varw));
        v = normrnd(0, sqrt(varv));
        c = normrnd(mu_c, sqrt(varc));
        
        if mod(t,k)==0
            xideal = [xideal a*xideal(1,t)+w+un(1,t-k+1)]; %+1 bc of Matlab 1-indexing
        else
            xideal = [xideal a*xideal(1,t)+w];
        end
        
        yn = [yn c*xideal(t+1)+v];
        un = [un -d*yn(t+1)];
        
    end
    
    X = [X; xideal];
    Y = [Y; yn];
    U = [U; un];
end

subplot(2,2,1), plot(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xsq),'b','LineWidth',2)
subplot(2,2,2), title('MSE ($$X^2$$)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xsq(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','fig1.pdf');

%%
%Ydelay + control

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 1.01;
mu_c = 1;
varc = 1;
varw = 1; %system noise
varv = 1; %observation noise

n = 50;
M = 1000;

k = 4; %delay + 1
d = a^k*mu_c/(mu_c^2+varc);
%d = a^k/c;

X = [];
Y = [];
U = [];

for m=1:M
    mu_x = 0;
    var_x = 1;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    c = normrnd(mu_c, sqrt(varc));
    yn = [0];
    un = -d*yn;
    
    %time passing
    for t=1:(n-1)
        mu_x = a*mu_x;
        
        w = normrnd(0,sqrt(varw));
        v = normrnd(0, sqrt(varv));
        c = normrnd(mu_c, sqrt(varc));
        
        if mod(t,k+1)==0 %action happens 1 timestep after
            xideal = [xideal a*(xideal(1,t)+un(1,t))+w]; %+1 bc of Matlab 1-indexing
        else
            xideal = [xideal a*xideal(1,t)+w];
        end
        
        if t<k
            yn = [yn 0];
        else
            yn = [yn c*xideal(t+1-k)+v];
        end
        
        un = [un -d*yn(t+1)];
        
    end
    
    X = [X; xideal];
    Y = [Y; yn];
    U = [U; un];
end

subplot(2,2,1), plot(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

subplot(2,2,2), plot(0:t, mean(Xsq),'b','LineWidth',2)
subplot(2,2,2), title('MSE ($$X^2$$)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xsq(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(varc) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','fig5.pdf');
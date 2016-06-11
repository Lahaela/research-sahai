% U delay with control, no additive observation noise

set(0,'defaulttextinterpreter','latex')

a = 1.5;
c = 1;
varw = 1; %system noise

n = 50;
M = 1000;
k = 2; %delay

X = [];

for m=1:M
    mu_x = 1;
    var_x = 1;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    yn = c*xideal;
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        
        if t < k
            xideal = [xideal a*xideal(1,t)+w];
        else
            xideal = [xideal a*xideal(1,t)+w-xideal(1,t-k+1)]; %+1 bc of Matlab 1-indexing
        end
        
        yn = [yn c*xideal(t+1)];
        
        mu_x = a*mu_x;
        
    end
    
    X = [X; xideal];
    %Y = [Y; yn];
end

subplot(2,2,1), plot(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xsq.^2),'b','LineWidth',2)
subplot(2,2,2), title('Mean Squared Error (X)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xsq(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(c) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','udelay_7.pdf');

%%
%Udelay + control + additive noise ***IN PROGRESS***

set(0,'defaulttextinterpreter','latex')

a = 1;
c = 1;
varw = 0; %system noise
varw = 0;

n = 30;
M = 1000;
k = 2; %delay

X = [];
Y = [];

for m=1:M
    mu_x = 0;
    var_x = 1;
    secondmom_x = var_x + mu_x^2;
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    yn = c*xideal + normrnd(0, sqrt(varv));
    un = mu_x + (c*secondmom_x/(c^2*secondmom_x + varv))*(yn-c*mu_x);
    
    %time passing
    for t=1:(n-1)
        mu_x = a*mu_x;
        
        w = normrnd(0,sqrt(varw));
        v = normrnd(0, sqrt(varv));
        
        if t < k
            xideal = [xideal a*xideal(1,t)+w];
            secondmom_x = a^2*secondmom_x+var_w;
        else
            xideal = [xideal a*xideal(1,t)+w-xideal(1,t-k+1)]; %+1 bc of Matlab 1-indexing
            secondmom_x = ; %involved, need to find recursive pattern
        end
        
        yn = [yn c*xideal(t+1)+v];
        
        
    end
    
    X = [X; xideal];
    Y = [Y; yn];
end

subplot(2,2,1), plot(0:t,mean(abs(X)), 'b','LineWidth',2)
subplot(2,2,1), title('Mean Abs(X)'), xlabel('n = time'), ylabel('Magnitude')

Xsq = X.^2;

%subplot(2,2,2), plot(0:t,var(Xhattilde),'b','LineWidth',2)
subplot(2,2,2), plot(0:t, mean(Xsq.^2),'b','LineWidth',2)
subplot(2,2,2), title('Mean Squared Error (X)'), xlabel('n = time'), ylabel('Magnitude')

subplot(2,2,3), hist(Xsq(:,n),20)
subplot(2,2,3), title(['Squared Error of $$\hat{X}$$(' num2str(n) ')'])

[vec, data] = cdfld(Xsq(:,n));
subplot(2,2,4), loglog(vec, 1-data, 'b','LineWidth',2)
subplot(2,2,4), title(['Loglog CCDF of $$\tilde{X}^2$$(' num2str(n) ')']), xlabel('e = Magnitude'), ylabel(['$$P(\tilde{X}^2($$' num2str(n) '$$) <$$ e)'])

suptitle(['A = ' num2str(a) '; V = ' num2str(varv) '; W = ' num2str(varw) '; C = ' num2str(c) '; M = ' num2str(M)])

% set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
% print('-dpdf','-r100','udelay_6.pdf');
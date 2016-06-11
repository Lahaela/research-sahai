%%
%Nonlinear observation Linear estimator y=cx2(n)

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.1; %Test 1.1
c = 1;
varv = 0; %additive observation noise
varw = 0;%system noise

n = 50;
M = 5000;
Xhattilde = [];% estimate

for m=1:M
    mu_x = 10;
    var_x = 1;
    secondmom_x = mu_x^2 + var_x;
    thirdmom_x = mu_x^3 + 3*mu_x*var_x;
    fourthmom_x = mu_x^4 + 6*mu_x^2*var_x + 3*var_x^2;
    
    mu_y = c*secondmom_x; %!change
    secondmom_y = c^2*fourthmom_x + varv; %!change
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    y = c*xideal.^2 + normrnd(0,sqrt(varv)); %note *xideal to be a matrix !change

    xhat = (thirdmom_x*c - mu_x*mu_y)/(secondmom_y-mu_y^2)*(y-mu_y) + mu_x; %cov(x,y)/var(y) * (y-E[y]) + E[x]; !change
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        v = normrnd(0,sqrt(varv));
        xideal = [xideal a*xideal(1,t)+w];
        y = [y c*(xideal(1,(t+1)))^2+v];
        
        mu_x = a*mu_x;
        secondmom_x = a^2*secondmom_x + varw;
        var_x = secondmom_x - mu_x^2;
        thirdmom_x = mu_x^3 + 3*mu_x*var_x;
        fourthmom_x = mu_x^4+6*mu_x^2*var_x+3*var_x^2;

        mu_y = c*secondmom_x;
        secondmom_y = c^2*fourthmom_x + varv;
        
        xhat = [xhat ((thirdmom_x*c - mu_x*mu_y)/(secondmom_y-mu_y^2)*(y(1,t+1)-mu_y) + mu_x)];

    end
    
    Xhattilde = [Xhattilde; xideal-xhat];
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

%%
%Nonlinear observation Linear estimator no memory y = cx(n) check

set(0,'defaulttextinterpreter','latex')

clc;
clear all;

a = 1.1; %Test 1.1
c = 1;
varv = 10; %additive observation noise
varw = 10;%system noise

n = 50;
M = 5000;
Xhattilde = [];% estimate

for m=1:M
    mu_x = 1;
    var_x = 1;
    secondmom_x = mu_x^2 + var_x;
    thirdmom_x = mu_x^3 + 3*mu_x*var_x;
    fourthmom_x = mu_x^4 + 6*mu_x^2*var_x + 3*var_x^2;
    
    mu_y = c*mu_x; %!change
    secondmom_y = c^2*secondmom_x + varv; %!change
    
    %initialize
    xideal = normrnd(mu_x,sqrt(var_x),[1,1]); %x[0]
    y = c*xideal + normrnd(0,sqrt(varv)); %note *xideal to be a matrix !change

    xhat = (secondmom_x*c - mu_x*mu_y)/(secondmom_y-mu_y^2)*(y-mu_y) + mu_x; %cov(x,y)/var(y) * (y-E[y]) + E[x]; !change
    
    %time passing
    for t=1:(n-1)
        w = normrnd(0,sqrt(varw));
        v = normrnd(0,sqrt(varv));
        xideal = [xideal a*xideal(1,t)+w];
        y = [y c*(xideal(1,(t+1)))+v];
        
        mu_x = a*mu_x;
        secondmom_x = a^2*secondmom_x + varw;
        var_x = secondmom_x - mu_x^2;
        thirdmom_x = mu_x^3 + 3*mu_x*var_x;
        fourthmom_x = mu_x^4+6*mu_x^2*var_x+3*var_x^2;

        mu_y = c*mu_x;
        secondmom_y = c^2*secondmom_x + varv;
        
        xhat = [xhat ((secondmom_x*c - mu_x*mu_y)/(secondmom_y-mu_y^2)*(y(1,t+1)-mu_y) + mu_x)];

    end
    
    Xhattilde = [Xhattilde; xideal-xhat];
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
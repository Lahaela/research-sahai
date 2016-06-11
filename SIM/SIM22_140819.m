%YDelay

clc
clear all
close all

set(0,'defaulttextinterpreter','latex')

a = 2;

varx = 1;
varw = 1;
varv = 0;
mu_c = 1;
varc = 0;

M = 10000;
n = 10;
k = 2;
j = k+1;

X = [];

for m = 1:M
    x = normrnd(0,varx,[1,1]);
    y = [0];
    for t = 1:(n-1)
        w = normrnd(0,varw);
        
        flag = k;
        u = 0;
        sign = -1;
        for i = 0:(t-k)
            if flag == k
                sign = -1;
                flag = flag+1;
            elseif flag == k+1;
                    sign = 1;
                    flag = 1;
            elseif (0 < flag < k);
                sign = 0;
                flag = flag+1;
            end
                
            u = u+sign*a^(j+i)*y(t-i);
        end

        x = [x a*x(t)+w+u];
        
        if t>=k
            c = normrnd(mu_c,varc);
            y = [y c*x(t-k+1)+normrnd(0,varv)];
        else
            y = [y 0];
        end
        
        
    end
    X = [X; x];
end

subplot(1,3,1), plot(0:t,mean(abs(X)),'LineWidth',2)
subplot(1,3,1), title('Abs(X)'), xlabel('Time')

Z = X.^2;

subplot(1,3,2), plot(0:t,mean(Z),'LineWidth',2)
subplot(1,3,2), title('$$X^2$$'), xlabel('Time')

[vec, data] = cdfld(Z(:,n));
subplot(1,3,3), plot(log(vec),log(1-data),'LineWidth',2)
ccdfslope = polyfit(log(vec(1,9000:9900)),log((1-data(9000:9900,1))'),1);
subplot(1,3,3), title(['CCDF of X(' num2str(n) ')']), xlabel(['Slope (9000:9900) = ' num2str(ccdfslope(1))])


suptitle(['A = ' num2str(a) ' ; K = ' num2str(k) '; W = ' num2str(varw) '; V = ' num2str(varv) '; C = ' num2str(varc) '; M = ' num2str(M)])
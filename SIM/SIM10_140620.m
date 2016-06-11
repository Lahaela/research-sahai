%%
% Empirical Variance of X[n]

clc
clear all

%a = 2;
A = [1 2 5 10];
colors = {'r', 'b', 'g', 'k', 'm'};
c=1;
varv = 1;
varw = 1;
n = 1000;

X_n = [];

for a=A
    for i=1:100
        %initialize
        x = 0.5; %x[0]
        y = c*x+normrnd(0,varv);
        xhat = y*c*1/(c^2+varv);

        varx = 0;
        vary = c^2+varv;

        %time passing
        for t=1:(n-1)
            x = [x a*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
            y = [y c*x(1,(t+1))+normrnd(0,varv)];

            varx = [varx a^2*varx(1,t)*varv/vary(1,t) + varw];
            vary = [vary c^2*varx(1,t+1) + varv];
            xhat = [xhat y(1,t+1)*c*varx(1,t+1)/vary(1,t+1)];
        end

        X_n = [X_n x'];
    end

    Variance = var(X_n');
    i = find(a==A);
    subplot(2,2,i), plot(Variance,'b')
    subplot(2,2,i), hold on, plot(varx,'Color','r','LineWidth',2)
    subplot(2,2,i), title(['A = ' num2str(a)])
end
% h = legend('W = 0.5','W = 1','W = 5','W = 10','Location','Best');
% xlabel('m = values of x')
% ylabel('P(x < m)')
suptitle('Var(X) Control Problem Varying A')
%Control Attempt 1
clc
clear all

a = 1; %Test 2,0.5,10
c = 1;
varv = 1; %channel noise
varw = 1;
n = 100;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);

varx = 1;
vary = c^2+varv;

%time passing
for t=1:(n-1)
    x = [x a*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
    y = [y c*x(1,(t+1))+normrnd(0,varv)];
    
    varx = a^2*(1-c^2*varx/vary)^2*varx + varw + (a*c*varx/vary)^2*varv;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
end

xtilde = x-xhat;

subplot(1,3,1), plot(0:t,x,'k')
subplot(1,3,1), hold on, plot(0:t,y,'b')
subplot(1,3,1), hold on, plot(0:t,xhat,'r')
xlabel('n=time')
legend1 = legend('x','xhat','xhatmem','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(1,3,2), hold on, plot(0:t,xtilde,'k')
%subplot(1,2,2), hold on, plot(0:t, 1./sqrt(1:500),'r')
xlabel('n=time')
legend2 = legend('xhat-x');
subplot(1,3,2), title('Error')

subplot(1,3,3), hist(xtilde)
subplot(1,3,3), title(['Error over timesteps n= ' num2str(n)])

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

sum(abs(xhat-x))/n
sum(abs(x))/n

%'Error'
%a*((1-c^2*varx/vary)^2*varx + (c*varx/vary)^2*varv)

%%
%Varying A

clc
clear all

A = [0.1 0.5 1 2];
colors = {'r', 'b', 'g', 'k', 'm'};
c = 1;
varv = 1; %channel noise
varw = 1;
n = 500;

for i=1:length(A)
    %initialize
    x = normrnd(0,1,[1,1]); %x[0]
    y = c*x+normrnd(0,varv);
    xhat = y*c*1/(c^2+varv);

    varx = 1;
    vary = c^2+varv;

    %time passing
    for t=1:(n-1)
        x = [x A(i)*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
        y = [y c*x(1,(t+1))+normrnd(0,varv)];

        varx = A(i)^2*(1-c^2*varx/vary)^2*varx + varw + (A(i)*c*varx/vary)^2*varv;
        vary = c^2*varx + varv;
        xhat = [xhat y(1,t+1)*c*varx/vary];
    end

    xtilde = x-xhat;

    %subplot(2,2,find(a==A)), hist(xtilde);
    %subplot(2,2,find(a==A)), title(['a = ' num2str(a)])
    %subplot(2,2,find(a==A)), xlim([-3 3])
    %set(get(gca,'child'),'FaceColor','none','EdgeColor',Edges(1,find(a==A)));
    
    xtilde_vec = linspace(min(xtilde),max(xtilde));
    data = zeros(length(xtilde_vec),1);
    for a = 1:length(xtilde_vec)
        data(a,1) = sum(xtilde < xtilde_vec(a));
    end
    data = data./length(xtilde); 
    plot(xtilde_vec,data,'Color',colors{i},'LineWidth',2.5); hold on;
end

h = legend('A=1','A=2','A=5','A=10','Location','Best');

%%
%Varying C

clc
clear all

a = 1;
colors = {'r', 'b', 'g', 'k', 'm'};
C = [1 2 5 10];
%V = [0.5 1 5 10]; %channel noise
varv = 1;
varw = 1;
n = 1000;

for c=C
    %initialize
    x = normrnd(0,1,[1,1]); %x[0]
    y = c*x+normrnd(0,varv);
    xhat = y*c*1/(c^2+varv);

    varx = 1;
    vary = c^2+varv;

    %time passing
    for t=1:(n-1)
        x = [x a*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
        y = [y c*x(1,(t+1))+normrnd(0,varv)];

        varx = a^2*(1-c^2*varx/vary)^2*varx + varw + (a*c*varx/vary)^2*varv;
        vary = c^2*varx + varv;
        xhat = [xhat y(1,t+1)*c*varx/vary];
    end

    xtilde = abs(x-xhat);
    
    xtilde_vec = linspace(min(xtilde),max(xtilde));
    data = zeros(length(xtilde_vec),1);
    for i = 1:length(xtilde_vec)
        data(i,1) = sum(xtilde < xtilde_vec(i));
    end
    data = data./length(xtilde); 
    plot(xtilde_vec,data,'Color',colors{find(c==C)},'LineWidth',2.5); hold on;
end

h = legend('C = 1','C = 2','C = 5','C = 10','Location','Best');
xlabel('m = magnitude of error')
ylabel('P(abs(xtilde) < m)')
title('CDF of Error X-Xhat Control Problem')
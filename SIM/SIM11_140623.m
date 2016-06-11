%Control Attempt 2
clc
clear all

a = 2; %Test 2,0.5,10
c = 1;
varv = 1; %channel noise
varw = 1;
n = 100;
P = [1 2 5 10];

for p=P
%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);

varx = 1;
vary = c^2+varv;

%time passing
for t=1:(n-1)
    x = [x a*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
    y = [y (1+normrnd(0,1/p))*x(1,(t+1))+normrnd(0,varv)];
    
    %varx = a^2*(1-c^2*varx/vary)^2*varx + varw + (a*c*varx/vary)^2*varv;
    varx = a^2*varx*varv/vary + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
end

xtilde = xhat-x;
i = find(p==P);

subplot(4,2,2*i-1), plot(0:t,xtilde,'k')
xlabel('n=time')
legend2 = legend('xhat-x');
subplot(4,2,2*i-1), title(['Error, p = ' num2str(p)])

subplot(4,2,2*i), hist(xtilde)
subplot(4,2,2*i), title(['Error over timesteps n= ' num2str(n) '; p = ' num2str(p)])

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

end

set(gcf,'PaperUnits','inches','PaperSize',[12,12],'PaperPosition',[0 0 12 12]);
print('-dpdf','-r100','fig1');

%%
%Control Attempt 2
clc
clear all

a = 2; %Test 2,0.5,10
c = 1;
varv = 1; %channel noise
varw = 1;
n = 50;
P = [2 3 5 10];
colors = {'r','b','k','g','m'};

for p=P
%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = y*c*1/(c^2+varv);

varx = 1;
vary = c^2+varv;

%time passing
for t=1:(n-1)
    x = [x a*(x(1,t)-xhat(1,t))+normrnd(0,varw)];
    y = [y (1+normrnd(0,1/p))*x(1,(t+1))+normrnd(0,varv)];
    
    %varx = a^2*(1-c^2*varx/vary)^2*varx + varw + (a*c*varx/vary)^2*varv;
    varx = a^2*varx*varv/vary + varw;
    vary = c^2*varx + varv;
    xhat = [xhat y(1,t+1)*c*varx/vary];
    
end

xtilde = xhat-x;

xtilde_vec = linspace(min(xtilde),max(xtilde));
data = zeros(length(xtilde_vec),1);
for i = 1:length(xtilde_vec)
    data(i,1) = sum(xtilde < xtilde_vec(i));
end
data = data./length(xtilde); 

i = find(p==P);
plot(xtilde_vec,data,'Color',colors{i},'LineWidth',2); hold on;

xlabel('n=time')
legend('p = 2','p = 3','p = 5','p = 10','Location','Best');

title('CDF Multiplicative Noise Varying P')

end

set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
print('-dpdf','-r100','fig2');
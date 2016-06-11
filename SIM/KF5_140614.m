%%
%Predicting the future
clc
clear all

a = 1; %Test 2,0.5,10
c = 1;
varv = 0.09; %channel noise
varw = 0.04;
n = 50;

%initialize
x = normrnd(0,1,[1,1]); %x[0]
y = c*x+normrnd(0,varv);
xhat = [0];
xhatmem = [0];

varx = 1;
vary = c^2+varv;

%sn = a^2+varw - 2*a^2*c^2/(c^2+varv) + a^2*c^4/(c^2+varv)^2 + a^2*c^2/(c^2+varv)^2*varv;

%time passing
for t=1:(n-1)
    w = normrnd(0,varw);
    x = [x a*x(1,t)+w];
    v = normrnd(0,varv);
    y = [y c*x(1,(t+1))+v];
    
    sn = varx -2*c*varx^2/vary + c^2*varx^2/vary^2*(c^2*varx+varv);
    %sn = a^2*sn+varw; %generating sn for this iteration
    kn = a*sn*c/(c^2*sn+varv);
    xhatmem = [xhatmem a*xhatmem(1,t)+kn*(y(1,t)-c*xhatmem(1,t))];
    %sn = (1-kn*c)*sn; %generate new sigma_n for next iteration
    %sn = a^2*sigman+varw;
    xhat = [xhat y(1,t)*a*c*varx/vary];
    
    varx = a^2*varx + varw;
    vary = c^2*varx + varv;
    
end

%hold all
subplot(1,2,1), stem(0:t,x,'k','filled')
%subplot(1,2,1), hold on, stem(0:t,y,'b','filled')
subplot(1,2,1), hold on, stem(0:t,xhat,'r')
subplot(1,2,1), hold on, stem(0:t, xhatmem, 'g','filled')
xlabel('n=time')
legend1 = legend('x','y','xhat','Location','Best');

%subplot(1,2,2), plot(0:t,abs(y/c-x),'g')
subplot(1,2,2), hold on, plot(0:t,abs(xhat-x),'r')
subplot(1,2,2), hold on, plot(0:t,abs(xhatmem-x),'b')
xlabel('n=time')
legend2 = legend('xhat-x','xhatmem-x');
subplot(1,2,2), title('Error')

suptitle(['a= ' num2str(a) '; c= ' num2str(c) '; varv= ' num2str(varv) '; varw= ' num2str(varw)])

sum(abs(xhat-x))/n
sum(abs(xhatmem-x))/n
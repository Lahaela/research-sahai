%Theoretical curves U/Y delay

clc
clear all
close all
set(0,'defaulttextinterpreter','latex')

% x = k, y = a
k  = [1:0.1:10]; %1 = no delay
gammap = 1-2.^(-k+1); %probability packet goes through
mu_c = 1;
varc = 1;

temp = (mu_c^2+varc)./((1-gammap).*mu_c^2+varc);
a = temp.^(1./(2.*k));

plot(k,a,'m','Linewidth',2)
xlabel('Delay K'), ylabel('Threshold for A')
title('Theoretical Effect of Delay on Tolerable A, Drop = $$2^{-k}$$, C$$\sim$$N(1,1)')

set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
print('-dpdf','-r100','theoretical_xkya_p2-k.pdf');
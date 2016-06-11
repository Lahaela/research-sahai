%Theoretical Curves A, D, R, Beta

clc
%clear all
%close all
figure()

set(0,'defaulttextinterpreter','latex');

D = 0:1:1000;
%D = 0:1:10;
b = 1; %rate of information
code = (D+1)*b; %number of bits in code
R = 1/3*(D+1)*b; %number of bits of information

mu_c = 1;
%varc = 1;
varc = 2.^(-2*R);

dropp = 2.^(-(code - R +1));
%test = [D; dropp]
%dropp = 2.^(-(k));

% gammap = 0:0.01:1;
% dropp = 1-gammap;

z = (mu_c^2+varc)./(dropp.*mu_c^2+varc);
a = z.^(1./(2*(D+1)));
% k = log(z)./(2*log(a));

plot(D,a,'m','LineWidth',2)
xlabel('Delay d','FontSize',24), ylabel('Tolerable growth a','FontSize',24) %font size 24
%title('Theoretical Effect of Delay on A, Drop = $$2^{-(1*D - 0.5*msg)}$$')

% set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]); %also do 6 x 6
% print('-dpdf','-r100','150307_forgireeja_stem.pdf');
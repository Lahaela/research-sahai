%Theoretical curve of alpha vs A

%clc
hold all
%clear all
%close all
figure()

set(0,'defaulttextinterpreter','latex');

D = 50;
b = 1; %bits per unit time across channel
code = (D+1)*b; %number of bits in code
alph = 0.01:0.01:1; %rate n/(n+k)
R = alph*(D+1)*b; %number of bits of information

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

%plot(alph,a,'LineWidth',2)
plot(R,a,'LineWidth',2)
xlabel('R = number of info bits','FontSize',24), ylabel('Tolerable growth A','FontSize',24) %font size 24
%ylim([1 1.5])
title('Theoretical Effect of Alpha on A, Various Delays', 'FontSize', 20)

% set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
% print('-dpdf','-r100','150424_alphvsa_all.pdf');

%%
%Plot of max
%Strategy is: try all alphas available, pick best

%clc
%clear all
%close all
%figure()
%figure()
hold all

set(0,'defaulttextinterpreter','latex');

%D independent values
%alph = 0.01:0.01:1;
b = 1; % bits per unit time

mu_c = 1;
size_R = [];
track_varc = [];
track_dropp = [];
track_a = [];

for D = 0:49
    code = floor((D+1)*b); %number of bits in codeword
    if code==0
        alph = [0 0];
    else
        alph = 0:1/code:1; %rate n/(n+k)
    end
    alph = alph(2:end);
    R = alph*(D+1)*b; %number of bits of information
    %varc = 2.^(-2*R);
%     for c=1:(D+1)
%         if varc(c) < 2^(-32)
%             varc(c) = 2^(32);
%         end
%     end
    %varc = 0.75.^(2*R);
    %varc = 1./R;
    varc = 3;

    dropp = 2.^(-(code - R + 1));
    %dropp = (0.5).^(code - R + 1);
    %dropp = 1./R;
    %dropp = 0.5;

    z = (mu_c^2+varc)./(dropp.*mu_c^2+varc);
    a = z.^(1./(2*(D+1)));
%     for a_idx=1:length(a)
%         if (mod(a_idx, length(alph)/(D+1)) ~= 0)
%             a(a_idx) = 0;
%         end
%     end
    [max_val, max_idx] = max(a);
    track_a = [track_a max_val];
    max_R = max_idx/length(alph)*code;
    max_alph = max_idx/length(alph);
    size_R = [size_R max_R];
    track_varc = [track_varc (2^(-2*max_R))];
    %track_varc = [track_varc 1/max_R];
    track_dropp = [track_dropp (2^(-(code - max_R +1)))];
    %scatter(D, max_val, 'b','filled','LineWidth',2);
%     if (D < 20)
%         scatter(max_alph, max_val, 'b','filled','LineWidth',2);
%     elseif (D < 40)
%         scatter(max_alph, max_val, 'r','filled','LineWidth',2);
%     elseif (D < 60)
%         scatter(max_alph, max_val, 'g','filled','LineWidth',2);
%     elseif (D < 80)
%         scatter(max_alph, max_val, 'm','filled','LineWidth',2);
%     else
%         scatter(max_alph, max_val, 'y','filled','LineWidth',2);
%     end
end

%plot(alph,a,'LineWidth',2)
% xlabel('Alpha = Rate of Info','FontSize',24), ylabel('Tolerable growth A','FontSize',24) %font size 24
% ylim([1 1.5])
% title('Theoretical (Max) Alpha vs A, D = 0 $$\rightarrow$$ 99', 'FontSize', 20)

% figure()
% plot(0:D, size_R,'LineWidth',2)
% xlabel('Delay','FontSize',24), ylabel('R = number of info bits','FontSize',24)
% title('Strategy: Best Alpha (D+1) Choices','FontSize',20)
% 
% figure()
% plot(0:D, track_varc,'LineWidth',2)
% xlabel('Delay','FontSize',24), ylabel('Quantitative Noise','FontSize',24)
% title('Strategy: Best Alpha (D+1) Choices','FontSize',20)
% 
% figure()
% plot(0:D, track_dropp,'LineWidth',2)
% xlabel('Delay','FontSize',24), ylabel('Drop Probability','FontSize',24)
% title('Strategy: Best Alpha (D+1) Choices','FontSize',20)

%figure()
plot(0:D, track_a,'LineWidth',2)
xlabel('Delay','FontSize',24), ylabel('Tolerable Growth a','FontSize',24)
title('Strategy: Optimal R (Encoding)','FontSize',20)

% set(gcf,'PaperUnits','inches','PaperSize',[6,6],'PaperPosition',[0 0 6 6]);
% print('-dpdf','-r100','150430_Dvsa.pdf');
% 
% x = 0:D;
% hold all
% plot(x, 1./(x+1),'LineWidth',2)
% plot(x, exp(-(x+1)),'LineWidth',2)
% plot(x,2.^(-(x+1)),'LineWidth',2)
% legend('System','1/(D+1)','e^{-(D+1)}','2^{-(D+1)}')

%legend('B = 0.1', 'B = 0.2', 'B = 0.3', 'B = 0.4', 'B = 0.5', 'B = 0.6', 'B = 0.7', 'B = 0.8', 'B = 0.9', 'B = 1', 'Location','Best')
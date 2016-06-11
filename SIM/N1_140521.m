clc;
clear all;


%% 2 Player Coin Bids

A = 10; % reward
C2 = 5;
trials = 1000;

coldisp = [1 3 5 7 9 2 4 6 8 10];

for C1=1:10 %lower value = stronger player
    out = zeros(1,trials);
    for coin=1:trials
        bid1 = A/C1*rand(1); %random bid interval [0,A/C1]
        bid2 = A/C2*rand(1);
        out(1,coin) = (bid1 > bid2);
    end
    heads = sum(out==1); %number of times player1 wins
    subplot(5,2,coldisp(1,C1)), bar(1:2,[heads, trials-heads],1)
    subplot(5,2,coldisp(1,C1)), set(gca, 'xticklabel', ['Player 1'; 'Player 2'])
    subplot(5,2,coldisp(1,C1)), ylim([0, 1000])
    title(subplot(5,2,coldisp(1,C1)),['Cost1 = ' num2str(C1) ' Maxbid1 = ' num2str(A/C1)])
end

suptitle('2 Player Coin Bids')

%% 3 Player Coin Bids

trials = 1000;

A = 10; % reward
C1 = 1; %strong
C2 = 5;
C3 = 10;

out = zeros(1,trials);
for coin=1:trials
    bids = [A/C1*rand(1), A/C2*rand(1), A/C3*rand(1)];
    out(1,coin) = find(max(bids)==bids);
end
numwins = [sum(out==1), sum(out==2), sum(out==3)];
bar(1:3,numwins,1)
set(gca, 'xticklabel', ['Player 1'; 'Player 2'; 'Player 3'])
ylim([0, 1000])
title(['Maxbid1 = ' num2str(A/C1) '; Maxbid2 = ' num2str(A/C2) '; Maxbid3 = ' num2str(A/C3)])
suptitle('3 Player Coin Bids')

%% Multiplayer Coin Bids

trials = 1000;

A = 10;
C = [1 5 10]; %position indicates player

out = zeros(1,trials);
for coin=1:trials
    bids = zeros(1,length(C));
    for i=1:length(C)
        bids(1,i) = A/C(1,i)*rand(1);
        out(1,coin) = find(max(bids)==bids);
    end
end
numwins = zeros(1,length(C));
for i=1:length(C)
    numwins(1,i) = sum(out==i);
end

bar(1:length(C),numwins,1)
ylim([0, 1000])
suptitle('Multiplayer Coin Bids')

xlabels = cell(1,length(C));
for i=1:length(C)
    xlabels(1,i) = {strcat('Player ',num2str(i))};
end

% xls = zeros(1,length(C));
% for i=1:length(C)
%     xls = [xls; 'Player ' num2str(i)];
% end

set(gca,'xticklabel',xlabels)


%% Nash Equilibrium

trials = 1000;

A = 10;
C = [1 5 10]; %lower number = stronger player, position indicates player identity

[c, ix] = sort(C,2);
c = c(:,1:2);
ix = ix(:,1:2); % to preserve identity

out = zeros(1,trials);
for coin=1:trials
    bids = [A/c(1,1)*rand(1) A/c(1,2)*rand(1)];
    out(1,coin) = (bids(1,1)>bids(1,2));
end

numwins = sum(out==1);
bar(1:2,[numwins, trials-numwins],1)
ylim([0 1000])
title('Multiplayer with Nash Equilibrium')
set(gca,'xticklabel',['Player ' num2str(ix(1,1)); 'Player ' num2str(ix(1,2))])

%% Theorem 1

trials = 1000;

A = 10;
C = [1 2 10];

[c, ix] = sort(C,2);
c = c(:,1:2);
ix = ix(:,1:2); % to preserve identity

bid1 = A/c(1,2).*rand(trials,1);
p = (c(1,2)-c(1,1))/c(1,2);
bid2 = zeros(trials,1);
for i=1:trials
    if rand() <= p
        bid2(i,1) = 0;
    else
        bid2(i,1) = A/c(1,2)*rand(1);
    end
end

out = bid1>bid2;
numwins = sum(out==1);
bar(1:2,[numwins, trials-numwins],1)
ylim([0 1000])
set(gca,'xticklabel',['Player ' num2str(ix(1,1)); 'Player ' num2str(ix(1,2))])
title('Theorem 1 simulator')
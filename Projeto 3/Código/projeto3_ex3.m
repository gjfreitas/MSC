%% Task 3.1
clear all
close all
clc

rng(2023) % fix the seed to obtain always the same results

N = 1000;
mu = [1.6, 2, 2.6];

l_max = 1000;
l=zeros(length(0:N),2);

r = zeros(length(0:N),2,length(mu)); %initialization the particle's positions 


for k=1:length(mu)
    x=0;
    y=0;
    for i=2:length(0:N)    % i is the number of jumps from 1 to N
        r_rand = rand(1);
        ang_rand = rand(1);

        l(i,k) = ( 1 - r_rand*(1- l_max^(1-mu(k))) )^(1/(1-mu(k)));
        angle = 2*pi*ang_rand;

        mov_x = l(i,k) * cos(angle);
        mov_y = l(i,k) * sin(angle);

        x = mov_x + x; 
        y = mov_y + y;

        r(i,:,k) = [x y]; % new position of the particle

    end

end


% Isotropic 2D-random walks having a fixed length of jumps, l = 1

x = 0;
y = 0;

trajectories=zeros(length(0:N),2);

step_size = 1;

for i = 2:length(0:N)
    % Generate random direction
    direction = randi(4); % generates a random number [1,2,3,4] with equal probability
    switch direction
        case 1 % Move up
            y = y + step_size;
        case 2 % Move right
            x = x + step_size;
        case 3 % Move down
            y = y - step_size;
        case 4 % Move left
            x = x - step_size;
    end
    trajectories(i,:) = [x y];
end


save('Task3_1.mat')

%% Plots

clear all
close all
clc


if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end

load Task3_1.mat

col=['r' 'b' 'k', 'g'];

figure(1);
% hold on
for i = 1:length(mu)
    plot(r(:,1,i), r(:,2,i),[col(i) '-'],'LineWidth',1.5);
    hold on
end
plot(trajectories(:,1), trajectories(:,2),[col(end) '-'],'LineWidth',1.5)
xlabel('X');
ylabel('Y');
legend('\mu = 1.6', '\mu = 2', '\mu = 2.6', 'Fixed (l = 1)', 'Location','northwest');
grid on
hold off;
filename = 'Imagens/Task3_1_all_XY.eps';
print('-depsc', filename)



figure(2);
% hold on
for i = 1:length(mu)
    plot(0:N, r(:,1,i),[col(i) '-.'],'LineWidth',1.5);
    hold on
end
plot(0:N, trajectories(:,1),[col(end) '-'],'LineWidth',1.5)
xlabel('t');
ylabel('X');
legend('\mu = 1.6', '\mu = 2', '\mu = 2.6', 'Fixed (l = 1)');
grid on
hold off;
filename = 'Imagens/Task3_1_all_Y0.eps';
print('-depsc', filename)


figure(3);
% hold on
for i = 1:length(mu)
    plot(0:N, r(:,2,i),[col(i) '-.'],'LineWidth',1.5);
    hold on
end
plot(0:N, trajectories(:,2),[col(end) '-'],'LineWidth',1.5)
xlabel('t');
ylabel('Y');
legend('\mu = 1.6', '\mu = 2', '\mu = 2.6', 'Fixed (l = 1)', 'Location','northwest');
grid on
hold off;
filename = 'Imagens/Task3_1_all_X0.eps';
print('-depsc', filename)



figure(4);
% hold on
for i = 1:length(mu)
    plot3(0:N, r(:,1,i), r(:,2,i),[col(i) '-'],'LineWidth',1.5);
    hold on
end
plot3(0:N, trajectories(:,1),trajectories(:,2),[col(end) '-.'],'LineWidth',1.5)
xlabel('t');
ylabel('X');
zlabel('Y');
legend('\mu = 1.6', '\mu = 2', '\mu = 2.6', 'Fixed (l = 1)');
grid on
hold off;
filename = 'Imagens/Task3_1_all_3D.eps';
print('-depsc', filename)


figure(5)
plot(trajectories(:,1),trajectories(:,2),[col(end) '-'],'LineWidth',1.5)
title('Fixed (L = 1)')
xlabel('X')
ylabel('Y')
grid on
axis square
filename = 'Imagens/Task3_1_fig1.eps';
print('-depsc', filename)


figure(6)
plot(r(:,1,1), r(:,2,1),[col(1) '-'],'LineWidth',1.5);
title('\mu = 1.6')
xlabel('X')
ylabel('Y')
grid on
axis square
filename = 'Imagens/Task3_1_fig2.eps';
print('-depsc', filename)

figure(7)
plot(r(:,1,2), r(:,2,1),[col(2) '-'],'LineWidth',1.5);
title('\mu = 2')
xlabel('X')
ylabel('Y')
grid on
axis square
filename = 'Imagens/Task3_1_fig3.eps';
print('-depsc', filename)

figure(8)
plot(r(:,1,3), r(:,2,1),[col(3) '-'],'LineWidth',1.5);
title('\mu = 2.6')
xlabel('X')
ylabel('Y')
grid on
axis square
filename = 'Imagens/Task3_1_fig4.eps';
print('-depsc', filename)



%% Task 3.2
clear all
close all
clc

if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end


mu = [1.6 2 2.6];
l_max = 1000;
N = 1E7;

l = zeros(N+1,length(mu));
for i=1:length(mu)
    for j=1:N
        r = rand();
        l(j,i) = (1-r*(1-l_max^(1-mu(i))))^(1/(1-mu(i))); % fórmula que se quer estudar
    end
end

l_theory = linspace(1,l_max,N+1);
C = zeros(1,length(mu));
P = zeros(N+1,length(mu));
for i=1:length(mu)
    C(i) = (mu(i)-1)/(1-l_max^(1-mu(i))); % Constante de normalização
    P(:,i) = C(i).*l_theory.^(-mu(i));   % Levy Flights distribution

    figure
%     sgtitle(['\mu=', num2str(mu(i)), newline,'l_{max}=',num2str(l_max)])
    
    h = histogram(l(:,i),'BinWidth',1,'Normalization','probability');
    hold on
    plot(l_theory,P(:,i),'r.')
    xlabel('$l$', 'Interpreter', 'latex'), ylabel('P')
    xlim([0 10]),ylim([0 1])
    legend('P_{experimental}','P_{Lévy}')
    
    set(gcf,'position',[10,10,800,400])

    if i == 1
        filename = 'Imagens/Task3_2_1.eps';
        print('-depsc', filename)
    elseif i == 2
        filename = 'Imagens/Task3_2_2.eps';
        print('-depsc', filename)
    else
        filename = 'Imagens/Task3_2_3.eps';
        print('-depsc', filename)
    end
end








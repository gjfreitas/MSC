%% Task 2.1
clear all
close all
clc

tic

N = 50000;
xb = -30;
dt = 10;
D = 1/4;

step_size = 1;

t = 50000;

Fc = zeros(length(0:t),1);    % First passage count
Sc = Fc;    % Survival count


for i = 1:N
    % Start at origin
    x = 0;
    y = 0;
    
    % Take random steps
    for j = 2:length(0:t)
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

        if x==xb
            Fc(j) = Fc(j)+1;
            break
        else
            Sc(j) = Sc(j)+1;
        end
    end
end


idx=1;
for i= 0:dt:length(0:t)-dt
    soma = 0;
    soma2 = 0;
    for j=1:dt
        soma = soma + Fc(i+j);
        soma2 = soma2 + Sc(i+j);
    end
    N_f(idx) = soma;
    N_s(idx) = soma2;
    idx=idx+1;
end

F = N_f./(N*dt);
S = (N_s./(N*dt));

t_dt = 1:dt:t-2;

F_theor = abs(xb-0)./(sqrt(4*pi*D*t_dt.^3)) .* exp(-((xb-0)^2)./(4*D.*t_dt));
S_theor = erf(abs(xb-0)./(2*sqrt(D.*t_dt)));

toc

figure
subplot(1,2,1)
plot(t_dt, F,'.')
hold on
grid on
plot(t_dt, F_theor, 'k-','LineWidth', 1.5)
title('First passage Probability, F(t)')
xlabel('t')
ylabel('F')
legend('F', 'F_{theoretical}')

subplot(1,2,2)
plot(log10(t_dt), log10(F),'.')
hold on
grid on
plot(log10(t_dt), log10(F_theor),'k-','LineWidth',1.5)
title('First passage Probability, F(t) (log-log)')
xlabel('log_{10}t')
ylabel('log_{10}F')
legend('log_{10}(F)', 'log_{10}(F_{theoretical})')

figure()
plot(t_dt(2:end), S(2:end),'r--','LineWidth',3)
hold on
grid on
plot(t_dt, S_theor,'k-','LineWidth',1.5)
title('Survival Probability, S(t)')
legend('S','S_{theoretical}')


save('Task2.mat')

%% Plots
clear all
close all
clc

if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end

load Task2.mat

figure(1)
plot(t_dt, F,'.')
hold on
grid on
plot(t_dt, F_theor, 'k-','LineWidth', 1.5)
% title('First passage time Probability, F(t)')
xlabel('t')
ylabel('F(t)')
legend('F', 'F_{teórico}')
filename = 'Imagens/Task2_fig1.eps';
print('-depsc', filename)

figure(2)
plot(log10(t_dt), log10(F),'.')
hold on
grid on
plot(log10(t_dt), log10(F_theor),'k-','LineWidth',1.5)
% title('First passage time Probability, F(t)')
xlabel('log_{10}t')
ylabel('log_{10}F(t)')
legend('log_{10}(F)', 'log_{10}(F_{teórico})', 'Location','southeast')
filename = 'Imagens/Task2_fig2.eps';
print('-depsc', filename)

figure(3)
plot(t_dt(2:end), S(2:end),'r--','LineWidth',3)
hold on
grid on
plot(t_dt, S_theor,'k-','LineWidth',1.5)
% title('Survival Probability, S(t)')
xlabel('t')
ylabel('S(t)')
legend('S','S_{teórico}')
filename = 'Imagens/Task2_fig3.eps';
print('-depsc', filename)



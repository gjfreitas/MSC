%% Task 1.1
clear all
close all
clc

rng(0) % fix the seed to obtain always the same results

% Parameters
num_steps = 100; % Number of steps
step_size = 1; % Step size
num_trajectories = 3; % Number of trajectories to simulate

% Initialize trajectories
trajectories = zeros(num_trajectories, num_steps+1, 2);

% Simulate trajectories
for i = 1:num_trajectories
    % Start at origin
    x = 0;
    y = 0;
    trajectories(i,1,:) = [x y];
    
    % Take random steps
    for j = 2:num_steps+1
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
        trajectories(i,j,:) = [x y];
    end
end

save('Task1_1.mat')


%% Plots - Task 1.1
clear all
close all
clc

load Task1_1.mat

if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end

col=['k','b','r'];

figure(1);
hold on;
for i = 1:num_trajectories
    plot(0:num_steps, trajectories(i,:,1),[col(i) '-.'], 'MarkerSize', 2,'LineWidth',1.5);
end
xlabel('t');
ylabel('X');
legend('Trajetória 1', 'Trajetória 2', 'Trajetória 3', 'Location','northwest');
grid on
hold off;
filename = 'Imagens/Task1_1_fig1.eps';
print('-depsc', filename)

figure(2);
hold on;
for i = 1:num_trajectories
    plot(0:num_steps, trajectories(i,:,2),[col(i) '-.'], 'MarkerSize', 2,'LineWidth',1.5);
end
xlabel('t');
ylabel('Y');
legend('Trajetória 1', 'Trajetória 2', 'Trajetória 3', 'Location','northwest');
grid on
hold off;
filename = 'Imagens/Task1_1_fig2.eps';
print('-depsc', filename)


figure(3);
% hold on
for i = 1:num_trajectories
    plot3(0:num_steps, trajectories(i,:,1), trajectories(i,:,2),[col(i) '-o']);
    hold on
end
xlabel('t');
ylabel('X');
zlabel('Y');
legend('Trajetória 1', 'Trajetória 2', 'Trajetória 3');
grid on
hold off;
filename = 'Imagens/Task1_1_fig3.eps';
print('-depsc', filename)

figure(4);
% hold on
for i = 1:num_trajectories
    plot(trajectories(i,:,1), trajectories(i,:,2),[col(i) '-o']);
    hold on
end
xlabel('X');
ylabel('Y');
legend('Trajetória 1', 'Trajetória 2', 'Trajetória 3');
grid on
hold off;
filename = 'Imagens/Task1_1_fig4.eps';
print('-depsc', filename)

%% Task 1.2
clear all
close all
clc


t_even = 50000;
t_odd = 50001;

N_vals = [1E1, 1E2, 1E3, 1E4, 1E5, 1E6];

for n=1:length(N_vals)
    N = N_vals(n);
    disp(['N -> ',num2str(N)])
    % N = 1E5;

    end_point_even = nan(N,2);
    end_point_odd  = nan(N,2);
    
    end_points = nan(N*2,2);
    
    tic
    parfor rep=1:N    
        end_point_even(rep,:) = random_walk_2d(t_even);
        end_point_odd(rep,:)  = random_walk_2d(t_odd);
    end
    
    end_point = [end_point_even;end_point_odd];
    [counts, number] = groupcounts(end_point);
    
    % Calculate the probability at each site
    counts = counts/N;
    
    % Repeat for even and odd times, and average:
    x = number{1};
    y = number{2};
    Pmean = 0.5*counts;
    toc
    
    % Check the normalization, averages and normalization
    disp(['Norm = ', num2str(sum(Pmean)), '; Expected = 1'])
    
    avg_x = sum(x.*Pmean);
    disp(['Avg(x) = ', num2str(avg_x), '; Expected = 0'])
    
    avg_y = sum(y.*Pmean);
    disp(['Avg(y) = ', num2str(avg_x), '; Expected = 0'])
    
    sd = sum((x.^2+y.^2).*Pmean);
    disp(['Variance = ', num2str(sd), '; Expected = ', num2str((t_odd + t_even)/2)]);

    fprintf('\n')
end

% Calculate the theoretical distribution function
t_mean = (t_odd + t_even)/2;
lim = 2*sqrt(t_odd);
[X,Y] = meshgrid(-lim:lim);
const = 1/(pi*t_mean);
Ptheor = const*exp(-(X.^2 +Y.^2)/t_mean);

save('task1_2.mat')

%% Plots - Task 1.2
clear all
close all
clc

if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end

load task1_2.mat

figure(1)
plot3(x,y,Pmean,'bo')
xlabel('x'),ylabel('y'),zlabel('$\hat{P}$', 'Interpreter', 'latex')
grid on

filename = 'Imagens/Task1_2_fig1.eps';
print('-depsc', filename)

figure(2)
mesh(X,Y,Ptheor)
xlabel('x'),ylabel('y'),zlabel('P_{teórico}')
grid on
filename = 'Imagens/Task1_2_fig2.eps';
print('-depsc', filename)

figure(3)
plot3(x,y,Pmean,'bo'),view(0,0)
xlabel('x'),ylabel('y'),zlabel('$\hat{P}$', 'Interpreter', 'latex')
grid on
filename = 'Imagens/Task1_2_fig3.eps';
print('-depsc', filename)

figure(4)
mesh(X,Y,Ptheor),view(0,0)
xlabel('x'),ylabel('y'),zlabel('P_{teórico}')
grid on
filename = 'Imagens/Task1_2_fig4.eps';
print('-depsc', filename)




function endpoint = random_walk_2d(num_steps)
    step_size = 1;

    % Start at origin
    x = 0;
    y = 0;
    % Take random steps
    for j = 2:num_steps+1
        % Generate random direction
        direction = randi([1 4]); % generates a random number [1,2,3,4] with equal probability
        if direction == 1 % Move up
                y = y + step_size;
        elseif direction == 2 % Move right
                x = x + step_size;
        elseif direction == 3 % Move down
                y = y - step_size;
        elseif direction == 4 % Move left
                x = x - step_size;
        end
    end
    endpoint = [x,y];
end

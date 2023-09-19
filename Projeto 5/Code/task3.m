%% Task 3
clear all
close all
clc

rng(0)

tic
% Parameters for simulations:
N = 10000; c = 50; m = 100;

% 8. Repeat steps 1-7 𝑚 times and average the degree distribution 𝑃(𝑞) over the 𝑚 realizations of the network
parfor n = 1:m
    rho(n) = degree_dist_function(N,c);
end


% Show that the Pearson coefficient is small.
pearson_coeff = mean(rho);

execution_time = toc;

save('Task3.mat')

disp(['Elapsed time is ', num2str(execution_time),' seconds (With Parallel Computing Toolbox).'])
disp(['Peorson coefficient = ', num2str(pearson_coeff)])

function [rho] = degree_dist_function(N,c)
    % 1. Initialize a 𝑁 × 𝑁 matrix 𝑨 with all entries set to 0.
    A = zeros(N,N);
    
    % 2. For each edge generate two integers 𝑖 and 𝑗 uniformly at random between 1 and 𝑁.
    e = 0;
    L = N/2 * c;
    while e ~= L % 5. Repeat steps 2-4 𝐿 times
        i = randi([1, N]);
        j = randi([1, N]);
        
        % 3. If 𝑖 = 𝑗 or 𝐴𝑖𝑗 = 1 repeat step 2, otherwise proceed to the next step
        if A(i, j) ~= 1 && i ~= j
            % 4. Update matrix 𝑨 with the new edge by setting 𝐴𝑖𝑗 = 𝐴𝑗𝑖 = 1.
            A(i, j) = 1;
            A(j, i) = 1;
            e = e + 1;
        end
    end
    
    % 6. Calculate degrees, 𝑞𝑖 = ∑ 𝐴𝑖𝑗
    q = zeros(N, 1);
        
    for i = 1:N
        for j = 1:N
            q(i) = q(i) + A(i,j);
        end
    end
    
    % 7. Calculate the number of vertices with degree 𝑞, i.e., 𝑁(𝑞), and then find degree distribution 𝑃(𝑞) = N(q) / N.
    N_q = zeros(N, 1); % Number of vertices with degree q
    for i = 1:N
        N_q(q(i)) = N_q(q(i)) + 1;
    end
    P = N_q / N; % Degree distribution

    q_mean = 1/N*sum(q);
    B = 1/N/q_mean.*sum(q.*(q-1));
    q_mean_2 = 1/N*sum(q.^2);
    q_mean_3 = 1/N*sum(q.^3);
    

    Q = q_mean_2/q_mean;
    sigma2 = q_mean_3/q_mean - (q_mean_2)^2/(q_mean)^2;
    
    soma = 0;
    for i=1:N
        for j=1:N
            soma = soma + A(i,j)*(q(i)-Q)*(q(j)-Q);
        end
    end
    rho = soma/(N*q_mean*sigma2);
end

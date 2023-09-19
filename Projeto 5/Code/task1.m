%% Task 1
clear all
close all
clc

rng(0)

tic
% Parameters for simulations:
N = 1E4; c = 50; m = 100;

% 8.1 Repeat steps 1-7 𝑚 times
parfor n = 1:m
    [PP, BB, QQ, Q2, Q3] = degree_dist_function(N,c);
    P_q(:,n) = PP;
    B_all(n) = BB;
    qq(:,n) = QQ;
    q2_mean(n) = Q2;
    q3_mean(n) = Q3;
end
% 8.2 Average the degree distribution 𝑃(𝑞) over the 𝑚 realizations of the network
meanP = mean(P_q,2);


% 9. Calculate the mean degree 〈𝑞〉
q_mean = mean(qq,2);

% 10. Calculate the branching coefficient
B = mean(B_all);

% 11. Calculate second and third moments, <q^2> and <q^3>
avg_q2 = mean(q2_mean);
avg_q3 = mean(q3_mean);

% Show that 𝐵/〈𝑞〉 ≈ 1.
show_that = B/q_mean;


% Compare with the theorical results

% P(q) theorical, P(q) = e^-c * c^q/q!
P = zeros(N, 1);
for i = 1:N
    P(i) = exp(-c) * c.^(i)./ (factorial(i));
end

%<q> téorico
avg_q_theor = c;

% B téorico
B_theor = c;

execution_time = toc;

save('Task1.mat')

%% Draw plots
clear all
close all
clc

if ~exist('Imagens', 'dir')
    mkdir('Imagens')
end

load Task1.mat

disp(['Elapsed time is ', num2str(execution_time),' seconds (With Parallel Computing Toolbox).'])
disp(['〈𝑞〉 got: ', num2str(q_mean), ', Expected: ', num2str(avg_q_theor)])
disp(['𝐵 got: ', num2str(B), ', Expected: ', num2str(B_theor)])
disp(['𝐵/〈𝑞〉 expected ≈ 1. Got: ', num2str(show_that)])
disp(['<q^2> obtained: ', num2str(avg_q2)])
disp(['<q^3> obtained: ', num2str(avg_q3)])

figure(1)
plot(1:N,meanP,'b-','LineWidth', 1.5)
hold on
plot(1:N, P, 'ro','LineWidth', 1.5);
xlabel('q');
ylabel('P(q)');
% title('Degree Distribution of Random Graph');
axis([0 100 0 0.06])
legend('Numerical', 'Theoretical')

filename = 'Imagens/Task1.eps';
print('-depsc', filename)


%% Function

function [P, B, q_mean, q_mean_2, q_mean_3] = degree_dist_function(N,c)
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

    q_mean = 1/N.*sum(q);
    B = 1/N/q_mean.*sum(q.*(q-1));
    q_mean_2 = 1/N*sum(q.^2);
    q_mean_3 = 1/N*sum(q.^3);
end

%% Changing the parameter N
clear all
close all
clc

rng(0)

tic
% Parameters for simulations:
c = 50; m = 100;
n_mat = [1E2, 1E3, 1E4];
avgP = nan(n_mat(end), length(n_mat));
avg_var = nan(3,1);
avg_mean = nan(3,1);

for nn = 1:length(n_mat)
    clear P_q B_all qq q2_mean q3_mean
    N = n_mat(nn);
    parfor n = 1:m
        [PP, BB, QQ, Q2, Q3, V1] = degree_dist_function(N,c);
        P_q(:,n) = PP;
        B_all(n) = BB;
        qq(:,n) = QQ;
        q2_mean(n) = Q2;
        q3_mean(n) = Q3;
        var(n) = V1;
    end
    avgP(1:N,nn) = mean(P_q,2);
    avg_var(nn) = mean(var);
    avg_mean(nn) = mean(qq,2);
end

% Compare with the theorical results
% P(q) theorical, P(q) = e^-c * c^q/q!
for i = 1:N
    P(i) = exp(-c) * c.^(i)./ (factorial(i));
end

col=['k','g','b'];
close all
figure(1);
hold on
for i = 1:length(n_mat)
    plot(1:n_mat(i), avgP(1:n_mat(i),i),[col(i) '-'],'LineWidth',1.5);
end
plot(1:N, P,'ro','LineWidth',1.5);
axis([0 100 0 0.1])
xlabel('q');
ylabel('P(q)');
legend('Numerical (N = 10^2)','Numerical (N = 10^3)', 'Numerical (N = 10^4)','Theorical', 'Location','northwest')

axes('Position',[.7 .7 .2 .2])
box on
hold on
for i = 2:length(n_mat)
    plot(1:n_mat(i), avgP(1:n_mat(i),i),[col(i) '-'],'LineWidth',1.5);
end
plot(1:N, P,'ro','LineWidth',1.5);
axis([45 55 0.045 0.06])

% Check that when N-> inf, the distribution approaches a Poisson, using the var/mean result
check_poisson = avg_var./avg_mean;

save('Task1_extra_N.mat')

%% Changing the parameter c
clear all
close all
clc

rng(0)

tic
% Parameters for simulations:
N = 1E4; m = 100;
c_mat = [25, 50, 75];
avgP = nan(c_mat(end), length(c_mat));
P = nan(c_mat(end), length(c_mat));

for nn = 1:length(c_mat)
    clear P_q B_all qq q2_mean q3_mean
    c = c_mat(nn);
    parfor n = 1:m
        [PP, BB, QQ, Q2, Q3] = degree_dist_function(N,c);
        P_q(:,n) = PP;
        B_all(n) = BB;
        qq(:,n) = QQ;
        q2_mean(n) = Q2;
        q3_mean(n) = Q3;
    end
    avgP(1:N,nn) = mean(P_q,2);

    % Compare with the theorical results
    % P(q) theorical, P(q) = e^-c * c^q/q!
    for i = 1:N
        P(i,nn) = exp(-c) * c.^(i)./ (factorial(i));
    end
end


col=['k','g','b'];
close all
figure(1);
hold on
for i = 1:length(c_mat)
    plot(1:N, avgP(:,i),[col(i) '-'],'LineWidth',1.5);
%     plot(1:N, P(:,i),'ro','LineWidth',1.5);
end
axis([0 100 0 0.1])
xlabel('q');
ylabel('P(q)');
legend('Numerical (c = 25)','Numerical (c = 50)', 'Numerical (c = 75)')


save('Task1_extra_C.mat')

%%
function [P, B, q_mean, q_mean_2, q_mean_3, var] = degree_dist_function(N,c)
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

    var_sum = 0;
    for i=1:N
        var_sum = var_sum + (q(i)-q_mean)^2;
    end
    var = var_sum/(N-1);

end

%% Task 2
clear all
close all
clc


rng(0)

tic
% Parameters for simulations:
N = 1E4; c = 50; m = 100;

parfor nn = 1:m
    [CC, PP, NPT, NTR] = t2_degree_dist_function(N,c)
    P_q(:,nn) = PP;
    C_all(nn) = CC;
    n_pt(nn) = NPT;
    n_tr(nn) = NTR;
end

meanP = mean(P_q,2);

C = mean(C_all);
C_theor = c/N;

mean_ntr = mean(n_tr);
mean_npt = mean(n_pt);

execution_time = toc;

save('Task2.mat')

%% Plot
clear all
close all
clc

% if ~exist('Imagens', 'dir')
%     mkdir('Imagens')
% end

load Task2.mat

disp(['Elapsed time is ', num2str(execution_time),' seconds (With Parallel Computing Toolbox).'])
disp(['C expected ≈ 〈𝑞〉/N = ', num2str(c/N),'. Got: ', num2str(C)])

% figure(2)
% q=1:100;
% P_th=exp(-c)*(c.^q)./factorial(q);
% bar(q,meanP(q))
% hold on 
% plot(P_th,'LineWidth',2)
% set(gcf,'color','w');
% xlabel('q'); ylabel('P(q)');
% legend('Numerical','Theoretical');
% 
% filename = 'Imagens/Task2.eps';
% print('-depsc', filename)


function [C, P, n_pt, n_tr] = t2_degree_dist_function(N,c)
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
    
    figure(1)
    colormap('hot')
    imagesc(A)
    set(gcf,'color','w');

  
    % 6. Calculate degrees, 𝑞𝑖 = ∑ 𝐴𝑖𝑗
    q = sum(A,1);

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
    
    qq1 = 0;
    for i = 1:N-1
        qq1 = qq1 + i*(i-1)*P(i);
    end

        
    % n_pt is number of possible triangles
    n_pt = 1/6 * N * qq1;
    
    % n_tr is the number of triangles in the network
    n_tr=0;
    for i=1:N
        for j=1:N
            if A(i,j)~=0
                for k=1:N
                    n_tr=n_tr+(A(i,j)*A(j,k)*A(k,i));
                end
            end
        end
    end

    n_tr=n_tr/6;
    C=n_tr/n_pt;

end


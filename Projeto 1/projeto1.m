%% Task 1.2

clear all
close all
clc


M = 100;

avg_exact = (M+1)/2;
sd_exact = (M^2-1)/12;

%N_mat = [1E1, 1E2, 1E3, 1E4, 1E5, 1E6];
N_mat = [1E2, 1E4, 1E6];

avg = zeros(length(N_mat),1); sd = avg; err_avg = avg; err_sd = avg;

for idx = 1:length(N_mat)
    N = N_mat(idx);
    x=randi(M,[1 N]); %generating N random numbers with unif distr.
    avg(idx) = sum(x)/N;
    sd(idx) = sum((x-avg(idx)).^2)/N;

    err_avg(idx) = abs(avg(idx)-avg_exact)/avg_exact*100;
    err_sd(idx) = abs(sd(idx)-sd_exact)/sd_exact*100;

end


figure(1)
plot(N_mat,err_avg,'r.-')
xlabel('N'); ylabel('error(%) mean')

figure(2)
plot(N_mat,err_sd,'r.-')
xlabel('N'); ylabel('error(%) sd')

%% Task 1.3

clear all
close all
clc

%N_mat = [1E1, 1E2, 1E3, 1E4, 1E5, 1E6];
N_mat = [1E2, 1E4, 1E6];

avg_z1 = zeros(3,1); avg_z2 = avg_z1; err_avg = avg_z1;

for idx = 1:length(N_mat)
    N = N_mat(idx);
    x = zeros(N,1); y = x; z = x;

    for i=1:N
      x(i) = rand(1);
      y(i) = rand(1);
      z(i) = x(i)*y(i);
    end
  avg_x = mean(x);
  avg_y =mean(y);
  avg_z1(idx) = mean(z); % <z> = 1/N sum(z)
  avg_z2(idx) = avg_x*avg_y; % <z> = <x><y>

  err_avg(idx)=(abs(avg_z1(idx) - avg_z2(idx))/ avg_z1(idx))*100;
end

figure(1)
plot(N_mat,err_avg,'r.-')
xlabel('N'); ylabel('error(%) avg')


%% Task 2.1

clear all
close all
clc

N = 1E6;

for i=1:N
    x(i) = rand(1);
end


figure('units','normalized','outerposition',[0 0 0.7 0.6])
subplot(121)
hist(x,30),title('(a) - Histograma 10⁵ valores rand(1)')
ylabel('Frequency'),xlabel('x'),ylim([0 4000])

dx = 0.005;
[M,x_m,bins] = histcounts(x,1/dx);

P = M./(N*dx);
normalize = sum(P*dx); %checking normalization (=1)

subplot(122)
plot(P,'k.-'),ylim([0 1.2])
ylabel("P(x)"),xlabel("x"),title(['(b) - PDF of x (dx = ', num2str(dx),')'])

%% Task 2.2
clear all
close all
clc

N_mat = [1E2, 1E4, 1E6];

for n = 1:length(N_mat)
    N = N_mat(n);

    for i=1:N
        x(i) = rand(1);
    end
    
    avg = mean(x);
    sd = var(x);
    
    avg_exact(n) = 1/2;
    sd_exact(n) = 1/12;
    
    err_avg(n) = abs(avg - avg_exact)/avg_exact * 100; 
    err_sd(n) = abs(sd - sd_exact)/sd_exact * 100;
end

figure(1)
plot(N_mat,err_avg,'r.-')
xlabel('N'); ylabel('error(%) avg')

%% Task 2.3
clear all
close all
clc

N_mat = [1E2, 1E4, 1E6];
dx = 0.005;

t_mean_x = zeros(length(N_mat),1); t_var_x = t_mean_x; % theoretical results

mean_x = zeros(length(N_mat),1); std_x = zeros(length(N_mat),1); % numerical results

%anallytical results:
for n=1:length(N_mat)
    N = N_mat(n);
    x = sqrt(rand(N,1));  

    [M,X_m,bins] = histcounts(x,(1/dx));

    P_n(:,n) = M./(N*dx);
    for k = 0:1/dx-1
        x_k(k+1) = (k+0.5)*dx;  %x axis
    end
    
    % theoretical results
    g = 2*sqrt(x_k);
    t_mean_x(n) = mean(g);
    t_var_x(n) = var(g);

    normalize(n)=sum(P_n(:,n).*dx);
    
    % numerical results
    mean_x(n)=2*sum((x_k.*P_n(:,n)').*dx); %why adjust to multiplication to 2?
    std_x(n)=std(x);

    if n==1
        figure('units','normalized','outerposition',[0 0 0.4 1])
    end
    while 1
        if n==1
            subplot(3,1,1)
        elseif n==2
           subplot(3,1,2) 
        elseif n==3
           subplot(3,1,3) 
        else
            break
        end
        plot(x_k, P_n(:,n), 'k.-'); hold on
        plot(sqrt(x_k), 2*sqrt(x_k),'k-')
        title(['Probability density Function, N=',num2str(N)])
        legend('Numerical PDF','Anallytical PDF');
        xlabel('x'); ylabel('P');
        break
    end
end


figure(2)
plot([-1 7],[t_mean_x(end) t_mean_x(end)],'LineWidth',3)
hold on
plot(log10(N_mat),mean_x,'ko')
xlim([1.5 6.5]),ylim([1.1 1.4])
title('Averaged values of x')
legend('Theoretical mean','Numerical mean')
xlabel('log_{10}(N)'),ylabel('<x>')


acc_mean = abs(t_mean_x-mean_x)/t_mean_x * 100;

%% Task 3.1
clear all
close all
clc

n_mat = [10,100,1000];
N = 1E6;
dy=0.005;
cols=['r','k','b'];

avg_Y = zeros(length(n_mat),1); var_Y = avg_Y; sd_Y = avg_Y;

for i=1:length(n_mat)
    n = n_mat(i);

    x = zeros(n,1);
    Y = zeros(N,1);

    for N_idx = 1:N
        x = rand(1,n);
        Y(N_idx) = sum(x)/n;
    end

    % Find the bin placement based on the Y values
    [M,Y_m,bins] = histcounts(Y,(1/dy)-1);
    size_P = 1/dy-1;
    
    P = zeros(size_P,1);
    for k = 1:size_P
        P(k) = M(k)/N;
        y(k) = (k + 0.5)*dy;
    end
    
    %checking normalization condition
    normalize(i) = sum(P);

    if i==1
        figure('units','normalized','outerposition',[0 0 0.7 0.6])
    end

    subplot(121)
    plot(y,P,[cols(i),'-'])
    hold on
    ylabel("P"),xlabel("Y"),title("Probability Density of x"),legend(['n=' num2str(n_mat(1))],['n=' num2str(n_mat(2))],['n=' num2str(n_mat(3))])
    ylim([0 0.025])

    avg_Y(i) = mean(Y_m);
    var_Y(i) = var(Y_m);
    sd_Y(i) = (std(Y_m)^2)/n;

end

subplot(122)
plot(log10([10,100,1000]),var_Y,'ko-')
hold on, grid on
plot(log10([10,100,1000]),sd_Y,'k.-','MarkerSize',15)
legend('var(Y)','\sigma^2/n')
xlabel('log_{10}(n)')
ylabel('Y')
xlim([0.7 3]),ylim([-0.02 0.1])
title('Variance and \sigma^2/n of Y')

%% Task 4
clear all
close all
clc

M = 9; % nº de caixas
N = 21; % nº de bolas
K_mat = [1E2,1E4,1E6];

N_tr=zeros(N+1,length(K_mat));

for i = 1:length(K_mat)
    K = K_mat(i);

     n=zeros(K,1);
    for j = 1:K
        x=randi(M,[1 N]); % N valores inteiros de 1 até M
        n(j) = sum(x==3); % Ponto 2
        N_tr(n(j)+1,i)=N_tr(n(j)+1,i)+1;
    end
end

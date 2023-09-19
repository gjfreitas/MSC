%% Task 1.1
clear all
close all
clc


tf=50;
t0 = 0;
t=t0:1:tf;

p = 0.5; % q = p

pos = zeros(length(t),3);
pos(1,1:3)=0; % começa em 0

figure
for n=1:3
    for i=1:length(t)-1
        r = rand();

        if r < p
            S = -1;

        elseif r > p
            S = 1;
        end

        pos(i+1,n) = S + pos(i,n); % nova posição
    end
end

plot(t,pos,'o-')
xlim([0 50])
xlabel('t')
ylabel('x(t)')
title('1-D Random Walks (t_f=50)')
legend('trajectory 1','trajectory 2','trajectory 3', 'Location','northwest')

%% Task 1.2

clear all
close all
clc
tf_matrix = [40,41; 400,401; 4000,4001];

p = 0.5; % q = p

N=50000;
x=-500:500;

P_mean=zeros(length(x),3);

for g=1:3

    P = zeros(length(x),2);

    tf_list = tf_matrix(g,:);

    for n = 1:N
        for j = 1:length(tf_list)
            tf = tf_list(j);
            t=0:1:tf;
            pos=zeros(length(t),1);
            
            for i=1:length(t)-1    % i is the number of steps from 1 to 10
                r = rand();
        
                if r < p
                    S = -1;
        
                elseif r > p
                    S = 1;
                end
        
                pos(i+1) = S + pos(i); % nova posição

                if i == length(t)-1
                    xi = pos(i+1);
                    idx = xi + abs(min(x)) + 1;
                    P(idx,j) = P(idx,j)+1;
                end
            end
        end
    end
    
    for i=1:length(P)
        P_mean(i,g) = (P(i,1)+P(i,2))/2;
    
    end
    P_mean(:,g) = P_mean(:,g)./N;

end
figure(1)
plot(x,P_mean(:,1),'r',x,P_mean(:,2),'g',x,P_mean(:,3),'b','LineWidth',1.5)
title('<P(x,t)> of 1-D Random Walks')
legend(['t = ' num2str(tf_matrix(1,1)) ', ' num2str(tf_matrix(1,2))] ...
        ,['t = ' num2str(tf_matrix(2,1)) ', ' num2str(tf_matrix(2,2))] ...
        ,['t = ' num2str(tf_matrix(3,1)) ', ' num2str(tf_matrix(3,2))])
xlabel('x')
ylabel('Probabilities')
xlim([-200 200])
grid on


tf = (tf_matrix(3,1)+tf_matrix(3,2))/2;
P_theoretical = 1/sqrt(2*pi*tf).*exp(-(x.^2)/(2.*tf));

figure(2)
grid on, hold on
xlim([-200,200])
plot(x, P_mean(:,3), x, P_theoretical,'-','LineWidth',1.5)
title('Theoretical P(x,t) and <P(x,t)> (t=4000,4001)')
legend('<P(x,t)>','P_{t}')
xlabel('x')
ylabel('Probabilities')

err1 = immse(P_theoretical',P_mean(:,1)); 
err2 = immse(P_theoretical',P_mean(:,2)); 
err3 = immse(P_theoretical',P_mean(:,3));

disp(['Erro 1:',num2str(err1)]);
disp(['Erro 2:',num2str(err2)]);
disp(['Erro 3:',num2str(err3)]);

%% Task 2
clear all
close all
clc

delta = 0.015;
p = 0.5-delta;
q = 0.5+delta;

tf_matrix = [40,41; 400,401; 4000,4001];

N = 50000; % nr of trajectories
x = -500:500;
P_mean=zeros(length(x),3);

for g=1:3

    P = zeros(length(x),2);
    tf_list = tf_matrix(g,:);

    for n = 1:N
        for j = 1:length(tf_list)
            tf = tf_list(j);
            t=0:1:tf;
            pos=zeros(length(t),1);
            
            for i=1:length(t)-1    % i is the number of steps from 1 to 10
                r = rand();
        
                if r < p
                    S = -1;
        
                elseif r > p
                    S = 1;
                end
        
                pos(i+1) = S + pos(i); % nova posição

                if i == length(t)-1
                    xi = pos(i+1);
                    idx = xi + abs(min(x)) + 1;
                    P(idx,j) = P(idx,j)+1;
                end
            end
        end
    end

    for i=1:length(P)
        P_mean(i,g) = (P(i,1)+P(i,2))/2;
    
    end
    P_mean(:,g) = P_mean(:,g)./N;
end

plot(x,P_mean(:,1),'r',x,P_mean(:,2),'g',x,P_mean(:,3),'b','LineWidth',1.5)
title('<P(x,t)> of 1-D Asymmetric Random Walks')
legend(['t = ' num2str(tf_matrix(1,1)) ', ' num2str(tf_matrix(1,2))] ...
        ,['t = ' num2str(tf_matrix(2,1)) ', ' num2str(tf_matrix(2,2))] ...
        ,['t = ' num2str(tf_matrix(3,1)) ', ' num2str(tf_matrix(3,2))])
xlabel('x')
ylabel('Probabilities')
xlim([-200 350])
grid on

tf = (tf_list(1)+tf_list(2))/2;
P_theoretical = 1/sqrt(2*pi*tf).*exp(-((x-2*tf*delta).^2)./(2.*tf));

grid on, hold on
xlim([-200,350])
plot(x,P_mean(:,3),x,P_theoretical,'-','LineWidth',1.5)
title('Theoretical P(x,t) and <P(x,t)> (t=4000,4001)')
legend('<P(x,t)>','P_{t}')
xlabel('x')
ylabel('Probabilities')

% sum P =1 (checking 1st property in step 2 of the algorithm)
s1 = sum(P_mean(:,1));
s2 = sum(P_mean(:,2));
s3 = sum(P_mean(:,3));

%checking 2nd property in step 2 of the algorithm
s1_2 = sum(P_mean(:,1)).*(tf_matrix(1,1)+tf_matrix(1,2))*delta;
s2_2 = sum(P_mean(:,2).*(tf_matrix(2,1)+tf_matrix(2,2))*delta);
s3_2 = sum(P_mean(:,3).*(tf_matrix(3,1)+tf_matrix(3,2))*delta);

%checking 3rd property in step 2 of the algorithm
s1_3 = sum((P_mean(:,1)-(tf_matrix(1,1)+tf_matrix(1,2))*delta).^2);
s2_3 = sum((P_mean(:,2)-(tf_matrix(2,1)+tf_matrix(2,2))*delta).^2);
s3_3 = sum((P_mean(:,3)-(tf_matrix(3,1)+tf_matrix(3,2))*delta).^2);

%mean squared errors
err1 = immse(P_theoretical',P_mean(:,1)); 
err2 = immse(P_theoretical',P_mean(:,2)); 
err3 = immse(P_theoretical',P_mean(:,3)); 

%maximum values and comparing to the theory
max1 = max(P_mean(:,1));
max2 = max(P_mean(:,2));
max3 = max(P_mean(:,3));
max1t = 1/sqrt(2*pi*40.5);
max2t = 1/sqrt(2*pi*400.5);
max3t = 1/sqrt(2*pi*4000.5);


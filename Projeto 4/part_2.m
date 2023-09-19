close all
clear all
clc

N = 100;
t = 50000;
delta_T = 0.05;
T = 0.1:delta_T:10;

M = zeros(length(T),m);
MT = zeros(length(T),1);

parfor temp=1:length(T)
average_m_local = nan(1,t/interval-1);
sigma_old = ones(1,N);
sigma = sigma_old;
E_old = -J/N*sum(sum(sigma(1:end-1).*sigma(2:end)));

counter=1;
m = nan(1,t);
m_position = 1;
for a=1:t
    pos = round(1 + (N-1).*rand(1));
    sigma(pos) = -1.*sigma_old(pos);
    E_new = -J.*sum(sigma(1:end-1).*sigma(2:end))-H.*sum(sigma);
    delta_E = E_new-E_old;
    if delta_E <= 0
        sigma_old = sigma;
        E_old = E_new;
    else
        w = exp(-delta_E/T(temp));
        r = rand(1);
        if r<w
            sigma_old = sigma;
            E_old = E_new;
        end
    end
    m(a) = 1/N*sum(sigma_old);
    if mod(a,100) == 0 && a>100
        average_m_local(m_position) = 1/a*sum(m(1:a));
        m_position = m_position+1;
    end
end
average_m(temp,:) = average_m_local;
end
close all;
clear;
clc

%% Data
load("iddata-13.mat")
uid = id.InputData;
yid = id.OutputData;
Tsid = id.Ts;
tid = id.Tstart:Tsid:200;

uval = val.InputData;
yval = val.OutputData;
Tsval = val.Ts;
tval = val.Tstart:Tsval:150;

na = 1:5;
nb = 1:5;
nk = 1;
m = 1:5;
%% Alocare spatiu
tic
fprintf("Prealocare spatiu...\n")
%Prealocarea spatiului pentru regresori si a modelului simulat(pentru MSE)
Phi_id = cell(max(na), max(nb), max(m));
Phi_val = cell(max(na), max(nb), max(m));
y_sim = cell(max(na), max(nb), max(m));
y_sim_id=cell(max(na),max(nb),max(m));
MSE_pred = struct('id', zeros(max(na), max(nb), max(m)), 'val', zeros(max(na), max(nb), max(m)));
MSE_sim = struct('id', zeros(max(na), max(nb), max(m)), 'val', zeros(max(na), max(nb), max(m)));

for i = 1:max(na)
    for j = 1:max(nb)
        for k = 1:max(m)
            Phi_id{i,j,k} = calculPhi(yid, uid, i, j, nk, k);
            theta=pinv(Phi_id{i,j,k}) * yid;
            Phi_val{i,j,k} = calculPhi(yval, uval, i, j, nk, k);
            y_sim_id{i,j,k}=model_sim(yid, yid, uid, theta, i, j, nk,k);
        end
    end
end

fprintf("Spatiu alocat in:")
toc
fprintf("\n")
%% MSE
tic
fprintf("Calcul MSE...\n")
load("Date.mat");
load("y_sim_id.mat");
for i = 1:max(na)
    fprintf("i\n")
    for j = 1:max(nb)
        fprintf("j\n")
        for k = 1:max(m)
            fprintf("k\n")
            theta = pinv(Phi_id{i,j,k}) * yid;
            [MSE_pred_v, MSE_sim_v] = calculMSE(Phi_val{i,j,k},y_sim{i,j,k},yval, theta);
            MSE_pred.val(i,j,k) = MSE_pred_v;
            MSE_sim.val(i,j,k) = MSE_sim_v;
            [MSE_pred_id, MSE_sim_id] = calculMSE(Phi_id{i,j,k},y_sim_id{i,j,k},yid, theta);
            MSE_pred.id(i,j,k)=MSE_pred_id;
            MSE_sim.id(i,j,k)=MSE_sim_id;
        end 
    end
end
fprintf("MSE-uri calculate in:")
toc
fprintf('\n');
%% Heatmaps MSE Predictie
fprintf("Generare heatmaps...\n")
for k = 1:max(m)
    figure
    subplot(2,2,1);
    heatmap(nb, na, MSE_pred.id(:,:,k));
    xlabel('nb');
    ylabel('na');
    title(['Predictie MSE pe identificare(m = ', num2str(m(k)), ')']);
    colormap('hot');
    colorbar;

    % figure
    subplot(2,2,2);
    heatmap(nb, na, MSE_pred.val(:,:,k));
    xlabel('nb');
    ylabel('na');
    title(['Predictie MSE pe validare(m = ', num2str(m(k)), ')']);
    colormap('hot');
    colorbar;

    subplot(2,2,3);
    heatmap(nb, na, MSE_sim.id(:,:,k));
    xlabel('nb');
    ylabel('na');
    title(['Simulare MSE pe identificare(m = ', num2str(m(k)), ')']);
    colormap('hot');
    colorbar;

    subplot(2,2,4);
    heatmap(nb, na, MSE_sim.val(:,:,k));
    xlabel('nb');
    ylabel('na');
    title(['Simulare MSE pe validare(m = ', num2str(m(k)), ')']);
    colormap('hot');
    colorbar;
end
%%
%MODELUL ARX pentru na,nb si m egale cu 2(cel mai bun model posibil)
tic
[y_sim, y_pred] = calcul_ARX(yid, yval, uval,uid, 2, 2, nk, 2, Phi_id, Phi_val);
toc
%% Plots
% Data
figure;
subplot(211);
plot(tid(1:length(uid)), uid);
title("Identificare - Intrare");
subplot(212);
plot(tid(1:length(yid)), yid);
title("Identificare - Ieșire");

figure;
subplot(211);
plot(tval(1:length(uval)), uval);
title("Validare - Intrare");
subplot(212);
plot(tval(1:length(yval)), yval);
title("Validare - Ieșire");

% Models
figure;
plot(tval, yval, 'b', tval, y_pred.val, 'r--');
title("Predictie pe Validare");
legend("Sistem", "Model");

figure;
plot(tval, yval, 'b', tval, y_sim.val, 'r--');
title("Simulare pe Validare");
legend("Sistem", "Model");

figure;
plot(tid, yid, 'b', tid, y_pred.id, 'r--');
title("Predictie pe Identificare");
legend("Sistem", "Model");

figure;
plot(tid, yid, 'b', tid, y_sim.id, 'r--');
title("Simulare pe Identificare");
legend("Sistem", "Model");

figure;
subplot(211);
yval_pred = iddata(y_pred.val, uval, Tsval);
compare(val, yval_pred);
title("PREDICTIE");

subplot(212);
yval_sim = iddata(y_sim.val, uval, Tsval);
compare(val, yval_sim);
title("SIMULARE");

figure;
subplot(211);
yid_pred = iddata(y_pred.id, uid, Tsid);
compare(id, yid_pred);
title("PREDICTIE");

subplot(212);
yid_sim = iddata(y_sim.id, uid, Tsid);
compare(id, yid_sim);
title("SIMULARE");

toc
%% ARX
function [y_sim, y_pred] = calcul_ARX(yid, yval, uval,uid, na, nb, nk, m, Phi_i, Phi_v)
    maxim = max(na, nb + nk - 1);
    y_pred.val = zeros(length(uval), 1);
    y_pred.val(1:maxim) = yid(1:maxim);

    Phi_id = Phi_i{na, nb, m};
    Phi_val = Phi_v{na, nb, m};
    theta = pinv(Phi_id) * yid;
    y_pred.id = Phi_id * theta;
    
    for k = maxim+1:length(yval)
        y_pred.val(k) = Phi_val(k, :) * theta;
    end
    
    y_sim.val = model_sim(yid, yval, uval, theta, na, nb, nk, m);
    y_sim.id=model_sim(yid,yid,uid,theta,na,nb,nk,m);
end
%Modelul ARX pentru simulare
function y_sim = model_sim(yid, y, u, theta, na, nb, nk, m)
    maxim = max(na, nb + nk - 1);
    y_sim = zeros(length(y), 1);
    y_sim(1:maxim) = yid(1:maxim);
    
    for k = maxim+1:length(y_sim)
        Phi_sim = calculPhi(y_sim(1:k), u(1:k), na, nb, nk, m);
        y_sim(k) = Phi_sim(k, :) * theta;
    end
end
%Calculul regresorilor
function Phi = calculPhi(y, u, na, nb, nk, m)
    N = length(y);
    NrTermeni = factorial(m + na + nb - 1) / (factorial(m) * factorial(na + nb - 1)) + 1 + (na + nb + nk) * m;
    Phi = zeros(N, NrTermeni);
    Phi(:, 1) = 1;

    variabile = zeros(1, na + nb);
    combinari = Calcul_Exponenti(na + nb, m);
    
    for k = 2:N
        index = 2;
        for i = 1:na
            if k > i
                Phi(k, index) = -y(k-i);
                variabile(i) = y(k-i);
                index = index + 1;
            end
        end
        for i = 1:nb
            if k > i + nk - 1
                Phi(k, index) = u(k-i-nk+1);
                variabile(na+i) = u(k-i-nk+1);
                index = index + 1;
            end
        end
        for grad = 2:m
            for i = 1:na
                if k > i
                    Phi(k, index) = -y(k-i)^grad;
                    index = index + 1;
                end
            end
            for i = 1:nb
                if k > i + nk - 1
                    Phi(k, index) = u(k-nk-i+1)^grad;
                    index = index + 1;
                end
            end
        end
        
        nrCombinari = length(combinari(:, 1));
        nrColoane = length(combinari(1, :));
        for i = 1:nrCombinari
            Phi(k, index) = 1;
            for j = 1:nrColoane
                if combinari(i, j) > 0 && variabile(j) ~= 0
                    Phi(k, index) = Phi(k, index) * variabile(j)^combinari(i, j); 
                end
            end
            index = index + 1;
        end
        for grad = 1:m
            Phi(k, index) = sum(variabile)^grad;
            index = index + 1;
        end
    end
end
%Calcularea combinatiilor exponentilor multinomului lui Newton
function comb_exponenti = Calcul_Exponenti(n, m)
    comb_exponenti = [];
    exp_curent = zeros(1, n);
    %initial gradul maxim se afla la primul element
    exp_curent(1) = m;
    while true
        comb_exponenti = [comb_exponenti; exp_curent];
        %presupunem ca nu avem exponenti nenuli
        ok = false;
        for i = n:-1:2
            %caut exponenti diferiti de 0 in ordine inversa
            if exp_curent(i-1) ~= 0
                %actualizez exponenti
                exp_curent(i-1) = exp_curent(i-1) - 1;
                exp_curent(i) = exp_curent(i) + 1;
                %s-a gasit un exponent nenul
                ok = true;
                %recalibrare exponenti
                if i < n
                    exp_curent(i+1:end) = 0;
                    exp_curent(i) = m - sum(exp_curent(1:i-1));
                end
                break;%opresc cautarea
            end
        end
        %daca nu exista expoenenti nenuli se opreste calculul
        if ~ok
            break;
        end
    end
end

%% MSE
function [MSE_pred, MSE_sim] = calculMSE(Phi,y_sim, yval, theta)
    y_pred = Phi * theta;
    e_pred = yval - y_pred;
    MSE_pred = sum(e_pred.^2) / length(yval);
    e_sim = yval - y_sim;
    MSE_sim = sum(e_sim.^2) / length(yval);
end
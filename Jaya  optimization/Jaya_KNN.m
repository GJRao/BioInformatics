
clear, clc, close;
% Number of k in K-nearest neighbor
opts.k = 5; 
% Ratio of validation data
ho = 0.20;
% Common parameter settings 
opts.N = 10;     % number of solutions
opts.T = 100;    % maximum number of iterations

% Load dataset
data = readtable('new_data.csv','PreserveVariableNames',true);

feat=data{:,2:177};
Y=data{:,end};
% Y=categorical(Y)
% us = unique(Y);
% %Applying SMOTE
% rng(42)
% [data1,Y1,new_feat,new_label] = smote(data1, [], 'Class', Y);

% Divide data into training and validation sets
% feat=new_feat(:,:);
% label=new_label;


label=Y;
HO = cvpartition(label,'HoldOut',ho); 
opts.Model = HO; 
% Perform feature selection 
FS     = jfs('ja',feat,label,opts);
% Define index of selected features
sf_idx = FS.sf;
FS;
% Accuracy  
Acc    = jknn(feat(:,sf_idx),label,opts); 
% Plot convergence
figure;
plot(FS.c); grid on;
xlabel('Number of Iterations'); 
ylabel('Fitness Value'); 
title('JA fitness graph');

fprintf("Acc %d::",Acc)
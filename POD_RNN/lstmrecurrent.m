% POF: 2020
% author: R.Gupta and S.R.Bukka 
% POD-RNN hybrid model 

%% make sure that you have run the script "pod.m" 
% note that predicted coefficients are already provided in the directory:
% (you can skip this-step or generate your own training and prediction)
% "y_pred_n_pres.mat" : time-series prediction (pressure) 
% "y_pred_n_velx.mat" : time-series prediction (x-velocity)

%% normalise the data 
a_10 = a1(1:5,:);
a_10 = a_10';
for i = 1:5
    am(:,i) = (a_10(:,i)-min(a_10(:,i)))./(max(a_10(:,i)-min(a_10(:,i))));
end
am = am(1:4500,:);

%% data organisation in train and test 
% snapshots for analysis = [501:5000]
% traning data = [501:3500]   % 3000 for training  
% testing data = [3501:4500]  % 1000 for testing 
% ts = 0.25 s
n_train = 3000;              
n_test = 1000; 
total = n_train + n_test; 

X_train = am(1:n_train-1,:);      % 1 till 2999
Y_train = am(2:n_train,:);        % 2 till 3000
X_test = am(n_train+1:total,:);   % 3001 till 4000  

X_train = X_train';
Y_train = Y_train';
X_test = X_test';

%% define LSTM-RNN network 
numFeatures = 5;
numResponses = 5;
numHiddenUnits = 512;

layers = [ ...
    sequenceInputLayer(numFeatures)               % input 
    lstmLayer(numHiddenUnits)                     % hidden 
    fullyConnectedLayer(numResponses)             % fc 
    regressionLayer];                             % regression 

options = trainingOptions('adam', ...             % optimizer 
    'MaxEpochs',3000, ...                         % epochs     
    'GradientThreshold',1, ...                    % max gradinet in back-prop    
    'InitialLearnRate',0.005, ...                 % initial learning rate 
    'LearnRateSchedule','piecewise', ...          % decay learing rate...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');

%% do some training! 
net_sbs = trainNetwork(X_train,Y_train,layers,options);
net_sbs = predictAndUpdateState(net_sbs,X_train);

%% prediction 
% take time-step 3500 and get 3501 
[net_sbs,YPred] = predictAndUpdateState(net_sbs,Y_train(:,end));

% do close loop-type prediction (3502-4500)
for i = 2:1000
    [net_sbs,YPred(:,i)] = predictAndUpdateState(net_sbs,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

% rmse of predicted and actual  
 rmse = X_test(1:5,1:end-1)-YPred(:,1:end-1);
for i = 1:5
    rmse_updated(i) = sqrt(mean(rmse(i,:).^2));
end

% denormalise the prediction 
y_pred_n = zeros(1000,5);
for i = 1:5
    y_pred_n(:,i) = YPred(i,:)'.*(max(a_10(:,i)-min(a_10(:,i)))) + min(a_10(:,i));
end

%% save 
save('y_pred_n.mat','y_pred_n'); 
  
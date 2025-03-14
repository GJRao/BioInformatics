% K-nearest Neighbor (9/12/2020)

function Acc = jknn(feat,label,opts)
% Default of k-value
%k = 7;

if isfield(opts,'k'), k = opts.k; end
if isfield(opts,'Model'), Model = opts.Model; end

% Define training & validation sets
trainIdx = Model.training;    testIdx = Model.test;
xtrain   = feat(trainIdx,:);  ytrain  = label(trainIdx);
xvalid   = feat(testIdx,:);   yvalid  = label(testIdx);


% Training model
 My_Model = fitcknn(xtrain,ytrain,'NumNeighbors',k); 

%My_Model = fitcknn(feat,label,'NumNeighbors',k); 

%Calculate the misclassification error and the classification accuracy on the training data.
trainError = resubLoss(My_Model);
trainAccuracy = 1-trainError;
fprintf('\n Train Accuracy: %g %%',100 * trainAccuracy);

%Create a partitioned model cvMdl. 
% Compute the 10-fold cross-validation misclassification error and classification accuracy.
cvMdl = crossval(My_Model); % Performs stratified 10-fold cross-validation
cvtrainError = kfoldLoss(cvMdl);
trainAccuracy = 1-cvtrainError;
fprintf('\n CV model Train Accuracy: %g %%',100 * trainAccuracy);

% Prediction with CV model
testError = loss(My_Model,xvalid,yvalid);
newAccuracy = 1-testError;
fprintf('\n new Accuracy: %g %%',100 * newAccuracy);

% Prediction
% [pred,score,cost]   = predict(My_Model,xvalid);

[pred,score]   = predict(My_Model,xvalid);
% Accuracy
%Acc      = sum(pred == yvalid) / length(yvalid);
Acc      = sum(strcmp(pred, yvalid)) / length(yvalid);
fprintf('\n Test Accuracy: %g %%',100 * Acc);



% score
% [~,score_knn] = resubPredict(My_Model);
% score_knn
C = confusionmat(pred , yvalid);
%[c_matrix,Result,RefereceResult]= confusion.getMatrix(pred , yvalid);
% loss(My_Model,xvalid,yvalid);
figure
cm = confusionchart(C);
title('Knn JA Confusion Matrix');
% pred
% yvalid
% score = score(:)';
classNames=My_Model.ClassNames;

rocObj = rocmetrics(yvalid,score,classNames);
%returns the optimal operating point of the ROC curve.
%imagesc(C);
%colorbar;
rocObj.AUC  ;
figure
plot(rocObj);
title('Knn JA ROC');
end



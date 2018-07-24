function ANNtrainParms = set_ann_parms

% ann training parameters
ANNtrainParms.verbose       = 0;
% trainParms.nodes           = 10;
ANNtrainParms.trainRatio    = 0.65;
ANNtrainParms.max_fail      = 50;
ANNtrainParms.epochs        = 1e3;
ANNtrainParms.trainFcn      = 'trainscg';
ANNtrainParms.performFcn    = 'mse';


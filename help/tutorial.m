%% Tutorial
% 
% Requirements:
%     - optimization toolbox (fminunc)

import regression.*

%% Load the raw data
rawData = load('exampleData.mat'); % run tutorial_exampleData to generate this
nTrials = rawData.nTrials; % number of trials
trial = rawData.trial;

% Settings
unitOfTime = 'ms';
binSize = 1; % TODO some continuous observations might need up/down-sampling if binSize is not 1!?


%% Build 'neuroGLM' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.

nGLM = neuroGLM(unitOfTime, binSize, rawData.param);

% 2-D eye position
bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 40, 4, nGLM.binfun);
bs.normalize

nGLM.addCovariateRaw(trial, 'eyepos', 'eyeposition', bs);

% Spike timing
nGLM.addCovariateSpiketrain(trial, 'hist', 'sptrain', 'History Filter')
nGLM.addCovariateSpiketrain(trial, 'coupling', 'sptrain2', 'Coupling from neuron 2')

% Instantaneous Raw Signal without basis
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, nGLM.binfun);
bs.normalize

nGLM.addCovariateRaw(trial, 'LFP', 'LFP', bs);

% Acausal Timing Event
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 8, nGLM.binfun);
offset = -200;
nGLM.addCovariateTiming(trial, 'saccade', 'saccade', [], bs, offset);

% Coherence: A box car that depends on the coh value
bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 200, 10, nGLM.binfun);
bs.normalize

stimHandle = @(trial) trial.coh * basisFactory.boxcarStim( nGLM.binfun(trial.dotson), nGLM.binfun(trial.dotsoff), nGLM.binfun(trial.duration));
nGLM.addCovariate(trial, 'cohKer', 'coh-dep dots stimulus', stimHandle, bs);




%% Compile the data into 'DesignMatrix' structure
trialIndices = 1:(nTrials-1); % use all trials except the last one
nGLM.compileDesignMatrix(trial, trialIndices);


%% Get the spike trains back to regress against
y = nGLM.getBinnedSpikeTrain(trial, 'sptrain', nGLM.dm.trialIndices);

%% Do some processing on the design matrix
nGLM.removeConstantCols;
colIndices = nGLM.getDesignMatrixColIndices('LFP');
nGLM.zscoreDesignMatrix([colIndices{:}]);

nGLM.addBiasColumn; % DO NOT ADD THE BIAS TERM IF USING GLMFIT

%% Least squares for initialization
tic
wInit = nGLM.dm.X \ y;
toc

%% Use matRegress for Poisson regression

fnlin = @expfun; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)neglogli_poiss(w, nGLM.dm.X, y, fnlin, 1); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wvar = diag(inv(hessian));

%% Visualize
ws = nGLM.combineWeights(wml);
wvar = nGLM.combineWeights(wvar);

fig = figure(2913); clf;
nCovar = numel(nGLM.covar);
for kCov = 1:nCovar
    label = nGLM.covar(kCov).label;
    subplot(nCovar, 1, kCov);
    plot(ws.(label).tr, ws.(label).data, ...
        ws.(label).tr, ws.(label).data+sqrt(wvar.(label).data), '--', ...
        ws.(label).tr, ws.(label).data-sqrt(wvar.(label).data), '--')
    title(label);
end


# Introduction

Neural recordings are accompanied by various other measurements, controlled manipulations, and other meta data.
For example, timing of a visual cue, time course of an auditory stimulus, behavioral measurements (button press, eye movement, EMG), time and magnitude of the reward, heart beat, on going EEG/ECoG/LFP measurements, and so on.
It is cumbersome to find a representation of these covariates for regression-style analysis, especially when they are events.
This is because events usually do not have an **instantaneous effect to the response (dependent) variable**.
This package allows you to **expand and transform your experimental variables to a feature space as a [design matrix](http://en.wikipedia.org/wiki/Design_matrix) such that a simple linear analysis could yield desired results**.

This tutorial will explain how to import your experimental data, and build appropriate features spaces to do fancy regression analysis easily.

# Loading your data

We assume the experimental variables are observations over time, and organized into **trials**.
If your data don't have a trial structure, you can put all your data in a single trial.
There are 4 types of variables that constitute data: *spike train*, *timing*, *continuous*, and *value*.
This framework uses string labels to address each variable later.

## Types of experimental variables
### Spike train

Each spike train is a sequence of spike timings from a single neuron.
The spike timings are relative to the beginning of the trial.

### Timing (events)

A typical trial based experiment may have a cue that indicates the beginning of a trial, cues that indicate waiting period, or presentation of a target at random times.
These covariates are best represented as events.
However, if two or more timings are perfectly correlated, it would be sufficient to include just one.
Note that many behaviors are also recorded as timing: reaction time, button press time, etc.

### Continuous

Continuous data are measured continuously over time.
For instance, eye position of the animal may be correlated with the activities of neurons of the study, so the experimentalist could carefully measure it throughout the experiment.
Note that the sampling rate should match the bin size of the analysis, otherwise up-sampling, or down-sampling (with appropriate filtering) is necessary.

### Value

Each trial can have a single value associated with it.
In many cases these are trial specific parameters such as strength of the stimulus, type of cue, or the behavioral category.
These values can be used to build a feature space, or  to include specific feature in trials only when certain conditions are met.


## Preparing the trial data structure

Experimental data needs to be organized into a structure ('trial') where you can need to add each of your experimental variables you have registered for the experiment as fields. Below are examples with randomly generated dummy data (see [tutorial_exampleData.m](../help/tutorial_exampleData.m)).

```matlab
trial = struct();
trial(nTrials).duration = 0; % preallocate

lambda = 0.1; % firing rate per bin

for kTrial = 1:nTrials
    duration = 150 + ceil(rand * 200);
    trial(kTrial).duration = duration;
    trial(kTrial).LFP = cumsum(randn(duration, 1));
    trial(kTrial).eyepos = cumsum(randn(duration, 2));
    trial(kTrial).dotson = ceil(rand * (duration - 100));
    trial(kTrial).dotsoff = trial(kTrial).dotson + ceil(rand * 10);
    trial(kTrial).saccade = 100 + ceil(rand * (duration - 100));
    trial(kTrial).coh = sign(rand - 0.5) * 2^ceil(rand*8);
    trial(kTrial).choice = round(rand);
    trial(kTrial).sptrain = sort(rand(poissrnd(lambda * 0.9 * duration), 1) * duration);
    trial(kTrial).sptrain2 = sort(rand(poissrnd(0.1 * lambda * duration), 1) * duration);

    trial(kTrial).sptrain = sort([trial(kTrial).sptrain; trial(kTrial).sptrain2 + 2]);
    trial(kTrial).sptrain(trial(kTrial).sptrain > trial(kTrial).duration) = [];

    trial(kTrial).meta = rand;
end
```

## The neuroGLM class

In the revisions to the repository, a dedicated neuroGLM class has been created by jcbyts to combine the definition of experimental objects, loading of trial data and specification of experimental design.

First, create an instance of the `neuroGLM` class:
```matlab
nGLM = neuroGLM(unitOfTime, binSize, rawData.param);
```
where `unitOfTime` is a string for the time unit that's going to be used consistently throughout (e.g., 's' or 'ms'), `binSize` is the duration of the time bin to discretize the timings.
`param` can be anything that you want to associate with the experiment structure for easy access later, since it will be carried around throughout the code.

Next add covariates associated with each experimental variable. We'll start here with the spike trains, because the basis functions for these variables are automatically defined. To add to the model, we therefore just pass the trial structure (`trial`), the type, label, and user friendly name of the variable. Note that the label must match the name of the field in the trial structure that you want to pass (e.g. `sptrain`):

```matlab
nGLM.addCovariateSpiketrain(trial, 'hist', 'sptrain', 'History Filter')
nGLM.addCovariateSpiketrain(trial, 'coupling', 'sptrain2', 'Coupling from neuron 2')
```

For other variables, we need to define the basis functions to fit using an instance of the `Basis` class:

```matlab
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 8, nGLM.binfun);
offset = -200;
nGLM.addCovariateTiming(trial, 'saccade', 'saccade', [], bs, offset);
```

<i>Add note about basis normalize</i>



## More on temporal basis functions
The `+basisFactory` package provides functions that generate basis function structures.
Instead of a boxcar, you can use raised cosine basis functions:
```matlab
bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 200, 10, expt.binfun);
```
The raised cosine functions are spaced such that they can represent linear functions as well as smoothly varying functions.

## Building the design matrix
To obtain the design matrix, we need to run the compile function:
```matlab
nGLM.compileDesignMatrix(trial, trialIndices);
```
where `trialIndices` are the trials to include in making the design matrix. This function is memory intensive, and could take a few seconds to complete. Once it does complete, the neuroGLM object will have a `dm` attribute containing the actual design matrix as `nGLM.dm.X`

If your design matrix is not very sparse (less than 10% sparse, for example), it's better to convert the design matrix to a full (dense) matrix for speed.
```matlab
dm.X = full(dm.X);
```

# Regression analysis
Once you have designed your features, and obtained the design matrix, it's finally time to do some analysis!

## Get the dependent variable
You need to obtain the response variable of the same length as the number of rows in the design matrix to do regression. For **point process** regression, where we want to predict the observed spike train from covariates, this would be a finely binned spike train concatenated over the trials of interest:
```matlab
y = nGLM.getBinnedSpikeTrain(trial, 'sptrain', nGLM.dm.trialIndices);
```

## Doing the actual regression

There are multiple ways to do the regression (least-squares, glmfit etc.) but the tutorial will use the fminunc function to optimize maximum likelihood:

```matlab
fnlin = @expfun; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)neglogli_poiss(w, nGLM.dm.X, y, fnlin, 1); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wvar = diag(inv(hessian));
```

# Post regression weight reconstruction
Result of regression is a weight vector (and sometimes additional associated statistics in a vector or matrix) in the feature space. Hence, the weight vector is as long as the number of columns in the design matrix. In order to obtain meaningful temporal weights back corresponding to each covariate, use
```matlab
ws = nGLM.combineWeights(wml);
wvar = nGLM.combineWeights(wvar);
```
It returns a structure that contains a time axis `ws.(covLabel).tr` and data `ws.(covLabel).data` for each `covLabel` in the design specification structure.

# Initialize problems
This explains how to initialize all class objects 

## Model
TODO

## Uncertainty
The "uncertainity" class helps to define object pertaining to different noises. As of now, we support "Uniform", "Gaussian", "Laplace" and more will be supported in the future. 

```
An example where we need no noise

modelUncertainty = datools.uncertainty.NoUncertainty;


Another example where we initialize an object which can be used to add Gaussian noise

modelUncertainity = datools.uncertainity.Gaussian('Covariance', R);
where `R`  is the observation error covariance
``` 

## Observation
TODO

## Tapering
TODO

## Filter
The most important part of using this package is to define a "filter" object. As of now, we support several ensemble, variational and gaussian mixture based filters/smoothers(provide reference to the folder).

All filters are derived from the base `DABase' class 

```
% an example

filter = datools.filter.ensemble.EnKF(model, ...
		'InitialEnsemble', ensembleGenerator(N), ...
		'Inflation', inflation, ...
		'Parallel', false)
```

Here `model` is an object of the model class that determines the model operator (`M`) and defines a solver to propogate the model. We have the flexibility for the users to provide the Tangent Linear/Adjoint of `M` (do we?)
 

# Initialize problems
This explains how to initialize all class objects 

## Model
TODO

## Uncertainty
TODO

## Observation
The observation class defines the observation operator(H). It also has the flexibility to accept the jacobian of H. The user needs to supply the jacobian of H for now.
Suppose we need to have indexed observation, we do the following:
```
% example for defining object of observation class

observation = datools.observation.Indexed(model.NumVars,...
			'Uncertainty', observationErrorModel,...
			'Indices', indicesObserved)
``` 
Here `model` is an object of the model class, `observationErrorModel' is an object of the error class, `indicesObserved` are the indices of the states observed.

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
 

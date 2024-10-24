# Contributing to DATools

General Guidelines for adding new filters/algorithms to DATools

## Submitting new pull requests

To request a pull request, please follow the steps below

* Create your own branch if you are starting new or modify the existing one
* Commit your changes to local branch and push to Github
* Create a draft pull request
* Do not merge or modify the master  branch

## Creating and running problems
TODO

## Style Conventions

In order to prevent confusions and maintain a consistent coding style, the following styles should be strictly followed.
The conventions match most of the MATLAB's own conventions. Failing to follow the conventions can delay in processing 
and approving new pull requests.

### Line Formatting

Four spaces are used for indentation where required. A line should be kept less than 120 characters.

### Variables

Variable names should follow camelCasing convention.

```
% Example
N = 100;
xValue = linspace(-5,5,N)
zValue = @(x, sigma, n) exp(-x) + sig*randn(1,n)
```

### Functions

Function names should be in camelCase. No special characters (for example '_') shall be used.

```
% Example
function val = myFun(args1, args2, ...)
    ...
end
```


### Structures

Structures should have camelCase property names 

```
% Example 
filter = struct('Type', 'Ensemble', 'Name', 'EnKF')
```


### Classes

Class names and properties should be in PascalCase. Acronyms should be capitalized. Methods will have camelCase as before. 

```
% Example
classdef ExampleClass < handle
    properties
        PropertyOne
        PropertyTwo
        PropertyACRONYM
    end

    methods
        function val = fun(...)
            ...
        end
    end
end

```



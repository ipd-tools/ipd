# Summarize ipd fit

Summarize ipd fit

## Usage

``` r
# S3 method for class 'ipd'
summary(object, ...)
```

## Arguments

- object:

  An object of class `ipd`.

- ...:

  Ignored.

## Value

An object of class `summary.ipd` containing:

- call:

  The model formula.

- coefficients:

  A glm-style table of estimates, SE, z, p.

- method:

  Which IPD method was used.

- model:

  Which downstream model was fitted.

- intercept:

  Logical; whether an intercept was included.

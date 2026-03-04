# ipd: S4 class for inference on predicted data results

ipd: S4 class for inference on predicted data results

## Slots

- `coefficients`:

  Numeric vector of parameter estimates.

- `se`:

  Numeric vector of standard errors.

- `ci`:

  Numeric matrix of confidence intervals.

- `coefTable`:

  Data frame summarizing `Estimate`, `Std. Error`, `z value`, and
  `Pr(>|z|)`.

- `fit`:

  The raw list returned by the helper function.

- `formula`:

  The formula used (class "formula").

- `data_l`:

  The labeled data (data.frame).

- `data_u`:

  The unlabeled data (data.frame).

- `method`:

  Character; which IPD method was used.

- `model`:

  Character; which downstream model was fitted.

- `intercept`:

  Logical; was an intercept included?

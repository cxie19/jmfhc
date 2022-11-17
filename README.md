# jmfhc

jmfhc is an R package to fit a new cure model with flexible hazard ratios 
between different covariate subgroups in longitudinal and cure-survival data. 
It is a joint model with individual random effects shared between cure-survival 
and longitudinal submodels and accommodates covariates’ non-proportional hazards 
structures easily. More specifically, this R package can fit a joint model using 
a flexible-hazards cure submodel (JMFHC) and a joint model using a proportional 
hazards submodel (JMPHC) as JMFHC's special case for longitudinal 
and right-censored survival data with a cure fraction. This package also provides 
functions to estimate a patient subgroup's cure rate and survival function 
based on the fitted model. <br />

## What is JMFHC

Our proposed JMFHC consists of a linear mixed-effects model as the longitudinal 
submodel and a flexible-hazards cure (FHC) model as the cure-survival submodel.
An FHC model is in the form of the promotion time cure model and an extension of 
the proportional hazards cure model. The FHC model can be fitted by the R package
fhc at https://github.com/cxie19/fhc.

For subject $i$, the observed values of the biomarker at measurement times 
$\boldsymbol{t_i}$ denoted as $\boldsymbol{b_i}$ are shown as 

$$
  \boldsymbol{b_i}(\boldsymbol{t_i})=
  \boldsymbol{D_i}\boldsymbol{\phi}+\boldsymbol{D_i}\boldsymbol{\alpha_i}+\boldsymbol{\epsilon_i},
$$

where $\boldsymbol{D_i}$ is a $n_i \times p^\ast$ design matrix for fixed effects 
with the first column as 1's and the reminding $(p^\ast-1)$ columns containing 
a biomarker's measurement time points using fractional polynomials (e.g., $\boldsymbol{t_i},\log(\boldsymbol{t_i}),\boldsymbol{t_i}^2$);
$\boldsymbol{\phi}$ is a $p^\ast$-length fixed-effect regression parameter vector;
$\boldsymbol{\alpha_i}$ is a $p^\ast$-length random-effect regression parameter 
vector following a multivariate normal distribution $N_{p\ast}(\boldsymbol{0},\boldsymbol{\Sigma})$ and the unstructured 
covariance matrix $\boldsymbol{\Sigma}$ containing elements 
of $\sigma_1,...,\sigma_{p^\ast}$ and $\rho_{jm}$, for $j,m = 1,...,p^\ast$ and 
$j \neq m$; 
and $\boldsymbol{\epsilon_i}$ is a $n_i$-length vector of measurement errors 
following a multivariate normal distribution $N_{n_i}(\boldsymbol{0},\boldsymbol{R_i}=\sigma_{\epsilon}^2\boldsymbol{I_{n_i}})$.
Here all $\boldsymbol{\alpha_i}$ and $\boldsymbol{\epsilon_i}$ are 
mutually independent.

The survival function in our joint model (i.e., JMFHC) for subject $i$ at time 
$T$ conditional on covariates $\boldsymbol{z_i}$ and $\boldsymbol{x_i}$ 
and random effects $\boldsymbol{\alpha_i}$ shared with the longitudinal submodel
is assumed to be

$$
    S(t|\boldsymbol{x_i},\boldsymbol{z_i},\boldsymbol{\alpha_i})=\exp\left[-e^{\beta_0}e^{\boldsymbol{z_i'\psi}}e^{\boldsymbol{\alpha_i'\eta}}\{F_0(t)\}^{\exp(\boldsymbol{x_i'\gamma})}\right],
$$

where $\beta_0$ is an unknown scalar; 
$\boldsymbol{\psi}$ and $\boldsymbol{\eta}$ are two vectors of unknown 
regression parameters with lengths $p$ and $p^\ast$ for baseline covariates 
$\boldsymbol{z_i}$ (long-term covariates) and random effects 
$\boldsymbol{\alpha_i}$, respectively;
$F_0(t)$ is a monotone increasing function with $F_0(0)=0$ 
and $\lim_{t\to\infty}F_0(t)=1$;
and $\boldsymbol{\gamma}$ is a vector of unknown regression parameters 
with length $q$ for short-term baseline covariates $\boldsymbol{x_i}$.
JMPHC has $\boldsymbol{\gamma}=\boldsymbol{0}$ in the cure submodel. 


## How to get started

Install the R package using the following commands on the R console:

```{r}
install.packages("devtools")
devtools::install_github("cxie19/jmfhc")
library(jmfhc)
```

An example data set called *jmfhc_dat* is provided in this package. It
is a longitudinal and cure-survival data set. Its documentation can be 
seen by using the following command.

```{r}
help(jmfhc_dat)
```

The function *jmfhc_point_est* is called to obtain point estimation of 
regression parameters and $F_0(t)$ in JMFHC or JMPHC. 
Its documentation can be seen by using the following command.

```{r}
help(jmfhc_point_est)
```

The function *jmfhc_se_est* is called to obtain standard error estimation of
regression parameters and $F_0(t)$ in a JMFHC or JMPHC. Its documentation can 
be seen by using the following command.

```{r}
help(jmfhc_se_est)
```

The function *estcure* is called to compute a cure rate for a patient subgroup 
with specified characteristics given a fitted JMFHC or JMPHC.
Its documentation can be seen by using the following command.

```{r}
help(estcure)
```

The function *est_surv_func* is called to obtain an estimated survival
function for a patient subgroup with specified characteristics given a 
fitted JMFHC or JMPHC. Its documentation can be seen by using the
following command.

```{r}
help(est_surv_func)
```

## Example
For example, we want to fit a JMFHC for the example data *jmfhc_dat*.
This joint model has repeatedly measured biomarker values as the outcome of the 
longitudinal submodel with measurement times as the 
covariate and treatment as the short- and long-term covariate in the cure 
submodel. These two submodels share individual random effects.
We call the function *fhcmodel*, and the following command is used.

```{r}
result_coef <- jmfhc_point_est(data=jmfhc_dat, 
                               event_time="event.time", event_status="event", 
                               id="patient.id", 
                               beta_variable="trt", gamma_variable="trt", 
                               fu_measure="measure", fu_time_variable="mes.times")
```
The running time of this point estimation (result_coef) is about 1.25 hours.<br /> 
With the estimated regression parameters from *result_coef*, 
we compute the cure rates for patient subgroups 
receiving treatment A and treatment B.
We call the function *estcure*, and the following 
commands are used.

```{r}
estcure(object=result_coef,z_value=0) # treatment A
estcure(object=result_coef,z_value=1) # treatment B 
```

The cure rates of patients receiving treatment A and treatment B are 41.48% and 
10.32%, respectively.<br /> 
Meanwhile, with the estimated regression parameters and $F_0(t)$ from 
*result_coef*,  we estimate the survival functions of these two subgroups 
of patients.

```{r}
# treatment A
survival_trt0 <- est_surv_func(object=result_coef,z_value=0,x_value=0,event_time="event.time") 
# treatment B
survival_trt1 <- est_surv_func(object=result_coef,z_value=1,x_value=1,event_time="event.time") 
```



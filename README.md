# rssc-state-space-models
Source code for:

Sansom, P. G., Williamson, D. B., & Stephenson, D. B. (2019). State space models for non-stationary intermittently coupled systems. Journal of the Royal Statistical Society: Series C (Applied Statistics), 68, 1259–1280. https://doi.org/10.1111/rssc.12354


The files included reproduce all the analysis and figures in both the main text
and the supplementary material. The main files are "ekf.R" which contains our
implementation of the Extended Kalman Filter, "make-model.R" which sets up the
observation and evolution matrices (F,G) and covariance matrices (V,W), 
"metropolis-hastings.R" performs adaptive Markov Chain Monte Carlo sampling of
the variance and intervention parameters, "trajectories.R" samples state
trajectories by backward sampling once MCMC sampling is complete for efficiency.


General files:
ekf.R                 - The Extended Kalman Filter
                        (filtering,smoothing,sampling,prediction,etc.)
era-interim.csv       - The NAO time series data
make-model.R          - Setup the Extended Kalman Filter
metropolis-hastings.R - Adaptive Markov Chain Monte Carlo
trajectories.R        - Sample state trajectories by backward-sampling

The mean intervention model:
nao-mean.R                 - MCMC sampling of the variance parameters etc.
nao-mean-bridge-sampling.R - Compute the marginal likelihood by bridge-sampling
nao-mean-trajectories.R    - Sample state trajectories by backward sampling
nao-mean-predictive.R      - Posterior predictive simulation

The autocorrelation intervention model:
nao-tvar.R                 - MCMC sampling of the variance parameters
nao-tvar-bridge-sampling.R - Compute the marginal likelihood by bridge-sampling
nao-tvar-trajectories.R    - Sample state trajectories by backward sampling
nao-tvar-predictive.R      - Posterior predictive simulation

Main figures and analysis:
acf.R      - Function to compute multiple ACFs
anova.R    - Perform the analysis of variance in Table 4
figure1.R  - Produce Figure 1
figure2.R  - Produce Figure 2
figure3.R  - Produce Figure 3
figure4.R  - Produce Figure 4
figure5.R  - Produce Figure 5
figure6.R  - Produce Figure 6
figure7a.R - Produce Figure 7a
figure7b.R - Produce Figure 7b
figure8a.R - Produce Figure 8a
figure8b.R - Produce Figure 8b

NB: By default figures 3,4,5 & 6 are produced for the mean and autocorrelation
intervention models, with both the standard and alternative priors. The output
can be limited by changing the "type" loop at the start of each file.


Supplementary figures and analysis:

The mean intervention model with alternative priors:
nao-mean-alt.R                 - MCMC sampling of the variance parameters
nao-mean-alt-bridge-sampling.R - Compute the marginal likelihood
nao-mean-alt-trajectories.R    - Sample state trajectories by backward sampling
nao-mean-alt-predictive.R      - Posterior predictive simulation

The autocorrelation intervention model with alternative priors:
nao-tvar-alt.R                 - MCMC sampling of the variance parameters
nao-tvar-alt-bridge-sampling.R - Compute the marginal likelihood
nao-tvar-alt-trajectories.R    - Sample state trajectories by backward sampling
nao-tvar-alt-predictive.R      - Posterior predictive simulation

The null model without any interventions:
nao-null.R                     - MCMC sampling of the variance parameters
nao-null-bridge-sampling.R     - Compute the marginal likelihood
nao-null-trajectories.R        - Sample state trajectories by backward sampling
nao-null-predictive.R          - Posterior predictive simulation

The null model without any interventions with alternative priors:
nao-null-alt.R                 - MCMC sampling of the variance parameters
nao-null-alt-bridge-sampling.R - Compute the marginal likelihood
nao-null-alt-trajectories.R    - Sample state trajectories by backward sampling
nao-null-alt-predictive.R      - Posterior predictive simulation

Simulation studies:
durbin-levinson.R            - Functions to perform Durbin-Levinson recursions
simulation-study-mean.R      - Mean intervention model
simulation-study-null.R      - Null model
simulation-study-tvar.R      - Fast evolving autocorrelation intervention model
simulation-study-tvar-good.R - Slow evolving autocorrelation intervention model

Posterior traces:
traces-mean.R - Mean intervention model (Supp. Fig. 4)
traces-null.R - Null model
traces-tvar.R - Autocorrelation intervention model (Supp. Fig. 5)

Posterior histograms:
posterior-mean.R     - Mean intervention model (Supp. Fig. 6)
posterior-mean-alt.R - Mean intervention model with alt. priors (Supp. Fig. 8)
posterior-tvar.R     - Autocorrelation intervention model (Supp. Fig. 7)
posterior-tvar-alt.R - Autocorrelation intervention model with alt. priors (9)

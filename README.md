# MCMC_unbiased_gibbs_sampler


## Presentation

Implementation of a Gibbs Sampler from the paper Unbiased Markov chain Monte Carlo with couplings by Pierre E. Jacob, John O’Leary , Yves F. Atchadé in Ocaml. 


  

## Normal Bayesian estimation

  

To launch the estimation, compile the launch_mcmc file:

```

ocamlfind ocamlopt -o launch_mcmc -linkpkg -package owl,owl-plplot gamma.ml randomx.ml mcmc.ml launch_mcmc.ml


```

And then launch it:

```

QT_QPA_PLATFORM=offscreen ./launch_mcmc 

```
  
  

## Coupled Bayesian estimation

  
To launch the estimation, compile the launch_mcmc_coupled file:

```

ocamlfind ocamlopt -o launch_mcmc_coupled -linkpkg -package owl,owl-plplot gamma.ml randomx.ml mcmc.ml launch_mcmc_coupled.ml

```

And then launch it:

```

QT_QPA_PLATFORM=offscreen ./launch_mcmc_coupled

```

Some changes are still to be done.
BCR analysis
--------------------------

## TODO
* [ ] vars_to_save feature with the main script

## Workflow
1. First, this is the tricky bit: Configure the analysis as a .m file. See a predefined one as a template.
2. Nice. Now edit the appropriate config line in `main_shockfringe.m`, which is the main analysis script.
3. Run it.

### Optional
* summary scripts are there but these were written particularly to generate figures for publication and is not intuitive

## Scratchpad
### Results
Analysis results can be saved automatically (see configs). In fact, these will be used by the summarising scripts.

Bare-bone results from main analysis to save (required by summary codes) are listed below.

#### summary_machdyn
We can predict how fast the BEC was zipping through the atom laser.

```matlab
varsummary={'r_cent','tof','g','nden_r','t0','c_const','pal_R','pal_nseq'};
```

#### summary_fringe
The fringe spacings can tell us about the physics of supersonic collision in superfluids.

```matlab
varsummary={'N_peak_max', 'lambda_ff', 'lambda_ff_err', 'Nal', 'Nal_err_tot', 'N0', 'N0_err_fit',...
        'v','c','lambda_nf', 'eff_al'};
```
```matlab	
loadvars={'dn1d','d_1d','ppeak','pal_R'};
```
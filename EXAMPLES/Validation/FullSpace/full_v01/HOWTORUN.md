# How to run the Validation

Notation `P` is input for point source as true model, and `E` is for external source model
1. run the `run_xmeshfem2D.sh` and check the location of external sources in `OUTPUT_FILES/gridfile.ps`.
> The red boxes should form closed surface.

2. make the station files for the input
3. run the `run_validation_pointsource.sh` to run the forward modeling as true model.
4. run the `run_validation_extsource.sh` to run the forward modeling with external source.
5. run the `plot_comparison.py` to plot the comparison
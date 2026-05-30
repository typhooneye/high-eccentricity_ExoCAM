# ExoCAM SourceMods

This folder contains modified `SourceMods` files for ExoCAM.

## What it does

- improves orbital calculations for high eccentricity
- uses `exo_porb` to set the model year length

## How to use

1. Copy these files into the matching `SourceMods` directories in your ExoCAM case.
2. Set the orbital parameters in `src.share/exoplanet_mod(...).F90`.
3. Rebuild ExoCAM.

## Month split

The calendar still uses 12 months. For `exo_porb = N` days:

- each month gets `floor(N/12)` days
- the first `mod(N,12)` months get 1 extra day

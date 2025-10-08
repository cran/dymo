# dymo 1.0.0

* Added a `NEWS.md` file to track changes to the package.


# dymo 2.0.0

Scaled returns only, stable and invertible.

State centering, removes bias, focuses DMD on true dynamics.

Eigenvalue clipping limits to prevent forecast explosion.

Truncated SVD, keeps only dominant modes, filters numerical noise.

Automatic rank selection, best rank chosen via rolling validation.

Per-feature conformal residuals, valid calibration, no scale mixing.

Joint (path-wise) sampling, preserves temporal coherence in simulated trajectories.

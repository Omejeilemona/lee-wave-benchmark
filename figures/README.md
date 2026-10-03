# Figures

The main figures are shown, with captions, in the top-level `README.md`. This file says which run each file
comes from (run numbers refer to `validation_log.md`).

| File | Run | What it shows |
|---|---|---|
| `summary_offset_vs_grid.png` | runs 3, 4, 8, 10 | Ratio of simulated to linear drag on four grids (f = 0) |
| `fig_w_lee_f0_nz256.png`, `fig_flux_lee_f0_nz256.png` | 3 | f = 0, dx 39 m, dz 3.9 m |
| `fig_w_lee_f0_nz512.png`, `fig_flux_lee_f0_nz512.png` | 4 | f = 0, dx 39 m, dz 1.95 m |
| `fig_w_lee_f0_nx1024_nz256.png`, `fig_flux_lee_f0_nx1024_nz256.png` | 8 | f = 0, dx 19.5 m, dz 3.9 m |
| `fig_w_lee_f0_nx1024_nz512.png`, `fig_flux_lee_f0_nx1024_nz512.png` | 10 | f = 0, dx 19.5 m, dz 1.95 m (finest grid) |
| `fig_w_lee_frot_nx1024_nz256.png`, `fig_flux_lee_frot_nx1024_nz256.png` | 9 | f = 1.4e-4 1/s, dx 19.5 m, dz 3.9 m |
| `fig_w_lee_f0_nz256_h40.png`, `fig_flux_lee_f0_nz256_h40.png` | 7 | f = 0, ridge height 40 m (Fr_L = 0.4), dx 39 m, dz 3.9 m |
| `flux_vs_theory_corrected.png` | 14 | Single-file demo: flux at z = -600 m over 18 h (Centered advection) |
| `internal_wave_generation_corrected_small.gif` | 14 | Single-file demo: vertical velocity over 18 h (reduced size) |

Naming: `fig_w_*` shows the vertical velocity field (top) and a profile against linear theory (bottom);
`fig_flux_*` shows the drag against time and the linear prediction.
`nz256`/`nz512` is the vertical grid (256 or 512 cells over 1 km) and `nx1024` is the horizontal grid (1024 cells
over 20 km); runs without `nx1024` used 512 horizontal cells.

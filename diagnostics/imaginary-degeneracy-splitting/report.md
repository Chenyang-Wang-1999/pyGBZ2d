# imaginary-degeneracy-splitting diagnostic report

Reference energy: `(1.648+0.0294j)`
Curve definition: solve `a2(mu1, mu2)=0`, then evaluate `a1(mu1, mu2)`.

## check_amoeba result

success=True, is_amoeba=False, classification=nearby_zero_plateau, plateau_check=found, mu1=0.11328125, mu2=0.0, is_continuum=False, zero_count=4, net_zero_count=0

## Resolution sensitivity at returned mu1

| N_points | mu2 | a1 | zero_count | net_zero_count | is_continuum |
|---:|---:|---:|---:|---:|---:|
| 301 | 0 | -1.30049346975e-14 | 4 | 0 | 0 |
| 401 | 0 | 0 | 0 | 0 | 0 |
| 601 | 0 | -0.0038120915485 | 4 | 4 | 0 |
| 901 | 0 | -0.00381209154862 | 4 | 4 | 0 |
| 1201 | 0 | -0.00381209154849 | 4 | 4 | 0 |

## a1(mu1) sweep

mu1 range: `[0.09, 0.13]`, samples: `81`, N_points: `901`
Minimum adjacent delta: `0.0`
Monotonic decreases below tolerance 1e-06: `0`

## Zero plateau test

Plateau predicate: `abs(a1) <= 1e-08` and `zero_count == 0`.

| sampled_start | sampled_end | refined_start | refined_end | width |
|---:|---:|---:|---:|---:|
| 0.1135 | 0.116 | 0.113282599424 | 0.11647195474 | 0.00318935531564 |

## Diagnostic conclusion

No monotonicity violation was found, and the zero plateau criterion indicates the energy is outside the amoeba spectrum.

## Files

CSV data: `diagnostics/imaginary-degeneracy-splitting/a1_vs_mu1.csv`

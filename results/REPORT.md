# Reproduced results

- Source: Cotfas, Deaconu & Cotfas (2021), Table 21, https://doi.org/10.1049/rpg2.12262
- Data: 34 published samples at STC; 25 training, 9 withheld for STC interpolation testing.
- STC discrete measured peak: 204.737 W; datasheet nameplate: 200 W.
- NOCT reference: manufacturer nominal Pmax 142 W at 800 W/m², 47°C **cell** temperature.

| Model | Train RMSE (A) | Held-out RMSE (A) | STC predicted Pmax (W) | NOCT predicted Pmax (W) | NOCT absolute Pmax error |
|---|---:|---:|---:|---:|---:|
| SDM | 0.01754 | 0.01050 | 203.65 | 147.12 | 3.61% |
| DDM | 0.01754 | 0.01050 | 203.65 | 147.12 | 3.61% |

Lowest held-out RMSE: **SDM**, 0.01050 A (interpolation within one STC curve).

Environmental sweep: **15 distinct irradiance/temperature settings per model, 30 predictions total**. Only the 800 W/m² / 47°C NOCT **rating** is used as an out-of-condition power comparison.

## Scientific limits

The STC test split is not validation against new weather conditions, and manufacturer NOCT is a specification rather than independent laboratory measurements. Double-diode parameters can be non-identifiable. Report actual results with these caveats.

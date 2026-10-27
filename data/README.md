# Dataset provenance

`kc200gt_stc_cotfas2021_table21.csv` contains 34 measured (voltage, current) numerical pairs as printed in **Table 21** of:

D. T. Cotfas, A. M. Deaconu, and P. A. Cotfas (2021), *Hybrid successive discretisation algorithm used to calculate parameters of the photovoltaic cells and panels for existing datasets*. IET Renewable Power Generation. https://doi.org/10.1049/rpg2.12262

- The research article republishes an earlier KC200GT curve. This repository attributes the accessible Table 21, not the original measurement campaign.
- Do not call it a new experiment or a synthetic curve.
- Metadata: KC200GT, irradiance **1000 W/m²**, cell temperature **25°C**.
- `sample_id` refers to the row order in the source table (1-based).
- Manufacturer nominal STC maximum power (**200 W**) is not identical to the largest discrete power in this published dataset (**204.737 W**). These are different reference descriptions.
- All 34 points are used only for validation plots; fitting uses 25 training rows, with 9 disjoint rows for interpolation testing.
- Manufacturer NOCT reference **142 W**, **800 W/m²**, **47°C cell temperature**: https://www.energymatters.com.au/images/kyocera/KC200GT.pdf

No measured irradiance/temperature sweep curves are claimed or generated here. The environmental sweep CSV is explicitly marked model predictions.

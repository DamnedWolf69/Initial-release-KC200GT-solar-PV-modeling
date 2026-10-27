"""CLI: fit, validate, compare NOCT, sweep conditions, plot, export outputs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .fit import fit_model
from .model import MODULE, Params, current, mpp

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data' / 'kc200gt_stc_cotfas2021_table21.csv'


def error_metrics(actual, predicted):
    a = np.asarray(actual)
    p = np.asarray(predicted)
    error = p - a
    return dict(RMSE_A=float(np.sqrt(np.mean(error**2))),
                MAE_A=float(np.mean(np.abs(error))),
                R2=float(1 - np.sum(error**2) / np.sum((a - np.mean(a))**2)))


def export_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False, float_format='%.12g')


def param_record(kind: str, p: Params):
    return {'Model': kind, 'Iph_A': p.iph_A, 'I01_A': p.i0_A[0],
            'I02_A': p.i0_A[1] if len(p.i0_A) > 1 else None,
            'n1': p.n[0], 'n2': p.n[1] if len(p.n) > 1 else None,
            'Rs_ohm': p.rs_ohm, 'Rsh_ohm': p.rsh_ohm}


def figures(data, labels, models, comparison, sweep, fig_dir):
    fig_dir.mkdir(parents=True, exist_ok=True)
    colors = {'SDM': '#1f77b4', 'DDM': '#d95f02'}
    v = np.linspace(0, 34.0, 400)
    def save(name):
        plt.tight_layout()
        plt.savefig(fig_dir / f'{name}.png', dpi=180)
        plt.savefig(fig_dir / f'{name}.svg')
        plt.close()

    plt.figure(figsize=(9, 5.5))
    for mask, marker, label in ((labels == 'train', 'o', 'Training (25)'),
                                (labels == 'test', 's', 'Held-out (9)')):
        plt.scatter(data[mask, 0], data[mask, 1], label=label, s=42,
                    marker=marker, facecolors='none', edgecolors='#222222', zorder=5)
    for kind, p in models.items():
        plt.plot(v, current(v, p), label=kind, color=colors[kind], linewidth=2)
    plt.xlabel('Voltage (V)'); plt.ylabel('Current (A)')
    plt.title('KC200GT - published STC I-V samples and diode fits')
    plt.xlim(0, 34); plt.ylim(-0.2, 8.7)
    plt.legend(); plt.grid(alpha=.22)
    save('01_stc_iv_fits')

    plt.figure(figsize=(9, 5.5))
    plt.scatter(data[:, 0], data[:, 0]*data[:, 1], s=26,
                label='Published sampled I-V', color='#3e3e3e')
    for kind, p in models.items():
        plt.plot(v, v*current(v, p), label=kind, color=colors[kind], linewidth=2)
    plt.scatter([MODULE.vmp_datasheet_V],
                [MODULE.vmp_datasheet_V * MODULE.imp_datasheet_A],
                label='Datasheet STC reference', marker='x', s=95, color='#a32a2a')
    plt.xlabel('Voltage (V)'); plt.ylabel('Power (W)')
    plt.title('STC P-V curve and maximum-power estimates')
    plt.legend(); plt.grid(alpha=.22)
    save('02_stc_pv_fits')

    plt.figure(figsize=(9, 5.5))
    for kind, p in models.items():
        pred = current(data[:, 0], p)
        mask = labels == 'test'
        plt.scatter(data[mask, 0], pred[mask] - data[mask, 1],
                    label=f'{kind} test residuals', s=55, color=colors[kind])
    plt.axhline(0, color='#222222', lw=1)
    plt.xlabel('Voltage (V)'); plt.ylabel('Prediction minus measured current (A)')
    plt.title('Held-out residuals - STC interpolation')
    plt.legend(); plt.grid(alpha=.22)
    save('03_test_residuals')

    plt.figure(figsize=(9, 5.5))
    p = models['SDM']
    for G in [200,400,600,800,1000]:
        voc=mpp(p,G,25)['Voc_V']
        vx=np.linspace(0,voc,260)
        plt.plot(vx,current(vx,p,G,25),label=f'{G} W/m2',linewidth=1.8)
    plt.xlabel('Voltage (V)'); plt.ylabel('Current (A)')
    plt.title('SDM predicted I-V: varying irradiance at 25 C cell temperature')
    plt.legend(); plt.grid(alpha=.22)
    save('04_irradiance_iv_predictions')

    plt.figure(figsize=(9, 5.5))
    for T in [25,50,75]:
        voc=mpp(p,1000,T)['Voc_V']
        vx=np.linspace(0,voc,260)
        plt.plot(vx,current(vx,p,1000,T),label=f'{T} C',linewidth=1.8)
    plt.xlabel('Voltage (V)'); plt.ylabel('Current (A)')
    plt.title('SDM predicted I-V: varying cell temperature at 1000 W/m2')
    plt.legend(); plt.grid(alpha=.22)
    save('05_temperature_iv_predictions')

    plt.figure(figsize=(9, 5.5))
    df=pd.DataFrame(sweep)
    for kind in ['SDM','DDM']:
        for T in [25,50,75]:
            sub=df[(df['Model'] == kind) & (df['Tcell_C']==T)].sort_values('G_Wm2')
            plt.plot(sub['G_Wm2'],sub['Pmax_W'],
                     linestyle='-' if kind=='SDM' else '--', marker='o',
                     color=['#2268ad','#00917a','#a65e2c'][[25,50,75].index(T)],
                     label=f'{kind}, {T} C', linewidth=1.6)
    plt.xlabel('Irradiance (W/m2)'); plt.ylabel('Predicted Pmax (W)')
    plt.title('Power sensitivity - 15 conditions x 2 models (predictions)')
    plt.legend(ncol=2); plt.grid(alpha=.22)
    save('06_pmax_sweep_predictions')

    plt.figure(figsize=(7.5, 5))
    kinds=['SDM','DDM']
    errs=[next(r['NOCT_Pmax_abs_error_pct'] for r in comparison if r['Model']==k) for k in kinds]
    bars=plt.bar(kinds, errs, color=[colors[k] for k in kinds])
    plt.bar_label(bars, fmt='%.2f%%', padding=3)
    plt.ylabel('Absolute Pmax error (%)')
    plt.title('NOCT comparison against 142 W datasheet rating')
    plt.grid(axis='y',alpha=.22)
    save('07_noct_rating_error')


def run(output_dir: Path, figures_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA).to_numpy(dtype=float)
    assert data.shape == (34,2)
    n = len(data)
    # Fixed, predetermined split, sample indices 1,5,9,...,33 (1-based)
    mask = np.arange(n)%4 == 0
    labels = np.where(mask,'test','train')
    manifest = [dict(sample_id=i+1,voltage_V=float(v),current_A=float(y),subset=str(labels[i]))
                for i,(v,y) in enumerate(data)]
    export_csv(output_dir / 'split_manifest.csv',manifest)
    models = {}; optimization={}; metrics=[]; parameters=[]; per_point=[]
    for kind in ('SDM','DDM'):
        p, fitinfo=fit_model(kind,data[~mask,0],data[~mask,1])
        models[kind]=p; optimization[kind]=fitinfo
        parameters.append(param_record(kind,p))
        pred=current(data[:,0],p)
        tr=error_metrics(data[~mask,1],pred[~mask])
        te=error_metrics(data[mask,1],pred[mask])
        stc=mpp(p)
        noct=mpp(p, MODULE.noct_irradiance, MODULE.noct_cell_temp_C)
        row=dict(Model=kind, Train_RMSE_A=tr['RMSE_A'],
                 Train_MAE_A=tr['MAE_A'], Test_RMSE_A=te['RMSE_A'],
                 Test_MAE_A=te['MAE_A'], Test_R2=te['R2'],
                 STC_Pmax_W=stc['Pmax_W'], STC_Vmp_V=stc['Vmp_V'],
                 STC_Imp_A=stc['Imp_A'], STC_Voc_V=stc['Voc_V'],
                 STC_FF_pct=stc['FF']*100,
                 STC_vs_datasheet_Pmax_pct=100*(stc['Pmax_W']-MODULE.vmp_datasheet_V*MODULE.imp_datasheet_A)/(MODULE.vmp_datasheet_V*MODULE.imp_datasheet_A),
                 STC_vs_measured_sparse_peak_pct=100*(stc['Pmax_W']-float(max(data[:,0]*data[:,1])))/float(max(data[:,0]*data[:,1])),
                 NOCT_Pmax_W=noct['Pmax_W'],NOCT_Vmp_V=noct['Vmp_V'],
                 NOCT_Imp_A=noct['Imp_A'],NOCT_Voc_V=noct['Voc_V'],
                 NOCT_Pmax_signed_error_pct=100*(noct['Pmax_W']-MODULE.noct_pmax_datasheet_W)/MODULE.noct_pmax_datasheet_W,
                 NOCT_Pmax_abs_error_pct=100*abs(noct['Pmax_W']-MODULE.noct_pmax_datasheet_W)/MODULE.noct_pmax_datasheet_W,
                 NOCT_Vmp_signed_error_pct=100*(noct['Vmp_V']-MODULE.noct_vmp_datasheet_V)/MODULE.noct_vmp_datasheet_V,
                 NOCT_Imp_signed_error_pct=100*(noct['Imp_A']-MODULE.noct_imp_datasheet_A)/MODULE.noct_imp_datasheet_A)
        metrics.append(row)
        for i, (voltage, measured) in enumerate(data):
            per_point.append(dict(Model=kind,sample_id=i+1,subset=labels[i],
                                  V_V=voltage, measured_I_A=measured,
                                  predicted_I_A=pred[i], residual_A=pred[i]-measured))

    sweep=[]
    for kind,p in models.items():
        for t in (25,50,75):
            for G in (200,400,600,800,1000):
                result=mpp(p,G,t)
                sweep.append(dict(Model=kind,G_Wm2=G,Tcell_C=t,**result))

    export_csv(output_dir/'model_comparison.csv',metrics)
    export_csv(output_dir/'fitted_parameters.csv',parameters)
    export_csv(output_dir/'pointwise_predictions.csv',per_point)
    export_csv(output_dir/'environment_sweep.csv',sweep)
    references={'data':'Cotfas, Deaconu & Cotfas (2021), Table 21, https://doi.org/10.1049/rpg2.12262',
                'manufacturer':'Kyocera KC200GT datasheet, https://www.energymatters.com.au/images/kyocera/KC200GT.pdf',
                'temperature_scaling':'PVPMC De Soto module model, https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/single-diode-equivalent-circuit-models/de-soto-five-parameter-module-model/'}
    record=dict(description='Reproducible two-model fit: 25 training, 9 interpolation-test observations',
                counts={'published_measured_points':34,'training_points':25,
                        'held_out_points':9,'environment_conditions_per_model':15},
                nominal_references={'STC_Pmax_W':200.0,
                                    'STC_datasheet_Vmp_times_Imp_W':MODULE.vmp_datasheet_V*MODULE.imp_datasheet_A,
                                    'NOCT_Pmax_rating_W':142.0,
                                    'NOCT_Vmp_times_Imp_W':MODULE.noct_vmp_datasheet_V*MODULE.noct_imp_datasheet_A,
                                    'measured_sparse_peak_W':float(max(data[:,0]*data[:,1]))},
                fitting=optimization, comparison=metrics, parameters=parameters,
                notes=['NOCT is an out-of-condition comparison to a manufacturer rating, NOT lab I-V validation.',
                       'The held-out STC set tests interpolation within a single measured curve, NOT unseen weather conditions.',
                       'Extra model parameters are not uniquely identifiable; check bounds and sensitivity.',
                       'Temperature correction is De Soto-inspired and extended to both DDM diodes under assumptions.',
                       'No temperature/irradiance swept I-V measurements are fabricated.'])
    (output_dir/'metrics.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    figures(data, labels, models, metrics, sweep, figures_dir)
    best=min(metrics,key=lambda row:row['Test_RMSE_A'])
    summary=(f"# Reproduced results\n\n"
             f"- Source: {references['data']}\n"
             f"- Data: 34 published samples at STC; 25 training, 9 withheld for STC interpolation testing.\n"
             f"- STC discrete measured peak: {max(data[:,0]*data[:,1]):.3f} W; datasheet nameplate: 200 W.\n"
             f"- NOCT reference: manufacturer nominal Pmax 142 W at 800 W/m², 47°C **cell** temperature.\n\n"
             f"| Model | Train RMSE (A) | Held-out RMSE (A) | STC predicted Pmax (W) | NOCT predicted Pmax (W) | NOCT absolute Pmax error |\n"
             f"|---|---:|---:|---:|---:|---:|\n")
    for row in metrics:
        summary+=(f"| {row['Model']} | {row['Train_RMSE_A']:.5f} | {row['Test_RMSE_A']:.5f} | "
                  f"{row['STC_Pmax_W']:.2f} | {row['NOCT_Pmax_W']:.2f} | "
                  f"{row['NOCT_Pmax_abs_error_pct']:.2f}% |\n")
    summary+=(f"\nLowest held-out RMSE: **{best['Model']}**, "
              f"{best['Test_RMSE_A']:.5f} A (interpolation within one STC curve).\n\n"
              f"Environmental sweep: **15 distinct irradiance/temperature settings per model, 30 predictions total**. "
              f"Only the 800 W/m² / 47°C NOCT **rating** is used as an out-of-condition power comparison.\n\n"
              f"## Scientific limits\n\nThe STC test split is not validation against new weather conditions, "
              f"and manufacturer NOCT is a specification rather than independent laboratory measurements. "
              f"Double-diode parameters can be non-identifiable. Report actual results with these caveats.\n")
    (output_dir/'REPORT.md').write_text(summary)
    print(summary)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    parser.add_argument('--figures',type=Path,default=ROOT/'figures')
    args=parser.parse_args()
    run(args.output,args.figures)

if __name__=='__main__':
    main()

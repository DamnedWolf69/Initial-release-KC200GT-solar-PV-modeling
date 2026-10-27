import numpy as np
import pandas as pd
import pytest
from kc200gt.model import MODULE, current, mpp, decode
from kc200gt.analysis import DATA, error_metrics


SDM = decode(np.array([8.19,-8.5,1.3,0.28,2.2]),'SDM')
DDM = decode(np.array([8.19,-9.0,-8.5,1.1,2.0,0.28,2.2]),'DDM')


@pytest.mark.parametrize('params',[SDM,DDM])
def test_current_descends_and_oc(params):
    peak = mpp(params)
    volt = np.linspace(0,peak['Voc_V'],100)
    iv = current(volt,params)
    assert np.all(np.diff(iv) <= 1e-6)
    assert abs(iv[-1]) < 1e-7
    assert 0 < peak['Vmp_V'] < peak['Voc_V']
    assert peak['Pmax_W'] == pytest.approx(peak['Vmp_V']*peak['Imp_A'],rel=1e-11)
    assert 0 < peak['FF'] < 1


@pytest.mark.parametrize('params',[SDM,DDM])
def test_physical_irradiance_and_temperature_trends(params):
    assert mpp(params,1000,25)['Pmax_W'] > mpp(params,600,25)['Pmax_W']
    assert mpp(params,1000,75)['Pmax_W'] < mpp(params,1000,25)['Pmax_W']


def test_public_data_and_split():
    dataset=pd.read_csv(DATA)
    assert dataset.shape==(34,2)
    assert dataset.voltage_V.is_monotonic_increasing
    assert (np.arange(34)%4==0).sum()==9
    assert np.max(dataset.voltage_V.to_numpy()*dataset.current_A.to_numpy()) == pytest.approx(204.73702938)


def test_rmse_contract():
    m=error_metrics(np.array([1.,2.,3.]),np.array([1.,2.,3.]))
    assert m['RMSE_A']==0
    assert m['MAE_A']==0
    assert m['R2']==pytest.approx(1)


def test_noct_rating_separate_from_voltage_times_current():
    assert MODULE.noct_pmax_datasheet_W==142.0
    assert MODULE.noct_vmp_datasheet_V*MODULE.noct_imp_datasheet_A == pytest.approx(142.216)


def test_fitted_sdm_regression():
    from kc200gt.fit import fit_model
    data = pd.read_csv(DATA).to_numpy(dtype=float)
    split = np.arange(len(data)) % 4 == 0
    params, info = fit_model('SDM', data[~split, 0], data[~split, 1])
    preds = current(data[split, 0], params)
    metrics = error_metrics(data[split, 1], preds)
    noct = mpp(params,800,47)
    assert info['n_starts'] == 3
    assert metrics['RMSE_A'] == pytest.approx(0.01050, abs=0.00015)
    assert noct['Pmax_W'] == pytest.approx(147.12, abs=0.20)

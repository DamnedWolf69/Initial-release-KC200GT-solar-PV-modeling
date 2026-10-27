% KC200GT reproducible SDM/DDM fitting + STC and NOCT checks.
% Requires MATLAB R2019b+ and Optimization Toolbox (lsqnonlin).
% Equivalent workflow to the tested Python runner in this repository.
% IMPORTANT: Python run was executed and tested; MATLAB is provided for
% native MATLAB replication and was not executable in the build environment.
%
% DATA: Cotfas et al., IET RPG 2021, Table 21, DOI 10.1049/rpg2.12262
% NOCT: Kyocera KC200GT datasheet, 800 W/m2, cell temperature 47 C.
% The NOCT 142 W value is a rounded rating, not independent measured I-V data.
clear; clc; close all;
repo = fileparts(fileparts(mfilename('fullpath')));
data = readmatrix(fullfile(repo,'data','kc200gt_stc_cotfas2021_table21.csv'));
assert(isequal(size(data),[34 2]),'Expected 34 sourced points.');
V = data(:,1); I = data(:,2);
% Same deterministic indices as Python: 1, 5, 9, ..., 33 are test rows.
test = false(34,1); test(1:4:end) = true;
train = ~test;

spec.Ns = 54; spec.k = 1.380649e-23; spec.q = 1.602176634e-19;
spec.alphaIsc = 0.00318; spec.Eg_eV = 1.12;
spec.TrefC = 25; spec.Gref = 1000;
spec.VmpDS = 26.3; spec.ImpDS = 7.61;
spec.PnoctDS = 142; spec.Gnoct = 800; spec.TnoctC = 47;

% Fitting coordinates: SDM [Iph, log10(I01), n1, Rs, log10(Rsh)]
% DDM [Iph, log10(I01), log10(I02), n1, n2, Rs, log10(Rsh)]
% A physically restricted DDM may collapse to a SDM; do not conceal this.
seedsSDM = [8.19 -9 1.1 .20 2.2; 8.19 -8.5 1.3 .30 2.3; 8.19 -7 1.7 .20 2.0];
seedsDDM = [8.19 -9 -9 1.05 2 .20 2.2; 8.19 -10 -8 1.1 2.2 .25 2.2; ...
            8.19 -8.8 -11 1.1 2.5 .30 2.0; 8.19 -9 -10 1.3 1.6 .25 2.4];
LB1=[7.8 -13 1 0 1]; UB1=[8.6 -4 2.5 .8 4];
LB2=[7.8 -13 -13 1 1.5 0 1]; UB2=[8.6 -4 -4 1.5 3 .8 4];
opts = optimoptions('lsqnonlin','Display','off','MaxIterations',1000,...
     'MaxFunctionEvaluations',16000,'FunctionTolerance',1e-11,...
     'StepTolerance',1e-11,'OptimalityTolerance',1e-11);

names = {'SDM','DDM'};
coords = cell(1,2);
coords{1} = multiStart('SDM',seedsSDM,LB1,UB1,V(train),I(train),spec,opts);
coords{2} = multiStart('DDM',seedsDDM,LB2,UB2,V(train),I(train),spec,opts);

output=fullfile(repo,'matlab','results_local');
if ~isfolder(output), mkdir(output); end
% Explicitly report both absolute error and signed error against NOCT rating.
T = table(); paramTable = table(); predictions = table(); sweep = table();
for j = 1:2
    kind = names{j}; x=coords{j};
    predicted = solarCurrent(V,x,kind,1000,25,spec);
    errTrain = sqrt(mean((predicted(train)-I(train)).^2));
    errTest = sqrt(mean((predicted(test)-I(test)).^2));
    maeTest = mean(abs(predicted(test)-I(test)));
    [Pstc,Vstc,Istc,Vocstc] = solarMPP(x,kind,1000,25,spec);
    [Pnoct,Vnoct,Inoct,Vocnoct] = solarMPP(x,kind,800,47,spec);
    NOCTpctSigned = 100*(Pnoct-spec.PnoctDS)/spec.PnoctDS;
    row = table(string(kind),errTrain,errTest,maeTest,...
        Pstc,Vstc,Istc,Vocstc,Pnoct,Vnoct,Inoct,Vocnoct,...
        NOCTpctSigned,abs(NOCTpctSigned),...
        'VariableNames',{'Model','Train_RMSE_A','Test_RMSE_A',...
        'Test_MAE_A','STC_Pmax_W','STC_Vmp_V','STC_Imp_A','STC_Voc_V',...
        'NOCT_Pmax_W','NOCT_Vmp_V','NOCT_Imp_A','NOCT_Voc_V',...
        'NOCT_Pmax_signed_error_pct','NOCT_Pmax_abs_error_pct'});
    T=[T; row]; %#ok<AGROW>
    if strcmp(kind,'SDM')
       pvals=[x(1),10^x(2),NaN,x(3),NaN,x(4),10^x(5)];
    else
       pvals=[x(1),10^x(2),10^x(3),x(4),x(5),x(6),10^x(7)];
    end
    parameterRow = array2table(pvals,'VariableNames',...
       {'Iph_A','I01_A','I02_A','n1','n2','Rs_ohm','Rsh_ohm'});
    parameterRow = addvars(parameterRow,string(kind),'Before',1,'NewVariableNames','Model');
    paramTable=[paramTable;parameterRow]; %#ok<AGROW>
    for i=1:numel(V)
        subset="train"; if test(i),subset="test";end
        one=table(string(kind),i,subset,V(i),I(i),predicted(i),...
             predicted(i)-I(i),'VariableNames',{'Model','sample_id','subset',...
             'V_V','measured_I_A','predicted_I_A','residual_A'});
        predictions=[predictions;one]; %#ok<AGROW>
    end
    for Tc=[25 50 75]
        for G=[200 400 600 800 1000]
            [Pm,Vm,Im,Voc]=solarMPP(x,kind,G,Tc,spec);
            sw=table(string(kind),G,Tc,Pm,Vm,Im,Voc,'VariableNames',...
              {'Model','G_Wm2','Tcell_C','Pmax_W','Vmp_V','Imp_A','Voc_V'});
            sweep=[sweep;sw]; %#ok<AGROW>
        end
    end
end
writetable(T,fullfile(output,'model_comparison_matlab.csv'));
writetable(paramTable,fullfile(output,'fitted_parameters_matlab.csv'));
writetable(predictions,fullfile(output,'pointwise_predictions_matlab.csv'));
writetable(sweep,fullfile(output,'environment_sweep_matlab.csv'));
disp('Measured points: 34; fit: 25; held out: 9');
disp(T);
fprintf('Measured sparse peak=%.6f W; nominal STC=%g W; NOCT rating=%g W\n',...
      max(V.*I),200,spec.PnoctDS);

vplot=linspace(0,34,400)';
figure('Color','w');hold on;
plot(V(train),I(train),'ko','DisplayName','train');
plot(V(test),I(test),'rs','DisplayName','test');
for j=1:2
    plot(vplot,solarCurrent(vplot,coords{j},names{j},1000,25,spec),...
        'LineWidth',1.8,'DisplayName',names{j});
end
xlabel('Voltage (V)');ylabel('Current (A)');
title('KC200GT - published STC I-V, model fits, train/test');
grid on;legend('Location','southwest');
saveas(gcf,fullfile(output,'01_STC_IV_MATLAB.png'));

figure('Color','w');hold on;
plot(V,V.*I,'ko','DisplayName','measured sampled power');
for j=1:2
    curve=solarCurrent(vplot,coords{j},names{j},1000,25,spec);
    plot(vplot,vplot.*curve,'LineWidth',1.8,'DisplayName',names{j});
end
xlabel('Voltage (V)');ylabel('Power (W)');grid on;
title('STC P-V validation');legend('Location','northwest');
saveas(gcf,fullfile(output,'02_STC_PV_MATLAB.png'));

figure('Color','w');hold on;
for j=1:2
    pred=solarCurrent(V(test),coords{j},names{j},1000,25,spec);
    plot(V(test),pred-I(test),'-o','DisplayName',names{j});
end
yline(0,'k-');xlabel('Voltage (V)');ylabel('I residual (A)');
title('Held-out STC residuals');grid on;legend();
saveas(gcf,fullfile(output,'03_Heldout_residuals_MATLAB.png'));

figure('Color','w');hold on;
for Tc=[25 50 75]
    rows = strcmp(sweep.Model,'SDM') & sweep.Tcell_C==Tc;
    plot(sweep.G_Wm2(rows),sweep.Pmax_W(rows),'-o','LineWidth',1.8,...
        'DisplayName',sprintf('%d C',Tc));
end
xlabel('Irradiance (W/m2)');ylabel('Predicted maximum power (W)');
title('SDM predicted irradiance and temperature sensitivity');
grid on;legend('Location','northwest');
saveas(gcf,fullfile(output,'04_Pmax_sweep_MATLAB.png'));
fprintf('Wrote MATLAB tables and PNG figures to %s\n',output);

function best=multiStart(kind,seeds,lb,ub,V,I,spec,opts)
    bestSSE=Inf;best=[];
    for si=1:size(seeds,1)
        fun=@(x) solarCurrent(V,x,kind,1000,25,spec)-I;
        [xc,sse]=lsqnonlin(fun,seeds(si,:),lb,ub,opts);
        if sse<bestSSE,bestSSE=sse;best=xc;end
    end
    fprintf('%s best training SSE: %.10g A^2\n',kind,bestSSE);
end

function p=adjusted(x,kind,G,Tc,spec)
    if strcmp(kind,'SDM')
        iph=x(1);i0=10.^x(2);nn=x(3);rs=x(4);rsh=10.^x(5);
    else
        iph=x(1);i0=10.^x(2:3);nn=x(4:5);rs=x(6);rsh=10.^x(7);
    end
    T=Tc+273.15; T0=spec.TrefC+273.15;
    egT=spec.Eg_eV*(1-0.0002677*(Tc-spec.TrefC));
    p.iph=G/spec.Gref*(iph+spec.alphaIsc*(Tc-spec.TrefC));
    p.i0=i0.*(T/T0)^3.*exp((spec.Eg_eV/T0-egT/T)./(nn*(spec.k/spec.q)));
    p.nn=nn;p.rs=rs;p.rsh=rsh*spec.Gref/G;
    p.a=nn*spec.Ns*spec.k*T/spec.q;
end

function I=solarCurrent(V,x,kind,G,Tc,spec)
    p=adjusted(x,kind,G,Tc,spec);
    sz=size(V);V=V(:); I=p.iph-V/p.rsh;
    for iter=1:90
        vd=V+I*p.rs;
        F=p.iph-I-vd/p.rsh;
        dF=-(1+p.rs/p.rsh)*ones(size(I));
        for d=1:numel(p.i0)
            z=max(-100,min(100,vd/p.a(d)));
            ez=exp(z);
            F=F-p.i0(d)*(ez-1);
            dF=dF-p.i0(d)*ez*p.rs/p.a(d);
        end
        dx=max(-4,min(4,F./dF));
        I=I-dx;
        if max(abs(dx))<1e-12,break;end
    end
    I=reshape(I,sz);
end

function [Pm,Vm,Im,Voc]=solarMPP(x,kind,G,Tc,spec)
    p=adjusted(x,kind,G,Tc,spec);
    Voc=fzero(@(v) openCircuitResidual(v,p),[0,80]);
    opts=optimset('TolX',1e-9);
    [Vm,negativePm]=fminbnd(@(v)-v*solarCurrent(v,x,kind,G,Tc,spec),0,Voc,opts);
    Pm=-negativePm;
    Im=Pm/Vm;
end

function f=openCircuitResidual(v,p)
    f=p.iph-v/p.rsh;
    for d=1:numel(p.i0)
        f=f-p.i0(d)*(exp(max(-100,min(100,v/p.a(d))))-1);
    end
end

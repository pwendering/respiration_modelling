% Rank reactions of in a day- and a night-specific version of the AraCore
% model by their flux

clear;clc

initCobraToolbox(false)

top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'fun_models.mat'));
[models, model_publications, ~] = deal(tmp.fun_models, tmp.fun_model_publications, tmp.fun_organisms_uniq);
model_orig = models{ismember(model_publications, 'ArnoldNikoloski2014')};
model = model_orig;
clear models model_publications tmp

% find CO2 ID in model
co2_form_match = ismember(lower(model.metFormulas), 'co2');
co2_ids = model.mets(co2_form_match);

%% find carboxylation and oxygenation reactions
v_c = {'RBC_h'};
v_o = {'RBO_h'};

%% constrain oxygenation to carboxylation ratio
s_co = [
    92 107 89 105,... % Parry et al. (1989) [https://doi.org/10.1093/jxb/40.3.317]
    82 74 89 93 61 66,... % Zhu et al. (1992) [https://doi.org/10.1104%2Fpp.98.2.764]
    99.9 100.8 96.2 96 98.7 98.4 97 99.2 98.5 94 100.8 95.6 93.1 99.7,...
    92.4 95.4 97.0 100.1 97.5 82.2 87.3 84.4]; % Hermida-Carrera et al. (2016) [https://doi.org/10.1104%2Fpp.16.01846]
f_mol_bar = 2417.20 / 104.90; % Walker et al. (2013) [https://doi.org/10.1111/pce.12166]
O = 210000;  % Farquhar et al. (1980), Planta
C = 230;  % Farquhar et al. (1980), Planta
% phi = 0.27; % Farquhar et al. (1980), Planta

phi_array = (1 ./ (f_mol_bar*s_co))*(O/C);
phi = mean(phi_array);
phi_tol = 2 * std(phi_array);


% fix ratio between RUBISCO carbocylation and oxygenation within a
% tolerated range
model = addMetabolite(model, 'phi_lb',...
    'csense', 'L');
model.S(findMetIDs(model, 'phi_lb'), findRxnIDs(model, [v_c v_o])) = [phi-phi_tol -1];

model = addMetabolite(model, 'phi_ub',...
    'csense', 'L');
model.S(findMetIDs(model, 'phi_ub'), findRxnIDs(model, [v_c v_o])) = [-phi-phi_tol 1];

% obtain carbon molar fraction in biomass for scaling
cmf = getCarbonMolarFraction(model, true);

%% Solve FBA and extract flux sum through CO2 producing (and consuming) reactions
% COBRA solver parameters
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'timeLimit', 180);

%% day model

% load day modifications (Arnold et al. 2015)
day_tab = readtable(fullfile('..', '..', 'models',...
    'ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'day_modifications_Arnold2015.txt'));
model_day = changeRxnBounds(model, day_tab.(1), day_tab.(2), 'l');
model_day = changeRxnBounds(model_day, day_tab.(1), day_tab.(3), 'u');


% split reversible reactions
model_day = convertToIrreversible(model_day);

% find indices of reactions that produce or consume CO2
co2_prod_idx = any(model_day.S(findMetIDs(model_day, co2_ids), :) > 0, 1);

% solve FBA problem, minimizing the sum of fluxes and without loops
fba_sol = optimizeCbModel(model_day, 'max', 'one', false);
model_flux = fba_sol.x;
model_flux(model_flux<1e-9) = 0;
day_rgr_opt = fba_sol.f;

% rank CO2-producing reactions by their flux values
tab = cell2table([...
    model_day.rxns(co2_prod_idx)...
    model_day.rxnNames(co2_prod_idx)...
    model_day.rxnECNumbers(co2_prod_idx)...
    [model_day.subSystems{co2_prod_idx}]'...
    printRxnFormula(model_day, model_day.rxns(co2_prod_idx), 0)...
    num2cell(model_flux(co2_prod_idx))...
    num2cell(model_flux(co2_prod_idx)/cmf/day_rgr_opt)...
    ],...
    'VariableNames', {'rxn_id', 'rxn_name', 'rxn_ec' 'rxn_subsystem',...
    'rxn_formula', 'pfba_flux', 'flux_per_int_carbon'});
tab = sortrows(tab, 'pfba_flux', 'descend');

writetable(tab, fullfile('..', '..', 'data', 'ranked_co2_prod_rxns_aracore_day.xlsx'))


%% night model
% load night modifications (Arnold et al. 2015)
night_tab = readtable(fullfile('..', '..', 'models',...
    'ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'night_modifications_Arnold2015.txt'));
model_night = changeRxnBounds(model_orig, night_tab.(1), night_tab.(2), 'l');
model_night = changeRxnBounds(model_night, night_tab.(1), night_tab.(3), 'u');
% allow for starch and glycine exchange
model_night = addExchangeRxn(model_night, 'starch5[h]', -1000, 1000);
model_night = changeRxnBounds(model_night, 'Ex_Gly_h', -1000, 'l');
% disable glucose uptake
model_night = changeRxnBounds(model_night, 'Ex_Glc', 0, 'b');

% split reversible reactions
model_night = convertToIrreversible(model_night);

% limit rgr by day growth rate
model_night.ub(findRxnIDs(model_night, 'Bio_opt')) = day_rgr_opt;

% find indices of reactions that produce or consume CO2
co2_prod_idx = any(model_night.S(findMetIDs(model_night, co2_ids), :) > 0, 1);

% solve FBA problem, minimizing the sum of fluxes and without loops
fba_sol = optimizeCbModel(model_night, 'max', 'one', false);
model_flux = fba_sol.x;
model_flux(model_flux<1e-9) = 0;

% rank CO2-producing reactions by their flux values
tab = cell2table([...
    model_night.rxns(co2_prod_idx)...
    model_night.rxnNames(co2_prod_idx)...
    model_night.rxnECNumbers(co2_prod_idx)...
    [model_night.subSystems{co2_prod_idx}]'...
    printRxnFormula(model_night, model_night.rxns(co2_prod_idx), 0)...
    num2cell(model_flux(co2_prod_idx))...
    num2cell(model_flux(co2_prod_idx)/cmf/day_rgr_opt)...
    ],...
    'VariableNames', {'rxn_id', 'rxn_name', 'rxn_ec', 'rxn_subsystem',...
    'rxn_formula', 'pfba_flux', 'flux_per_int_carbon'});
tab = sortrows(tab, 'pfba_flux', 'descend');

writetable(tab, fullfile('..', '..', 'data', 'ranked_co2_prod_rxns_aracore_night.xlsx'))
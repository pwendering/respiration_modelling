% Predict growth respiration during day and night by introducing
% appropriate adaptations into the AraCore model

clear;clc

initCobraToolbox(false)

% load AraCore model
top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'fun_models.mat'));
[models, model_publications, ~] = deal(tmp.fun_models, tmp.fun_model_publications, tmp.fun_organisms_uniq);
model = models{ismember(model_publications, 'ArnoldNikoloski2014')};
clear models model_publications tmp


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

% load night modifications (Arnold et al. 2015)
night_tab = readtable(fullfile('..', '..', 'models',...
    'ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'night_modifications_Arnold2015.txt'));
model_night = changeRxnBounds(model, night_tab.(1), night_tab.(2), 'l');
model_night = changeRxnBounds(model_night, night_tab.(1), night_tab.(3), 'u');
% allow for starch and glycine exchange
model_night = addExchangeRxn(model_night, 'starch5[h]', -1000, 1000);
model_night = changeRxnBounds(model_night, 'Ex_Gly_h', -1000, 'l');
% disable glucose uptake
model_night = changeRxnBounds(model_night, 'Ex_Glc', 0, 'b');

% load day modifications (Arnold et al. 2015)
day_tab = readtable(fullfile('..', '..', 'models',...
    'ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'day_modifications_Arnold2015.txt'));
model_day = changeRxnBounds(model, day_tab.(1), day_tab.(2), 'l');
model_day = changeRxnBounds(model_day, day_tab.(1), day_tab.(3), 'u');

% find CO2 ID in model
co2_form_match = ismember(lower(model.metFormulas), 'co2');
co2_ids = model.mets(co2_form_match);

%% Solve FBA and extract flux sum through CO2 producing (and consuming) reactions
% COBRA solver parameters
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'timeLimit', 180);

% day
model_day = convertToIrreversible(model_day);
co2_prod_idx = find(any(model_day.S(findMetIDs(model_day, co2_ids), :) > 0, 1));
co2_prod_idx = setdiff(co2_prod_idx, find(startsWith(model_day.rxns, {'Tr_', 'Im_'})));
co2_prod_rxns_day = model_day.rxns(co2_prod_idx);

fba_sol = optimizeCbModel(model_day, 'max', 'one', false);
day_rgr = fba_sol.f;
fba_sol.x(fba_sol.x<1e-9) = 0;
c_rxn_flux_day = fba_sol.x(co2_prod_idx);
day_growth_respiration = sum(c_rxn_flux_day);
day_growth_respiration_scaled = day_growth_respiration / fba_sol.f / cmf;

% night 
model_night = convertToIrreversible(model_night);
model_night.ub(findRxnIDs(model_night, 'Bio_opt')) = day_rgr;
co2_prod_idx = find(any(model_night.S(findMetIDs(model_night, co2_ids), :) > 0, 1));
co2_prod_idx = setdiff(co2_prod_idx, find(startsWith(model_night.rxns, {'Tr_', 'Im_'})));
co2_prod_rxns_night = model_night.rxns(co2_prod_idx);

fba_sol = optimizeCbModel(model_night, 'max', 'one', false);

fba_sol.x(fba_sol.x<1e-9) = 0;
c_rxn_flux_night = fba_sol.x(co2_prod_idx);
night_growth_respiration = sum(c_rxn_flux_night);
night_growth_respiration_scaled = night_growth_respiration / fba_sol.f / cmf;

%% Plot fluxes through CO2-producing reactions in day and night
[plot_rxns,ia,ib] = intersect(co2_prod_rxns_day, co2_prod_rxns_night);

names = model_day.rxnNames(findRxnIDs(model_day, plot_rxns));
for i = 1:numel(names)
    res = regexp(erase(plot_rxns{i}, {'_b','_f'}), '_(?<c>\w)$', 'names');
    cof = regexp(erase(plot_rxns{i}, {'_b','_f'}), '(?<c>NAD[P]*)', 'names');
    
    if ~isempty(cof)
        names{i} = strcat(names{i}, [' ' cof.c]);
    end
    
    if isempty(res)
        met_comps = unique(getCompartment(findMetsFromRxns(model_day, plot_rxns{i})));
        if numel(met_comps) > 1
            met_comps = setdiff(met_comps, 'c');
        end
        names{i} = strcat(names{i}, [' (' char(met_comps) ')']);
    else
        names{i} = strcat(names{i}, [' (' res.c ')']);
    end
end

bar(sqrt([c_rxn_flux_day(ia) c_rxn_flux_night(ib)]))
xticks(1:numel(plot_rxns))
xticklabels(names)
xtickangle(90)
ylabel('Flux [mmol/gDW/h]')
legend({'day', 'night'}, 'Box', 'off')
set(gca,...
    'Box', 'on',...
    'LineWidth', 1.3,...
    'FontSize', 16,...
    'FontName', 'Arial',...
    'TickLength', [0 0])
set(gcf, 'OuterPosition', 1000*[0.0830    0.0957    1.1120    0.6253])
exportgraphics(gca, fullfile('..', '..', 'figures', 'day_night_flux_diff_bar.png'),...
    'Resolution', 400)
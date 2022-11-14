% Minimize the ratio between the sum of fluxes through CO2-producing
% reactions and the product of relative growth rate and the molar fraction
% of carbon in the biomass reaction

clear;clc

initCobraToolbox(false)

top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'fun_models.mat'));
[models, model_publications, ~] = deal(tmp.fun_models, tmp.fun_model_publications, tmp.fun_organisms_uniq);
model = models{ismember(model_publications, 'ArnoldNikoloski2014')};
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

%% Solve FBA and extract flux sum through CO2 producing (and consuming) reactions
% COBRA solver parameters
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'timeLimit', 180);

% split reversible reactions
model = convertToIrreversible(model);

% find indices of reactions that produce or consume CO2
co2_prod_idx = find(any(model.S(findMetIDs(model, co2_ids), :) > 0, 1));
co2_cons_idx = find(any(model.S(findMetIDs(model, co2_ids), :) < 0, 1));

% solve FBA to obtain optimal relative growth rate
fba_sol = optimizeCbModel(model, 'max', 'one', true);
rgr_opt = fba_sol.f;
resp_no_opt = sum(fba_sol.x(co2_prod_idx));
opt_distr = fba_sol.x;
opt_distr(opt_distr<1e-9) = 0;
% add new objective that minimizes growth respiration
model_resp = changeObjective(model, model.rxns(co2_prod_idx), -1);
% carbon molar fraction
cmf = getCarbonMolarFraction(model, true);
% Charnes-Cooper transformation
model_cc = addMetabolite(model_resp, 't_const');
model_cc = addReaction(model_cc, 't', 'reactionFormula', '-> t_const');
model_cc.S(end,end) = 0;
model_cc.S(findMetIDs(model_cc, 't_const'), findRxnIDs(model_cc, 'Bio_opt')) = cmf;
model_cc.b(end) = 1;
model_cc.ub(:) = Inf;
% update upper bounds
for i=1:numel(model.rxns)-1
    model_cc = addMetabolite(model_cc, ['t_const_' model.rxns{i}],...
        'csense', 'L');
    model_cc.S(end, findRxnIDs(model_cc, model.rxns{i})) = 1;
    model_cc.S(end, findRxnIDs(model_cc, 't')) = -1000;
end
% solve LFP
cc_sol = optimizeCbModel(model_cc, 'max', 'one');
% calculate flux values
fluxes_resp_red = cc_sol.x / cc_sol.x(end);
% calculate growth respiration and determine relative growth rate
growth_resp_red = sum(fluxes_resp_red(co2_prod_idx));
rgr_resp_red = fluxes_resp_red(findRxnIDs(model_cc, 'Bio_opt'));

% distance of fluxes between distributions resulting from minimization and
% parsimonious FBA
d = fluxes_resp_red(1:numel(model.rxns)) - opt_distr;

% reduce reactions to plot to the ones that are in the top 5 percent of
% absolute flux distances
plot_idx = abs(d)>=prctile(abs(d),95);

plot_mat = [opt_distr fluxes_resp_red(1:numel(model.rxns))];
plot_mat = plot_mat(plot_idx,:);

% sort distances
[d_sort, sort_idx] = sort(d(plot_idx), 'descend');

% add compartment info to reaction names
names = model.rxnNames(findRxnIDs(model, model.rxns(plot_idx)));
rxns = model.rxns(plot_idx);
for i = 1:numel(names)
    res = regexp(erase(rxns{i}, {'_b','_f'}), '_(?<c>\w)$', 'names');
    cof = regexp(erase(rxns{i}, {'_b','_f'}), '(?<c>NAD[P]*)', 'names');
    
    if ~isempty(cof)
        names{i} = strcat(names{i}, [' ' cof.c]);
    end
    
    if isempty(res)
        met_comps = unique(getCompartment(findMetsFromRxns(model, rxns{i})));
        if numel(met_comps) > 1
            met_comps = setdiff(met_comps, 'c');
        end
        names{i} = strcat(names{i}, [' (' char(met_comps) ')']);
    else
        names{i} = strcat(names{i}, [' (' res.c ')']);
    end
end

% add direction to reaction names
direction = regexp(model.rxns(plot_idx), '_[fb]$', 'match');
direction(cellfun(@isempty, direction)) = {'f'};
direction = [direction{:}]';
direction = erase(direction, '_');
direction = cellfun(@(x)[' (' x ')'], direction, 'un', 0);

% sort y-labels
names = strcat(names, direction);
ylabel_sort = names(sort_idx);


figure
bar(d_sort, 'Horizontal', false)
xticks(1:sum(plot_idx))
xticklabels(ylabel_sort)
xtickangle(90)
ylabel('v^{red} - v^{WT} [mmol/gDW/h]', 'FontSize', 12)
set(gcf, 'OuterPosition', [353.6667   58.3333  704.6667  640.6667])
exportgraphics(gcf, fullfile('..', '..', 'figures', 'changed_rxns_resp_min.png'))

% find changed reactions
fprintf('Reaction ID\tReaction Name\tWT\tminimization\n')
for i = 1:numel(model.rxns)
    if xor(opt_distr(i),fluxes_resp_red(i))
        % reaction has zero flux in one of the conditions
        opt_val = log2(opt_distr(i));
        red_val = log2(fluxes_resp_red(i));
        fprintf('%s\t%s\t%.2f\t%.2f\n', model.rxns{i}, model.rxnNames{i}, opt_distr(i), fluxes_resp_red(i))
    end
end
        
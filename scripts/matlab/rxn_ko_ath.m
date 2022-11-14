% Perform simulation of single and double KOs with the AraCore model

clear;clc

initCobraToolbox(false)

% load models
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

% find indices of reactions that produce or consume CO2
co2_prod_idx = find(any(model.S(findMetIDs(model, co2_ids), :) > 0, 1));
co2_cons_idx = find(any(model.S(findMetIDs(model, co2_ids), :) < 0, 1));

growth_respiration = nan(numel(model.rxns));
net_respiration = nan(numel(model.rxns));
rgr = nan(numel(model.rxns));

n_rxns = numel(model.rxns);
rxns = model.rxns;

% perform single KOs to exclude lethal KOs
param = cplexoptimset;
param.threads = 1;
isLethal = false(size(model.rxns));
solver = getCobraSolver('LP');
environment = getEnvironment;
parfor i=1:numel(model.rxns)
    restoreEnvironment(environment);
    changeCobraSolver(solver, 'LP', 0, -1);
    
    tmp_model = changeRxnBounds(model, model.rxns{i}, 0, 'b');
    fba_sol = optimizeCbModel(tmp_model, 'max', 0, true, param);
    isLethal(i) = fba_sol.f < 1e-9;
end

parfor i=1:n_rxns
    
    if ~isLethal(i)
        tic
        fprintf('#%d: Simulating double KOs for %s\n', rxns{i}, i)
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        row_g = nan(1,n_rxns);
        row_n = nan(1,n_rxns);
        row_rgr = nan(1,n_rxns);
        for j=i:n_rxns
            
            if ~isLethal(j)
                % KO
                if i==j
                    tmp_model = changeRxnBounds(model, rxns{i}, 0, 'b');
                else
                    tmp_model = changeRxnBounds(model, rxns([i j]), 0, 'b');
                end
                
                % solve FBA problem, minimizing the sum of fluxes and without loops
                fba_sol = struct;
                try
                    fba_sol = optimizeCbModel(tmp_model, 'max', 'one', true, param);
                catch ME
                    fba_sol.f = -1;
                    disp(ME)
                end
                % calculate growth and net respiration
                if fba_sol.f > 1e-9
                    model_flux = fba_sol.x;
                    model_flux(model_flux<1e-9) = 0;
                    row_g(j) = sum(model_flux(co2_prod_idx));
                    row_n(j) = sum(model_flux(co2_prod_idx)) - sum(model_flux(co2_cons_idx));
                    row_rgr(j) = fba_sol.f;
                end
                
            end
            growth_respiration(i,:) = row_g;
            net_respiration(i,:) = row_n;
            rgr(i,:) = row_rgr;
        end
        fprintf('Time: %.2f\n', toc/60)
    end
    
end

% save('../../data/ko_results.mat', 'growth_respiration', 'net_respiration', 'rgr')
load('../../data/ko_results.mat')

% get default respiration value
fba_sol = optimizeCbModel(model, 'max', 'one', true, param);
growth_resp_default = sum(fba_sol.x(co2_prod_idx));
rgr_default = fba_sol.f;

% calculate carbon loss
cmf = getCarbonMolarFraction(model, true);
carbon_loss = growth_respiration ./ rgr / cmf;
carbon_loss_default = growth_resp_default / rgr_default / cmf;
carbon_loss_scaled = carbon_loss / carbon_loss_default;

for i = 1:n_rxns-1
    for j = i+1:n_rxns
        carbon_loss_scaled(j,i) = carbon_loss_scaled(i,j);
        growth_respiration(j,i) = growth_respiration(i,j);
        rgr(j,i) = rgr(i,j);
    end
end
    
av_cl_scaled = mean(carbon_loss_scaled, 'omitnan');
sd_cl_scaled = std(carbon_loss_scaled, 'omitnan');
av_gr = mean(growth_respiration, 'omitnan');
sd_gr = std(growth_respiration, 'omitnan');

% only take top 10 percent of reactions and those that are not all-lethal
nz_cols = any(carbon_loss_scaled);
gt_thesh_cols = any(av_cl_scaled>=prctile(av_cl_scaled, 90), 1);

% reduce arrays and matrices accordingly
carbon_loss_scaled_red = carbon_loss_scaled(nz_cols&gt_thesh_cols, nz_cols&gt_thesh_cols);
growth_respiration_red = growth_respiration(nz_cols&gt_thesh_cols, nz_cols&gt_thesh_cols);
sd_cl_scaled_red = sd_cl_scaled(nz_cols&gt_thesh_cols);
av_gr_red = av_gr(nz_cols&gt_thesh_cols);
sd_gr_red = sd_gr(nz_cols&gt_thesh_cols);

% sort values
[av_cl_scaled_red_sorted, sort_idx] = sort(av_cl_scaled(nz_cols&gt_thesh_cols), 'descend');
sd_cl_scaled_red_sorted = sd_cl_scaled_red(sort_idx);
av_gr_red_sorted = av_gr_red(sort_idx);
sd_gr_red_sorted = sd_gr_red(sort_idx);

lab_esc = strrep(rxns(nz_cols&gt_thesh_cols), '_', '\_');

% add compartment info to reaction names
names = model.rxnNames(nz_cols&gt_thesh_cols);
for i = 1:numel(names)
    res = regexp(lab_esc{i}, '_(?<c>\w)$', 'names');
    if isempty(res)
        met_comps = unique(getCompartment(findMetsFromRxns(model, strrep(lab_esc{i}, '\_', '_'))));
        if numel(met_comps) > 1
            met_comps = setdiff(met_comps, 'c');
        end
        names{i} = strcat(names{i}, [' (' char(met_comps) ')']);
    else
        names{i} = strcat(names{i}, [' (' res.c ')']);
    end
end
names{ismember(names, 'Proton sink/source for charge balancing of im-/export (H)')} = 'Proton sink/source (h)';
names{ismember(names, 'GCA/GCEA shuffle (h)')} = 'GCA/GCEA shuttle (h)';

carbon_loss_scaled_red_sort = carbon_loss_scaled_red(sort_idx, sort_idx);
plot_mat = carbon_loss_scaled_red_sort;
for i = 1:size(carbon_loss_scaled_red_sort,1)
    for j = i+1:size(carbon_loss_scaled_red_sort,2)
        plot_mat(j,i) = nan;
    end
end

%% Plot results
figure
tiledlayout(1,5)
nexttile([1,4])
alpha_data = ~isnan(plot_mat) + tril(ones(size(plot_mat)));
imagesc(plot_mat,'AlphaData',alpha_data)
colormap('summer')
cb = colorbar('Location', 'northoutside');
cb.Title.String = '\lambda_{ko} / \lambda_{wt}';
hold on
growth_respiration_red_sort = growth_respiration_red(sort_idx, sort_idx);

for i = 1:size(plot_mat,1)-1
    for j = i+1:size(plot_mat,2)
        if isnan(growth_respiration_red_sort(i,j))
            col = [1 1 1];
        elseif growth_respiration_red_sort(i,j) > growth_resp_default
            col = [.2 .3 .8];
        elseif growth_respiration_red_sort(i,j) == growth_resp_default
            col = 'black';
        else
            col = [.8 .4 .3];
        end
        scatter(i,j,75,'filled', 'MarkerFaceColor', col)
    end
end
xticks(1:sum(nz_cols))
xticklabels(lab_esc(sort_idx))
xtickangle(90)

yticks(1:sum(nz_cols))
yticklabels(names(sort_idx))

set(gcf, 'OuterPosition', [-3.6667   78.3333  946.6667  620.6667]);

h = zeros(3,1);
h(1) = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', [.2 .3 .8],...
    'MarkerEdgeColor', [0 0 0]);
h(2) = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', [.8 .4 .3],...
    'MarkerEdgeColor', [0 0 0]);
h(3) = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', [1 1 1],...
    'MarkerEdgeColor', [0 0 0]);
legend(h, {'R_G increase', 'R_G decrease', 'lethal'}, 'box', 'off', 'color', 'none',...
    'position', [0.83 0.2 0.126967095851216 0.0884375000000001])

nexttile([1,1])
cdata = zeros(numel(av_cl_scaled_red_sorted), 3);
cdata(av_gr_red_sorted>growth_resp_default,:) = repmat([.2 .3 .8], sum(av_gr_red_sorted>growth_resp_default), 1);
cdata(av_gr_red_sorted==growth_resp_default,:) = repmat([0 1 0], sum(av_gr_red_sorted==growth_resp_default), 1);
cdata(av_gr_red_sorted<growth_resp_default,:) = repmat([.8 .4 .3], sum(av_gr_red_sorted<growth_resp_default), 1);
bar(av_cl_scaled_red_sorted(numel(av_cl_scaled_red_sorted):-1:1),...
    'CData', cdata(numel(av_cl_scaled_red_sorted):-1:1,:),...
    'FaceColor', 'flat', 'Horizontal', true);
set(gca, 'Box', 'off')
axis tight
yticks(NaN)
hold on
errorbar(av_cl_scaled_red_sorted,numel(av_cl_scaled_red_sorted):-1:1,...
    [],[],...
    sd_cl_scaled_red_sorted,sd_cl_scaled_red_sorted,...
    'LineStyle', 'none', 'Color', 'k')
xlabel('Average \lambda_{ko} / \lambda_{wt}')

set(gca, 'color', 'none')

exportgraphics(gcf, fullfile('..', '..', 'figures', 'ko_growth_resp_top_10_pct.png'))
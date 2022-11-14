% Determine net CO2 assimilation, carbon use efficiency, maximum CO2
% release/assimilation, and growth respiration

clear;clc

initCobraToolbox(false)

top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'fun_models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.fun_models, tmp.fun_model_publications, tmp.fun_organisms_uniq);

% find CO2 ID in each model
co2_ids = repmat({''}, numel(models), 1);
for i=1:numel(models)
    model = models{i};
    if isfield(model, 'metFormulas')
        co2_form_match = ismember(lower(model.metFormulas), 'co2');
        if any(co2_form_match)
            co2_ids{i} = model.mets(co2_form_match);
        end
        
    end
    
    if ~isfield(model, 'metFormulas') || isempty(co2_form_match) || ~any(co2_form_match)
        co2_id_match = startsWith(lower(model.mets), {'co2', 'carbon_dioxide', 'carbon dioxide', 'carbon-dioxide'});
        if any(co2_id_match)
            co2_ids{i} = model.mets(co2_id_match);
        else
            co2_id_match = startsWith(lower(model.metNames), {'co2', 'carbon_dioxide', 'carbon dioxide', 'carbon-dioxide'});
            if any(co2_id_match)
                co2_ids{i} = model.mets(co2_id_match);
            else
                
            end
        end
    end
end

co2_ids{strcmp(model_publications, 'Prigent2014')} = 'META19193[cell]';
co2_ids{strcmp(model_publications, 'Chatterjee2017')} = {'met_237[DefaultCompartment]',...
    'CO2_str[DefaultCompartment]'};
co2_ids{strcmp(model_publications, 'Moreira2019')} = {'CARBON-DIOXIDE_c[c]',...
    'CARBON-DIOXIDE_m[c]', 'CARBON-DIOXIDE_p[c]', 'CARBON-DIOXIDE_x[c]'};
co2_ids{strcmp(model_publications, 'ShawCheung2019')} = {'CARBON-DIOXIDE_c[c]',...
    'CARBON-DIOXIDE_m[c]', 'CARBON-DIOXIDE_p[c]', 'CARBON-DIOXIDE_x[c]'};

%% find carboxylation and oxygenation reactions
v_c = {
    'RBC_h',...                                     ArnoldNikoloski2014
    'R00024',...                                    Botero2018
    'chl_Rubisco',...                               Chatterjee2017
    'R00024__plst'...                               Cunha2022
    'R00024_p',...                                  Simons2014
    'RBPCh',...                                     Gerlin2022
    'RBPCh',...                                     Imam2015
    'RUBISC_h',...                                  Levering2016
    'Ar0065',...                                    Klanchui2018
    'RBPCh'...                                      Loira2017
    'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p',...   Moreira2019
    'META45956',...                                 Prigent2014
    'NanoG0589',...                                 Shah2017
    'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p',...   ShawCheung2019
    'RBPCh',...                                     Zuniga2016
    };

v_o = {
    'RBO_h',...                                     ArnoldNikoloski2014
    'R03140',...                                    Botero2018
    'chl_RuBPOxid',...                              Chatterjee2017
    'R03140__plst'...                               Cunha2022
    'R03140_p',...                                  Simons2014
    'RBCh_1',...                                    Gerlin2022
    'RBCh',...                                      Imam2015
    'RUBISO_h',...                                  Levering2016
    'Ar0066',...                                    Klanchui2018
    'RBCh',...                                      Loira2017
    'RXN-961_p',...                                 Moreira2019
    'META51505',...                                 Prigent2014
    'NanoG0460',...                                 Shah2017
    'RXN-961_p',...                                 ShawCheung2019
    'RBCh',...                                      Zuniga2016
    };

% Klanchui2018: Ar0064 is carboxylation and oxygenation combined with a
% ratio

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

for i=1:numel(models)
    if ~contains(organisms_uniq{i}, {'Chlamydomonas', 'Chlorella'})
        model = models{i};
        % fix ratio between RUBISCO carbocylation and oxygenation within a
        % tolerated range
        model = addMetabolite(model, 'phi_lb',...
            'csense', 'L');
        model.S(findMetIDs(model, 'phi_lb'), findRxnIDs(model, [v_c(i) v_o(i)])) = [phi-phi_tol -1];
        
        model = addMetabolite(model, 'phi_ub',...
            'csense', 'L');
        model.S(findMetIDs(model, 'phi_ub'), findRxnIDs(model, [v_c(i) v_o(i)])) = [-phi-phi_tol 1];
        
        
        % check whether growth is still possible
        s = optimizeCbModel(model);
        if ~(s.f > 0)
            fprintf('Organism: %s\n', organisms_uniq{i})
            fprintf('Imposed oxygenation to carboxylation ration leads to infeasibility\n')
        else
            models{i} = model;
        end
    end
end

%% Solve FBA and extract flux sum through CO2 producing (and consuming) reactions
% COBRA solver parameters
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'timeLimit', 180);

sum_co2_prod = nan(numel(models), 1);
sum_co2_cons = nan(numel(models), 1);
max_co2_prod = nan(numel(models), 1);
max_co2_cons = nan(numel(models), 1);
obj_norm_one_loopless = nan(numel(models), 1);
model_flux = cell(numel(models), 1);

for i=1:numel(models)
    fprintf('Processing model #%d: %s\n', i, organisms_uniq{i})
    model = convertToIrreversible(models{i});
    
    % change default flux bounds if they are larger than 1000 mmol/gDW/h
    model.ub(model.ub > 1000) = 1000;
    
    % add constraint sense field if not present
    if ~isfield(model, 'csense')
        model.csense = repmat('E', numel(model.mets), 1);
    end
    
    % find indices of reactions that produce or consume CO2
    co2_prod_idx = any(model.S(findMetIDs(model, co2_ids{i}), :) > 0, 1);
    co2_cons_idx = any(model.S(findMetIDs(model, co2_ids{i}), :) < 0, 1);
    
    % solve FBA problem, minimizing the sum of fluxes and without loops
    fba_sol = optimizeCbModel(model, 'max', 'one', false);
    bio_opt = fba_sol.f;
    if bio_opt > 0
        % save sum of CO2-producing and consuming fluxes      
        sum_co2_prod(i) = sum(fba_sol.x(co2_prod_idx));
        sum_co2_cons(i) = sum(fba_sol.x(co2_cons_idx));
        model_flux{i} = fba_sol.x;
        obj_norm_one_loopless(i) = bio_opt;
        
        % maximize CO2 production
        tmp_model = model;
        co2_rxn_idx = any(tmp_model.S(findMetIDs(tmp_model, co2_ids{i}),:)~=0, 1);
        co2_exc_idx = sum(tmp_model.S~=0, 1) == 1 & co2_rxn_idx;
        exc_co2_met = findMetsFromRxns(tmp_model, tmp_model.rxns(co2_exc_idx));
        tmp_model = addReaction(tmp_model, 'co2_release',...
            'metaboliteList', exc_co2_met,...
            'stoichCoeffList', -ones(numel(exc_co2_met), 1),...
            'reversible', false);
        tmp_model = addMetabolite(tmp_model, 'fix_objective');
        tmp_model.S(end, tmp_model.c==1) = 1;
        tmp_model.b(end) = 0.99 * bio_opt;
        tmp_model.csense(end) = 'G';
        tmp_model = changeObjective(tmp_model, 'co2_release');
        fba_sol = optimizeCbModel(tmp_model, 'max', 'one', false);
        max_co2_prod(i) = fba_sol.f;
        
        % maximize CO2 fixation
        tmp_model = model;
        co2_rxn_idx = any(tmp_model.S(findMetIDs(tmp_model, co2_ids{i}),:)~=0, 1);
        co2_exc_idx = sum(tmp_model.S~=0, 1) == 1 & co2_rxn_idx;
        exc_co2_met = findMetsFromRxns(tmp_model, tmp_model.rxns(co2_exc_idx));
        tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns(co2_exc_idx),...
            0, 'b');
        tmp_model = addMetabolite(tmp_model, 'fix_objective');
        tmp_model.S(end, tmp_model.c==1) = 1;
        tmp_model.b(end) = 0.99 * bio_opt;
        tmp_model.csense(end) = 'G';
        tmp_model = addReaction(tmp_model, 'co2_uptake',...
            'metaboliteList', exc_co2_met,...
            'stoichCoeffList', ones(numel(exc_co2_met), 1),...
            'reversible', false);
        tmp_model = addReaction(tmp_model, 'co2_release',...
            'metaboliteList', exc_co2_met,...
            'stoichCoeffList', -ones(numel(exc_co2_met), 1),...
            'reversible', false);
        tmp_model = changeObjective(tmp_model, {'co2_uptake', 'co2_release'}, [1 -1]);
        
        fba_sol = optimizeCbModel(tmp_model, 'max', 'one', false);
        max_co2_cons(i) = fba_sol.f;
    end
end

save(fullfile('..', '..', 'data', 'co2_rxn_flux'), 'model_flux', 'sum_co2_prod', 'sum_co2_cons',...
    'obj_norm_one_loopless', 'max_co2_prod', 'max_co2_cons')

v_c_flux = nan(numel(models), 1);
v_o_flux = nan(numel(models), 1);
for i=1:numel(models)
    if ~isempty(model_flux{i})
        v_c_flux(i) = model_flux{i}(findRxnIDs(models{i}, v_c{i}));
        v_o_flux(i) = model_flux{i}(findRxnIDs(models{i}, v_o{i}));
    end
end

nz_idx = ~cellfun(@isempty, model_flux);
net_fixation = sum_co2_cons - sum_co2_prod;
net_fixation(net_fixation<1e-9) = 0;

growth_respiration = sum_co2_prod;
dark_resp = v_c_flux(nz_idx) - 0.5*v_o_flux(nz_idx) - net_fixation(nz_idx);
cue = 1 - (0.5*v_o_flux + dark_resp) ./ v_c_flux(nz_idx);

%% plot results
x_lab = strrep(organisms_uniq(nz_idx), '_', ' ');
x_lab_split = cellfun(@strsplit, x_lab, 'un', 0);
for i=1:numel(x_lab)
    x_lab{i} = ['{\it' strjoin(x_lab_split{i}([1 2]), ' ') '}'];
    if numel(x_lab_split{i}) > 2
        x_lab{i} = [x_lab{i} ' ' strjoin(x_lab_split{i}(3:end), ' ')];
    end
end

letter_coords = [-0.12, 0.965];
letter_fz = 10;
axis_fz = 8;
ylab_x_pos = -6.5;
bg = [1 1 1];
lw = 1.3;
tl = [0.01 0.01];
plot_pos = [0.3 0.6 0.1];

close all
figure

% CO2 fixation in FBA solution
subplot(5,1,1)
bar(net_fixation(nz_idx), 'FaceColor', [.3 .4 .6])
y_lim = get(gca, 'YLim');
ylabel({'A_{net} (pFBA)', '[mmol/gDW/h]'})
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'HorizontalAlignment','left',...
    'VerticalAlignment', 'middle', 'Position', [ylab_x_pos 0.5*y_lim(2)],...
    'Units', 'normalized')
xticks(NaN)
set(gca, 'FontSize', axis_fz, 'FontName', 'Arial', 'YScale', 'log',...
    'Box', 'off', 'Color', bg, 'linewidth', lw, 'TickLength', tl)
yticks = logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 4);
set(gca, 'YTick', logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 4),...
    'YTickLabel', strtrim(cellstr(num2str(round(log10(yticks(:))), '10^{%d}'))))
text(letter_coords(1), letter_coords(2), 'a', 'units', 'normalized',...
    'FontSize', letter_fz, 'FontWeight', 'bold', 'FontName', 'Arial')
axh = gca();
axh.Position([1 3 4]) = plot_pos;
axh.Position(2) = axh.Position(2) + 0.08;

% carbon use efficiency
subplot(5,1,2)
bar(cue, 'FaceColor', [.3 .4 .6])
y_lim = get(gca, 'YLim');
y_lim(2) = 1.3*max(cue(~isinf(cue)));
ylim(y_lim)
ylabel('CUE (pFBA)')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'HorizontalAlignment','left',...
    'VerticalAlignment', 'middle', 'Position', [ylab_x_pos 0.5*y_lim(2)],...
    'Units', 'normalized')
xticks(NaN)
set(gca, 'FontSize', axis_fz, 'FontName', 'Arial', 'YScale', 'log',...
    'Box', 'off', 'Color', bg, 'linewidth', lw, 'TickLength', tl)
yticks = logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 5);
set(gca, 'YTick', logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 5),...
    'YTickLabel', strtrim(cellstr(num2str(round(log10(yticks(:))), '10^{%d}'))))
text(letter_coords(1), letter_coords(2), 'b', 'units', 'normalized',...
    'FontSize', letter_fz, 'FontWeight', 'bold', 'FontName', 'Arial')
axh = gca();
axh.Position([1 3 4]) = plot_pos;
axh.Position(2) = axh.Position(2) + 0.08;

% maximum CO2 release
subplot(5,1,3)
bar(max_co2_prod(nz_idx), 'FaceColor', [.3 .4 .6])
y_lim = get(gca, 'YLim');
y_lim(2) = 1.3*max(max_co2_prod(nz_idx));
ylim(y_lim)

ylabel({'R_{net}^{max}', '[mmol/gDW/h]'})
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'HorizontalAlignment','left',...
    'VerticalAlignment', 'middle', 'Position', [ylab_x_pos 0.5*y_lim(2)],...
    'Units', 'normalized')
xticks(NaN)

set(gca, 'FontSize', axis_fz, 'FontName', 'Arial', 'YScale', 'log',...
    'Box', 'off', 'Color', bg, 'linewidth', lw, 'TickLength', tl)
yticks = logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 7);
set(gca, 'YTick', logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 7),...
    'YTickLabel', strtrim(cellstr(num2str(round(log10(yticks(:))), '10^{%d}'))))
text(letter_coords(1), letter_coords(2), 'c', 'units', 'normalized',...
    'FontSize', letter_fz, 'FontWeight', 'bold', 'FontName', 'Arial')
axh = gca();
axh.Position([1 3 4]) = plot_pos;
axh.Position(2) = axh.Position(2) + 0.08;

% maximum CO2 fixation
subplot(5,1,4)
bar(max_co2_cons(nz_idx), 'FaceColor', [.3 .4 .6])
y_lim = get(gca, 'YLim');
ylabel({'A_{net}^{max}', '[mmol/gDW/h]'})
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'HorizontalAlignment','left',...
    'VerticalAlignment', 'middle', 'Position', [ylab_x_pos 0.5*y_lim(2)],...
    'Units', 'normalized')
xticks(NaN)
y_lim(2) = 1.3*max(max_co2_cons(nz_idx));
ylim(y_lim)
set(gca, 'FontSize', axis_fz, 'FontName', 'Arial', 'YScale', 'log',...
    'Box', 'off', 'Color', bg, 'linewidth', lw, 'TickLength', tl)
yticks = logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 8);
set(gca, 'YTick', logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 8),...
    'YTickLabel', strtrim(cellstr(num2str(round(log10(yticks(:))), '10^{%d}'))))
text(letter_coords(1), letter_coords(2), 'd', 'units', 'normalized',...
    'FontSize', letter_fz, 'FontWeight', 'bold', 'FontName', 'Arial')
axh = gca();
axh.Position([1 3 4]) = plot_pos;
axh.Position(2) = axh.Position(2) + 0.08;

% growth respiration
subplot(5,1,5)
bar(growth_respiration(nz_idx), 'FaceColor', [.3 .4 .6])
y_lim = get(gca, 'YLim');
ylabel({'R_g', '[mmol/gDW/h]'})
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'HorizontalAlignment','left',...
    'VerticalAlignment', 'middle', 'Position', [ylab_x_pos 0.5*y_lim(2)],...
    'Units', 'normalized')
xticks(1:sum(nz_idx))
xticklabels(x_lab)
xtickangle(45)
y_lim(2) = 1.3*max(max_co2_cons(nz_idx));
ylim(y_lim)
set(gca, 'FontSize', axis_fz, 'FontName', 'Arial', 'YScale', 'log',...
    'Box', 'off', 'Color', bg, 'linewidth', lw, 'TickLength', tl)
yticks = logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 8);
set(gca, 'YTick', logspace(log10(y_lim(1)+1e-6),log10(y_lim(2)), 8),...
    'YTickLabel', strtrim(cellstr(num2str(round(log10(yticks(:))), '10^{%d}'))))
text(letter_coords(1), letter_coords(2), 'e', 'units', 'normalized',...
    'FontSize', letter_fz, 'FontWeight', 'bold', 'FontName', 'Arial')
axh = gca();
axh.Position([1 3 4]) = plot_pos;
axh.Position(2) = axh.Position(2) + 0.08;

set(gcf, 'OuterPosition', [589.0000   45.6667  464.0000  653.3333])
exportgraphics(gcf, fullfile('..', '..', 'figures', 'resp_anet_cue_models.png'))
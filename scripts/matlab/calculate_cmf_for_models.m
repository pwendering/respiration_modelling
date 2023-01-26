% calculate carbon molar fractions in the biomass reaction for the models 
% that are used in the analyses for Figure 5 and 6

clear;clc

%% Load models
top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'fun_models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.fun_models, tmp.fun_model_publications, tmp.fun_organisms_uniq);

%% calculate carbon molar fraction in biomass
cmf = nan(numel(models), 1);
for i = 1:numel(models)
    if isfield(models{i}, 'metFormulas')
        if isequal(model_publications{i}, 'ArnoldNikoloski2014')
            cmf(i) = getCarbonMolarFraction(models{i}, true, true);
        else
            try
                cmf(i) = getCarbonMolarFraction(models{i}, false, true);
            catch
                cmf(i) = NaN;
            end
        end
    end
end

%% Plot results as barplot
x_lab = strrep(organisms_uniq, '_', ' ');
x_lab_split = cellfun(@strsplit, x_lab, 'un', 0);
for i=1:numel(x_lab)
    x_lab{i} = ['{\it' strjoin(x_lab_split{i}([1 2]), ' ') '}'];
    if numel(x_lab_split{i}) > 2
        x_lab{i} = [x_lab{i} ' ' strjoin(x_lab_split{i}(3:end), ' ')];
    end
end

valid_cmf_idx = ~isnan(cmf);
bar(cmf(valid_cmf_idx))
xticks(1:sum(valid_cmf_idx))
xticklabels(x_lab(valid_cmf_idx))
xtickangle(90)
ylabel({'Carbon molar fraction', '[mmol/gDW]'})
set(gca,...
    'box', 'on',...
    'linewidth', 1.3,...
    'fontsize', 14,...
    'fontname', 'Arial')
set(gcf, 'OuterPosition', [337.6667   95.6667  476.6667  625.3333])
exportgraphics(gcf, fullfile('..', '..', 'figures', 'carbon_model_fraction_si.png'),...
    'Resolution', 400);



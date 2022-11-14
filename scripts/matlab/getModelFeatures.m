% Determined general characteristics of the models that are compared for
% Figure 5 and 6
% Number of 
% * genes
% * reactions
% * metabolites
% * compartments

top_model_dir = fullfile('..', '..', 'models');

tmp = load(fullfile(top_model_dir, 'models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.models, tmp.model_publications, tmp.organisms_uniq);

non_empty_idx = find(~cellfun(@isempty,models));

model_publications = model_publications(non_empty_idx);
organisms_uniq = organisms_uniq(non_empty_idx);

n_rxns = cellfun(@(x)numel(x.rxns),models(non_empty_idx));
n_mets = cellfun(@(x)numel(x.mets),models(non_empty_idx));
n_genes = cellfun(@(x)numel(x.genes),models(non_empty_idx));
n_genes(n_genes==0) = NaN;
n_comps = nan(size(models(non_empty_idx)));
for i = 1:numel(models(non_empty_idx))
    if isfield(models{non_empty_idx(i)}, 'comps')
        n_comps(i) = numel(models{non_empty_idx(i)}.comps);
    end
end
n_comps(isnan(n_comps)) = 1;
% info from paper (other changes were introduced later to the exported
% table)
n_comps(strcmp(model_publications, 'Chatterjee2017')) = 4;
n_comps(strcmp(model_publications, 'ShawCheung2019')) = 4;

tab = cell2table([model_publications strrep(organisms_uniq, '_', ' '),...
    num2cell([n_genes n_rxns n_mets n_comps])],...
    'VariableNames', {'AuthorYear', 'Species', 'Genes', 'Reactions',...
    'Metabolites', 'Compartments'});
    
writetable(tab, fullfile('..','..','data','model_features.xlsx'))

% Test if FBA can yield non-zero growth for all models
top_model_dir = fullfile('..', '..', 'models');

tmp = load(fullfile(top_model_dir, 'models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.models, tmp.model_publications, tmp.organisms_uniq);

non_empty_idx = find(~cellfun(@isempty,models));

for i=non_empty_idx'
    fprintf('Current model: %s (%s)\n', organisms_uniq{i}, model_publications{i})
    bio_rxn_idx = contains(models{i}.rxnNames, 'biomass', 'IgnoreCase', true) | ...
        contains(models{i}.rxns, 'biomass', 'IgnoreCase', true);
    
    if sum(bio_rxn_idx) == 1
        fprintf('--> biomass reaction found:\t%s\n', models{i}.rxns{bio_rxn_idx})
    elseif sum(bio_rxn_idx) > 1
        fprintf('--> multiple reaction matches found for keyword ''biomass'' at positions ')
        fprintf('%s\n', regexprep(num2str(find(bio_rxn_idx')),'\ +',','))
    else
        bio_met_idx = contains(models{i}.metNames, 'biomass', 'IgnoreCase', true) | ...
            contains(models{i}.mets, 'biomass', 'IgnoreCase', true);
        
        if sum(bio_met_idx) == 1
            bio_rxn = findRxnsFromMets(models{i}, models{i}.mets(bio_met_idx));
            fprintf('--> biomass reaction found via metabolite:\t%s\n', char(bio_rxn))
        elseif sum(bio_met_idx) > 1
            fprintf('--> multiple metabolite matches found for keyword ''biomass'' at positions ')
            fprintf('%s\n', regexprep(num2str(find(bio_met_idx')),'\ +',','))
        else
            fprintf('X no matches found for keyword ''biomass'' in both reactions and metabolites.\n')
        end
        
    end
end

bio_rxn_ids = {...
    {'Bio_opt'}...                    % A. thaliana (ArnoldNikoloski2014)                        
    {'RBS01'}...                      % S. tuberosum (Botero2018)
    {''}...                           % O. sativa (Chatterjee2017) [38 draw reactions]
    {'e_Biomass_Leaf__cyto'}...       % Q. suber (Cunha2022)
    {''}...                           % S. officinarum (DalMolin2010) [no biomass reaction found]
    {''}...                           % S. bicolor (DalMolin2010) [no biomass reaction found]
    {'Bio_Nplus'} ...                 % Z. mays (Simons2011)
    {''}...                           % S. italica (DalMolin2016) [no model available]
    {'BIOMASS_LEAF'}...               % S. lycopersicum (Gerlin2022)
    {'Biomasssynth_u'}...             % B. napus (Hay2014)
    {'Biomass_Chlamy_auto'}...        % C. reinhardtii (Imam2015)
    {'bof_c'}...                      % P. tricornutum (Levering2016)
    {'Biomass_Nanno_auto'}...         % N. salina (Loira2017)
    models{15}.rxns(contains(models{15}.rxns, 'biomass', 'IgnoreCase', true))...
    ...                               % G. max (Moreira2019) [109 draw reactions]
    {'BiomassShoot'}...               % M. truncatula (Pfau2018)
    {''}...                           % E. siliculosus (Prigent2014) [no biomass reaction found]
    {'BiomassRxn'}...                 % P. trichocarpa (SarkarMaranas2020)
    {'Biomass_NanoG_auto'}...         % N. gaditana (Shah2017)
    models{20}.rxns(contains(models{20}.rxns, 'biomass', 'IgnoreCase', true))...                           
    ...                               % S. viridis (ShawCheung2019) [40 draw reactions]
    {'Biomass_Cvu_auto-'}...          % C. vulgaris (Zuniga2016)
    };

%     {'R998'}...                       % H. vulgare (GrafahrendBelau2009)
%     {'Biomass'}...                    % Mentha x piperita (Johnson2017)
%     {'Ar0828'}...                     % A. plantensis C1 (Klanchui2018)

fba_solutions = nan(size(models));
for i=1:numel(models)
    if ~isempty(models{i})
         model = models{i};
         
         % set objective
         bio_idx = findRxnIDs(model,bio_rxn_ids{i});
         if bio_idx ~= 0
            model.c = zeros(size(model.S,2),1);
            model.c(findRxnIDs(model,bio_rxn_ids{i})) = 1;
         end
         
         if strcmp(model_publications{i},'Chatterjee2017')
             tmp = load(fullfile(top_model_dir, 'Chatterjee2017_Oryza_sativa',...
                 'data_files', 'biomass_rxns.mat'));
             bio_rxns = tmp.bio_rxns;
             model.c(findRxnIDs(model, bio_rxns)) = 1;
             model.lb(findRxnIDs(model, bio_rxns)) = 0;
             model.S(:,findRxnIDs(model, bio_rxns)) = -model.S(:,findRxnIDs(model, bio_rxns));
         elseif strcmp(model_publications{i},'Shameer2018')
             model = changeRxnBounds(model, {'Sucrose_tx1', 'Sucrose_tx2',...
                 'GLC_tx1', 'GLC_tx2'}, 0, 'b');
         elseif strcmp(model_publications{i},'Moreira2019')
             model.lb(findRxnIDs(model,bio_rxn_ids{i})) = 0;
             model.ub(findRxnIDs(model,bio_rxn_ids{i})) = 1000;
             model.S(:,findRxnIDs(model,bio_rxn_ids{i})) = -model.S(:,findRxnIDs(model,bio_rxn_ids{i}));
         elseif strcmp(model_publications{i},'ShawCheung2019')
             model.lb(findRxnIDs(model,bio_rxn_ids{i})) = 0;
             model.S(:,findRxnIDs(model,bio_rxn_ids{i})) = -model.S(:,findRxnIDs(model,bio_rxn_ids{i}));
         elseif strcmp(model_publications{i},'Hay2014')
             model = addExchangeRxn(model,{'no3[e]', 'ph[e]', 'co2[e]', 'o2[e]', 'nh4[e]', 'so4[e]', 'h2o[c]', 'pi[e]','h2s[p]', 'h[a]'});
         elseif strcmp(model_publications{i},'Prigent2014')
             % biomass reaction
             tmp = load(fullfile(top_model_dir, 'Prigent2014_Ectocarpus_siliculosus',...
                 'data_files', 'bio_rxn.mat'));
             [bio_mets, bio_coeffs] = deal(tmp.bio_mets, tmp.bio_coeffs);
             model = addReaction(model,...
                 'biomass',...
                 'metaboliteList', [bio_mets; {'biomass[cell]'}],...
                 'stoichCoeffList', [-bio_coeffs; 1],...
                 'objectiveCoef', 1,...
                 'reversible', false);
             model = addExchangeRxn(model, 'biomass[cell]', 0, 1000);
             
             % exchange reactions
             tmp = load(fullfile(top_model_dir, 'Prigent2014_Ectocarpus_siliculosus',...
                 'data_files', 'exc_mets.mat'));
             exc_mets = tmp.exc_mets;
             model = addExchangeRxn(model, exc_mets);
         elseif strcmp(model_publications{i},'Pfau2018')
             model = changeRxnBounds(model,...
                 {'R_TCE_CARBON-DIOXIDE',...
                 'R_TEC_CARBON-DIOXIDE',...
                 'R_TEC_CO+2',...
                 'R_TEC_CPD-3',...
                 'R_TEC_FE+2',...
                 'R_TEC_MG+2',...
                 'R_TEC_NITRATE',...
                 'R_TEC_OXYGEN-MOLECULE',...
                 'R_TEC_Pi',...
                 'R_TEC_SULFATE',...
                 'R_TCE_WATER',...
                 'R_TCE_PROTON',...
                 'R_TCE_HCO3',...
                 'R_TEH_Light',...
                 'R_TEC_WATER',...
                 'R_TCE_THIAMINE-PYROPHOSPHATE',...
                 'R_TCE_THF',...
                 'R_TCE_PYRIDOXAL_PHOSPHATE',...
                 'R_TCE_MYO-INOSITOL',...
                 'R_TCE_FAD',...
                 'R_TCE_CO-A',...
                 'R_TCE_CHOLINE',...
                 'R_TCE_ADENOSYL-HOMO-CYS',...
                 'R_TCE_L-ALPHA-ALANINE',...
                 'R_TCE_L-ASPARTATE'...
                 }, -1000, 'l');
         end
         
         % run FBA
         sol = optimizeCbModel(model);
         fba_solutions(i) = sol.f;
         
         models{i} = model;
    end
end

keep_idx = fba_solutions > 0;
fun_models = models(keep_idx);
fun_model_publications = model_publications(keep_idx);
fun_organisms_uniq = organisms_uniq(keep_idx);
save(fullfile(top_model_dir, 'fun_models'),'fun_models','fun_model_publications', 'fun_organisms_uniq')

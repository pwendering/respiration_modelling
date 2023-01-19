% compate sets of EC numbers of CO2-producing reactions with the Zea mays
% model

top_model_dir = fullfile('..', '..', 'models');
tmp = load(fullfile(top_model_dir, 'models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.models, tmp.model_publications, tmp.organisms_uniq);
non_empty_idx = find(~cellfun(@isempty,models));

% find CO2 ID in each model
co2_ids = repmat({''}, numel(models), 1);
for i=non_empty_idx'
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

% find CO2 producing reactions
co2_prod_rxns = repmat({''}, numel(models), 1);
dummy_models = models;
for i=non_empty_idx'
    model = models{i};
    co2_prod_idx = any(model.S(findMetIDs(model, co2_ids{i}), :) > 0, 1);
    co2_prod_rxns{i} = model.rxns(co2_prod_idx);
    model.rxns = model.rxns(co2_prod_idx);
    model.rxnNames = model.rxnNames(co2_prod_idx);
    if isfield(model, 'rxnECNumbers')
        model.rxnECNumbers = model.rxnECNumbers(co2_prod_idx);
    end
    if isfield(model, 'eccodes')
        model.eccodes = model.eccodes(co2_prod_idx);
    end
    if isfield(model, 'rxnNotes')
        model.rxnNotes = model.rxnNotes(co2_prod_idx);
    end
    
    dummy_models{i} = model;
end

ec_nums = cell(numel(non_empty_idx), 1);
for i=non_empty_idx'
    model = dummy_models{i};
    fprintf('Current model: %s (%s)\n', organisms_uniq{i}, model_publications{i})
    if isfield(model, 'rxnECNumbers')
        % COBRA style
        ec_f_name = 'rxnECNumbers';
        fprintf('--> EC field: %s\n', ec_f_name)
    elseif isfield(model, 'eccodes')
        % RAVEN style
        ec_f_name = 'eccodes';
        fprintf('--> EC field: %s\n', ec_f_name)
    else
        % try to find EC field
        f_names = fieldnames(model);
        potential_ec_fields = f_names(contains(f_names, 'ec', 'i', 1));
        if ~isempty(potential_ec_fields)
            fprintf('--> Potential EC number fields found: %s\n', strjoin(potential_ec_fields, ', '))
        end
        ec_f_name = '';
    end
    
    % if not EC field, try to find EC numbers in reaction names, IDs, and
    % notes
    if ~isempty(ec_f_name)
        tmp_ec = regexp(model.(ec_f_name), '\d+\.\d+\.\d+\.\d+', 'match');
        tmp_ec = [tmp_ec{:}]';
        ec_nums{i} = tmp_ec;
        clear tmp_ec
    else
        fprintf('Attempt to find EC numbers in reactions\n')
        
        tmp_ec = regexp(model.rxns, '\d+\.\d+\.\d+\.\d+', 'match');
        tmp_ec_from_rxns = [tmp_ec{:}]';
        fprintf('--> %d EC numbers found in ''rxnNames'' field.\n', numel(tmp_ec_from_rxns))
        
        tmp_ec = regexp(model.rxnNames, '\d+\.\d+\.\d+\.\d+', 'match');
        tmp_ec_from_rxnnames = [tmp_ec{:}]';
        fprintf('--> %d EC numbers found in ''rxnNames'' field.\n', numel(tmp_ec_from_rxnnames))
        
        if isfield(model, 'rxnNotes')
            tmp_ec = regexp(model.rxnNotes, '\d+\.\d+\.\d+\.\d+', 'match');
            tmp_ec_from_rxnnotes = [tmp_ec{:}]';
            fprintf('--> %d EC numbers found in ''rxnNotes'' field.\n', numel(tmp_ec_from_rxnnotes))
        else
            tmp_ec_from_rxnnotes = '';
        end
        
        max_n_ec = max([numel(tmp_ec_from_rxns) numel(tmp_ec_from_rxnnames) numel(tmp_ec_from_rxnnotes)]);
        
        if numel(tmp_ec_from_rxns) == max_n_ec && ~isempty(tmp_ec_from_rxns)
            ec_nums{i} = tmp_ec_from_rxns;
        elseif numel(tmp_ec_from_rxnnames) == max_n_ec && ~isempty(tmp_ec_from_rxnnames)
            ec_nums{i} = tmp_ec_from_rxnnames;
        elseif ~isempty(tmp_ec_from_rxnnotes)
            ec_nums{i} = tmp_ec_from_rxnnotes;
        end
        
    end
end
clear max_n_ec tmp_ec tmp_ec_from_rxns tmp_ec_from_rxnnames tmp_ec_from_rxnnotes ...
    f_names potential_ec_fields ec_f_name

%% translate reaction IDs to MNXref namespace
mnx_xref = readtable(fullfile('..', '..', 'id_translation', 'reac_xref.tsv'),...
    'FileType', 'text',...
    'CommentStyle', '#',...
    'Delimiter', '\t',...
    'ReadVariableNames', false...
    );
mnx_prop = readtable(fullfile('..', '..', 'id_translation', 'reac_prop.tsv'),...
    'FileType', 'text',...
    'CommentStyle', '#',...
    'Delimiter', '\t',...
    'ReadVariableNames', false...
    );

model_dirs = {
    'Botero2018_Solanum_tuberosum'
    'DalMolin2010_Saccharum_officinarum'
    'DalMolin2010_Sorghum_bicolor'
    'Imam2015_Chlamydomonas_reinhardtii'
    'Pfau2018_Medicago_truncatula'
    'Prigent2014_Ectocarpus_siliculosus'
    'Zuniga2016_Chlorella_vulgaris_UTEX_395'
    'Loira2017_Nannochloropsis_salina'
    };

sources = {'kegg', 'kegg', 'kegg', 'bigg', 'metacyc', 'metacyc', 'bigg', 'bigg'};

all_dirs = strcat(model_publications, strcat('_', organisms_uniq));
for i=1:numel(model_dirs)
    model_idx = ismember(all_dirs, model_dirs{i});
    content = ls(fullfile(top_model_dir, model_dirs{i}));
    if ~any(contains(cellstr(content), 'ec_resp.txt'))
        
        model = dummy_models{model_idx};
        if strcmp(model_publications{model_idx}, 'Prigent2014')
            model.rxns = erase(model.rxnNames, 'metacyc:');
        else
            model.rxns = strtok(model.rxns, {'[', '_'});
        end
        writeEC(model, sources{i}, mnx_xref, mnx_prop, fullfile(top_model_dir, model_dirs{i}))
    end
    
    fid = fopen(fullfile(top_model_dir, model_dirs{i}, 'ec_resp.txt'));
    ec = textscan(fid, '%s');
    fclose(fid);
    
    ec_nums{model_idx} = [ec{:}];
end
clear mnx_prop mnx_xref all_dirs ec fid model_dirs model_idx content sources
save('ec_models_resp', 'ec_nums')

%% remove incomplete EC numbers from comparison
for i=1:numel(ec_nums)
    if ~isempty(ec_nums{i})
        tmp_ec = cellfun(@(x)regexp(x, '\d+\.\d+\.\d+\.\d+', 'match'), ec_nums{i}, 'un', 0);
        ec_nums{i} = [tmp_ec{:}]';
    end
end

%% Compare EC number sets with the set of CO2-producing reactions
n_per_model = cellfun(@(x)numel(unique(x)), ec_nums);
ec_ref = ec_nums{n_per_model == max(n_per_model)};

ec_ovlp = nan(size(non_empty_idx));
ec_ji = nan(size(non_empty_idx));
ec_ovlp_pairwise = nan(numel(non_empty_idx));
ec_ji_pairwise = nan(numel(non_empty_idx));

for i=1:numel(non_empty_idx)
    if ~isempty(ec_nums{non_empty_idx(i)})
        tmp_ec1 = ec_nums{non_empty_idx(i)};
        
        n_is = numel(intersect(tmp_ec1, ec_ref));
        n_un = numel(union(tmp_ec1, ec_ref));
        
        ec_ovlp(i) = n_is;
        ec_ji(i) = n_is/n_un;
        
        for j=1:numel(non_empty_idx)
            
            if ~isempty(ec_nums{non_empty_idx(j)})
                tmp_ec2 = ec_nums{non_empty_idx(j)};
                
                n_is = numel(intersect(tmp_ec1, tmp_ec2));
                n_un = numel(union(tmp_ec1, tmp_ec2));
                
                ec_ovlp_pairwise(i,j) = n_is;
                ec_ji_pairwise(i,j) = n_is/n_un;
            end
        end
    end
end

%% plot results
%{
subplot(1,2,2)
nan_idx = isnan(ec_ji);

bar(ec_ji(~nan_idx),...
    'FaceColor', [.3 .4 .6],...
    'FaceAlpha', .8,...
    'EdgeColor', [.6 .6 .6],...
    'Horizontal', true)

xlabel('Jaccard Index')
yticks(NaN)
xticks(0:0.2:1)
set(gca,...
    'FontSize', 14,...
    'Box', 'off',...
    'TickLength',[0.01, 0.01])
%}
%% heatmap
labels = strrep(organisms_uniq(non_empty_idx), '_', ' ');
labels_split = cellfun(@strsplit, labels, 'un', 0);
for i=1:numel(labels)
    labels{i} = ['{\it' strjoin(labels_split{i}([1 2]), ' ') '}'];
    if numel(labels_split{i}) > 2
        labels{i} = [labels{i} ' ' strjoin(labels_split{i}(3:end), ' ')];
    end
end

% figure

% read phylogenetic tree based on rbcl sequences and re-order the labels
tree = phytreeread(fullfile('..','..','data','rbcl_seqs','rbcl_tree.phylo.io.nwk'));
tax_sorted_labels = regexprep(get(tree, 'LeafNames'), '^\d+_', '');
tax_sorted_labels = strrep(tax_sorted_labels, 'Chlorella_vulgaris',...
    'Chlorella_vulgaris_UTEX_395');
tax_sorted_labels(ismember(tax_sorted_labels, setdiff(organisms_uniq,organisms_uniq(non_empty_idx)))) = [];
reorder_idx = cellfun(@(x)find(contains(organisms_uniq(non_empty_idx),x)), tax_sorted_labels);

subplot(1,2,2)
heatmap(ec_ji_pairwise(reorder_idx,reorder_idx),...
    'Colormap', summer, 'XDisplayLabels', labels(reorder_idx),...
    'YDisplayLabels', repmat({''}, size(labels)))

annotation('textarrow',[.926,.926],[0.95,0.95],'string','JI', ...
      'HeadStyle','none','LineStyle','none','HorizontalAlignment',...
      'center');
exportgraphics(gcf, '../../figures/co2_perc_and_pairwise_ji_co2_ec.png', 'Resolution', 300)

function mnx = getMNX(rxn_ids, source, db_table)

if strcmp(source, 'kegg')
    pfx = 'kegg.reaction:';
elseif strcmp(source, 'metacyc')
    pfx = 'metacyc.reaction:';
elseif strcmp(source, 'bigg')
    pfx = 'bigg.reaction:';
else
    error('Unknown source')
end

tmp = cellfun(@(x)...
    db_table.Var2(ismember(db_table.Var1, [pfx x])),...
    rxn_ids,...
    'un', 0);
mnx = [tmp{:}]';

end

function ec = getEC(mnx_ids, db_table)
tmp = cellfun(@(x)...
    db_table.Var4(ismember(db_table.Var1, x)),...
    mnx_ids);
tmp_split = cellfun(@(x)...
    strsplit(x, ';'),...
    tmp,...
    'un', 0);
ec = strtrim([tmp_split{:}]');
end

function writeEC(model, source, mnx_xref, mnx_prop, model_dir)
% write ID to EC translation results to file
mnx = getMNX(model.rxns, source, mnx_xref);
writetable(cell2table(mnx),...
    fullfile(model_dir, 'mnx_ids_resp.txt'),...
    'WriteVariableNames', false)

ec = getEC(mnx, mnx_prop);
writetable(cell2table(ec),...
    fullfile(model_dir, 'ec_resp.txt'),...
    'WriteVariableNames', false)

end
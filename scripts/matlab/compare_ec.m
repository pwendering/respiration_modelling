% compare the set of EC numbers that catalyze CO2-producing reactions
% (BioCyc, BRENDA) with the sets of EC numbers in the genome-scale metabolic models

top_model_dir = fullfile('..', '..', 'models');

tmp = load(fullfile(top_model_dir, 'models.mat'));
[models, model_publications, organisms_uniq] = deal(tmp.models, tmp.model_publications, tmp.organisms_uniq);

non_empty_idx = find(~cellfun(@isempty,models));

ec_nums = cell(numel(non_empty_idx), 1);
for i=non_empty_idx'
    model = models{i};
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
% if only reaction IDs are available, find associated EC numbers using
% MNXref
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
    if ~any(contains(cellstr(content), 'ec.txt'))
        
        model = models{model_idx};
        if strcmp(model_publications{model_idx}, 'Prigent2014')
            model.rxns = erase(model.rxnNames, 'metacyc:');
        else
            model.rxns = strtok(model.rxns, {'[', '_'});
        end
        writeEC(model, sources{i}, mnx_xref, mnx_prop, fullfile(top_model_dir, model_dirs{i}))
    end
    
    fid = fopen(fullfile('models', model_dirs{i}, 'ec.txt'));
    ec = textscan(fid, '%s');
    fclose(fid);
    
    ec_nums{model_idx} = [ec{:}];
end
clear mnx_prop mnx_xref all_dirs ec fid model_dirs model_idx content sources
save('ec_models', 'ec_nums')

%% Compare EC number sets with the set of CO2-producing reactions
tab = readtable(fullfile('..', '..', 'data', 'full-enzyme-data-table.tsv'),...
    'FileType', 'text',...
    'ReadVariableNames', true);
ec_ref = tab.ec_number(tab.is_viridiplantae==1);
clear tab

ec_ovlp = nan(size(non_empty_idx));
ec_perc = nan(size(non_empty_idx));

for i=non_empty_idx'
    if ~isempty(ec_nums{i})
        tmp_ec = regexp(ec_nums{i}, '(\d+\.){3}\d+', 'match');
        tmp_ec = [tmp_ec{:}];
        
        n_is = numel(intersect(tmp_ec, ec_ref));
        n_co2_prod = numel(ec_ref);
        
        ec_ovlp(i) = n_is;
        ec_perc(i) = n_is/n_co2_prod;
        
    end
end

%% plot results
% first subplot, the second one will be added by 'compare_ec_resp_models'
subplot(1,2,1)
bar(100*ec_perc(non_empty_idx),...
    'FaceColor', [.3 .4 .6],...
    'FaceAlpha', .8,...
    'EdgeColor', [.6 .6 .6],...
    'Horizontal', true)

y_lab = strrep(organisms_uniq(non_empty_idx), '_', ' ');

y_lab_split = cellfun(@strsplit, y_lab, 'un', 0);
for i=1:numel(y_lab)
    y_lab{i} = ['{\it' strjoin(y_lab_split{i}([1 2]), ' ') '}'];
    if numel(y_lab_split{i}) > 2
        y_lab{i} = [y_lab{i} ' ' strjoin(y_lab_split{i}(3:end), ' ')];
    end
end
yticks(1:numel(non_empty_idx))
yticklabels(y_lab)
xticks(0:20:60)
xlabel('Percentage')

set(gca,...
    'Box', 'off',...
    'TickLength',[0.01, 0.01])
%'FontSize', 14,...



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
    fullfile(model_dir, 'mnx_ids.txt'),...
    'WriteVariableNames', false)

ec = getEC(mnx, mnx_prop);
writetable(cell2table(ec),...
    fullfile(model_dir, 'ec.txt'),...
    'WriteVariableNames', false)

end
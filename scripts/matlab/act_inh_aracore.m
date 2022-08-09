% Map activators and inhibitors of CO2-producing reactions (BioCyc) to the
% Arabidopsis core model (Arnold and Nikoloski 2014)
top_model_dir = fullfile('..', '..', 'models');

% Read AraCore model
aracore = readCbModel(fullfile(top_model_dir, 'ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'ArabidopsisCoreModel.xml'));
inchi_keys = readtable(fullfile('..', '..', 'data', 'ArabidopsisCoreModel-inchikeys.tsv'),...
    'FileType', 'text');

% read enzyme parameter file
enz_tab = readtable(fullfile('..', '..', 'data', 'full-enzyme-data-table.tsv'),...
    'FileType', 'text');

% translation of BRENDA subtrate names to InChI keys
tmp = load(fullfile('..', '..', 'data', 'substrate-inchi-names.mat');
brenda_inchi = tmp.inchiKeys;
brenda_names = tmp.substrateNames;
clear tmp

% first check the overall overlap of the set of metabolites in AraCore with
% the set of activators/inhibitors
tmp = enz_tab.activator(~cellfun('isempty', enz_tab.activator));
tmp_split = cellfun(@(x)strsplit(x, '|'), tmp, 'un', 0);
act_set = unique([tmp_split{:}]');

tmp = enz_tab.inhibitor(~cellfun('isempty', enz_tab.inhibitor));
tmp_split = cellfun(@(x)strsplit(x, '|'), tmp, 'un', 0);
inh_set = unique([tmp_split{:}]');

clear tmp

act_inchi = cellfun(@(x)brenda_inchi(ismember(brenda_names, x)), act_set, 'un', 0);
act_inchi = [act_inchi{:}]';
act_inchi = act_inchi(~cellfun('isempty', act_inchi));

inh_inchi = cellfun(@(x)brenda_inchi(ismember(brenda_names, x)), inh_set, 'un', 0);
inh_inchi = [inh_inchi{:}]';
inh_inchi = inh_inchi(~cellfun('isempty', inh_inchi));

ovl_act = intersect(inchi_keys.inchikey, act_inchi);
ovl_inh = intersect(inchi_keys.inchikey, inh_inchi);

% translate inchi keys back to metabolite names
met_names_act = cellfun(@(x)unique(inchi_keys.name(ismember(inchi_keys.inchikey, x))), ovl_act);
met_names_inh = cellfun(@(x)unique(inchi_keys.name(ismember(inchi_keys.inchikey, x))), ovl_inh);

brenda_names_act = cellfun(@(x)...
    unique(regexprep(brenda_names(ismember(brenda_inchi, x)), '^\d+ ', '')),...
    ovl_act,...
    'un', 0);
brenda_names_act = vertcat(brenda_names_act{:});

brenda_names_inh = cellfun(@(x)...
    unique(regexprep(brenda_names(ismember(brenda_inchi, x)), '^\d+ ', '')),...
    ovl_inh,...
    'un', 0);
brenda_names_inh = vertcat(brenda_names_inh{:});

fprintf('Jaccard index inhibitors and activators in AraCore: %.2g\n\n',...
    numel(intersect(ovl_act, ovl_inh)) / numel(union(ovl_act, ovl_inh)))

only_act = setdiff(met_names_act, met_names_inh);
only_inh = setdiff(met_names_inh, met_names_act);

disp('Exclusively activating compounds:')
disp(only_act)

disp('Exclusively inhibiting compounds:')
disp(only_inh)
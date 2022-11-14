% Read AraCore model
clear;clc
aracore = readCbModel(fullfile('..','..','models','ArnoldNikoloski2014_Arabidopsis_thaliana',...
    'ArabidopsisCoreModel.xml'));
aracore = convertToIrreversible(aracore);
aracore = improveAraCoreMetNames(aracore);

n_rxns = numel(aracore.rxns);

ath_lineage = strtrim(strsplit(['Arabidopsis; Camelineae; Brassicaceae; Brassicales; malvids;'...
    'rosids; Pentapetalae; Gunneridae; eudicotyledons; Mesangiospermae;'...
    'Magnoliopsida; Spermatophyta; Euphyllophyta; Tracheophyta; Embryophyta;'...
    'Streptophytina; Streptophyta; Viridiplantae; Eukaryota; cellular organisms'], ';'));

% EC: exact matches => substrate matches => if yes ==> multiple? ==> take
% the one that shares the most of lineage with A. th, then take the maximum
% per reaction

%% in vitro turnover numbers
param_tab = readtable(fullfile('..','..','data','brenda_sabio_kcat_comb_lineage.tsv'),...
    'FileType', 'text');
kcats_vitro = zeros(n_rxns, 1);
confidence = zeros(n_rxns, 1);
for i=1:n_rxns
    if ~isempty(aracore.rxnECNumbers{i})
        
        % find EC numbers associated with reaction
        tmp_ec = strtrim(strsplit(aracore.rxnECNumbers{i}, '+'));
        tmp_kcats = zeros(numel(tmp_ec), 1);
        
        for j=1:numel(tmp_ec)
            % find matches for EC number in kcat table
            match_idx = find(ismember(param_tab.ec_number, tmp_ec{j}));
            tmp_subs = param_tab.substrate(match_idx);
            tmp_model_subs = aracore.metNames(aracore.S(:,i)<0);
            
            % match substrates
            substrate_matches = findSubstrateMatches(tmp_subs, tmp_model_subs);
            if any(substrate_matches)
                confidence(i) = confidence(i) + 2;
                match_idx_sub = match_idx(substrate_matches);
            else
                match_idx_sub = match_idx;
            end
            
            % match lineage
            if numel(match_idx_sub) > 1
                tmp_lin = cellfun(@(x)strtrim(strsplit(x, ';')),...
                    param_tab.lineage(match_idx_sub),...
                    'un', 0);
                lin_common = cellfun(@(x)numel(intersect(ath_lineage, x)), tmp_lin);
                match_idx_sub_lin = match_idx_sub(lin_common == max(lin_common));
                if max(lin_common) >= 18
                    confidence(i) = confidence(i) + 1;
                end
            else
                match_idx_sub_lin = [];
            end
            
            if numel(match_idx_sub_lin) > 0
                tmp_kcats(j) = max(max(param_tab.kcat(match_idx_sub_lin), 0));
            elseif numel(match_idx_sub) > 0
                tmp_kcats(j) = max(max(param_tab.kcat(match_idx_sub), 0));
            elseif numel(match_idx) > 0
                tmp_kcats(j) = max(max(param_tab.kcat(match_idx), 0));
            end
            
        end
        kcats_vitro(i) = max(tmp_kcats);
    end
end

rxns_vitro = strrep(aracore.rxns, '_f', '');
rxns_vitro = strrep(rxns_vitro, '_b', '_rev');

%% in vivo turnover numbers
tmp = load(fullfile('..', '..', 'data', 'max_vivo_kcats_kueken2020.mat'));
rxnMaxKcatTable = tmp.rxnMaxKcatTable;
kcats_vivo = rxnMaxKcatTable.k_cat_max;

%% find respiration reactions
resp_idx = cellfun(@(x)...
    contains(x, {'glycolysis', 'pentose phosphate pathway',...
    'tricarboxylic acid cycle', 'oxidative phosphorylation'}),...
    aracore.subSystems);
resp_rxns = rxns_vitro(resp_idx);
conf_resp = confidence(resp_idx);

%% match in vivo and in vivo data
idx_resp_vitro = cellfun(@(x)find(ismember(rxns_vitro, x)), resp_rxns);
idx_resp_vivo = cellfun(@(x)find(ismember(rxnMaxKcatTable.ReactionNames, x)), resp_rxns);

%% plot results
markers = {'^', 's', 'o', '*'};
conf_labels = {'unspecific', 'lineage', 'substrate', 'lineage + substrate'};
colors = lines(numel(conf_labels));
for i=1:numel(resp_rxns)
    marker = markers{conf_resp(i)+1};
    scatter(log10(kcats_vitro(idx_resp_vitro(i))), log10(kcats_vivo(idx_resp_vivo(i))), 80, ...
        'filled', ...
        'Marker', marker,...
        'MarkerFaceColor', colors(conf_resp(i)+1,:),...
        'MarkerEdgeColor', [.2 .2 .2])
    hold on
end

line([-3 7], [-3 7], 'Color', [.7 .7 .7], 'linestyle', '--')
ylim([-3 7])
xlim([-3 7])


h = zeros(numel(markers), 1);
for i=1:numel(markers)
    h(i) = scatter(NaN,NaN,'filled', ...
        'Marker', markers{i},...
        'MarkerFaceColor', colors(i,:),...
        'MarkerEdgeColor', [.2 .2 .2]);
end
hold off

legend(h, conf_labels,...
    'location', 'southeast',...
    'box', 'off',...
    'Fontsize', 14)

xlabel('log_{10} {\it in vitro} k_{cat} [s^{-1}]', 'FontSize', 24, 'FontName', 'Arial')
ylabel('log_{10} {\it in vivo} k_{cat} [s^{-1}]', 'FontSize', 24, 'FontName', 'Arial')

set(gca, 'Box' ,'on', 'LineWidth', 1.5,...
    'FontName', 'Arial', 'FontSize', 14)

% exportgraphics(gca, 'kcats_resp_vivo_vitro_aracore.png')

%% write result table
kcat_resp_tab = cell2table([...
    aracore.rxnNames(resp_idx), aracore.rxnECNumbers(resp_idx),...
    num2cell([kcats_vitro(idx_resp_vitro), kcats_vivo(idx_resp_vivo), conf_resp])],...
    'VariableNames', {'name', 'ec_number', 'in_vitro', 'in_vivo', 'in_vitro_confidence'},...
    'RowNames', aracore.rxns(resp_idx));

writetable(kcat_resp_tab, fullfile('..', '..', 'data', 'kcats_resp_vivo_vitro_aracore.csv'))
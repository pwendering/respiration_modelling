% read metabolic models and save them as Matlab workspace

top_model_dir = fullfile('..', '..', 'models');
model_dir_content = struct2table(dir(top_model_dir));
all_model_studies = model_dir_content.name(model_dir_content.isdir);
all_model_studies = setdiff(all_model_studies, {'.', '..'});

years = regexp(all_model_studies,'\d+','once','match');
organisms = regexp(all_model_studies,'_','split','once');
organisms = cellfun(@(x)x(2),organisms);
organisms_uniq = unique(organisms,'stable');

% select model for each organisms by latest publication date
model_dirs = repmat({''}, numel(organisms_uniq), 1);
model_publications = repmat({''},numel(organisms_uniq),1);
for i=1:numel(organisms_uniq)
    org_idx = find(ismember(organisms, organisms_uniq(i)));
    org_years = str2double(years(org_idx));
    [~,select_idx] = max(org_years);
    model_dirs(i) = all_model_studies(org_idx(select_idx));
    model_publications(i) = strtok(all_model_studies(org_idx(select_idx)),'_');
end

% read models from file
model_file_endings = '(xml)|(sbml)|(mat)|(xlsx)|(xls)';
models = cell(numel(organisms_uniq), 1);
for i=1:numel(model_dirs)
    fprintf('Processing %s metabolic model (%s)...\n', organisms_uniq{i}, model_publications{i})
    
    % trying to find a model file with correct file extension
    tmp_files = {dir(fullfile(top_model_dir, model_dirs{i})).name};
    model_file = regexp(tmp_files, ['(.+\.)(' model_file_endings ')'], 'match');
    non_empty_idx = ~cellfun(@isempty,model_file);
    
    if sum(non_empty_idx) == 1
        model_file = char(model_file{non_empty_idx});
        file_type = char(regexp(model_file, model_file_endings,'match'));
    elseif sum(non_empty_idx) > 1
        fprintf('More than one model file detected.\n')
        continue
    else
        model_file = [];
    end
    
    if ~isempty(model_file)
        
        % attempt to import the model using the COBRA toolbox
        try
            models{i} = readCbModel(fullfile(top_model_dir, model_dirs{i}, model_file));
            fprintf('Model read successfully (.%s).\n', file_type)
        catch ME
            fprintf('%s: %s\n', organisms_uniq{i}, ME.message)
        end
    else
        fprintf('No valid model file found.\n')
    end
    
    if isempty(models{i})
        
        % attempt to import the model using the RAVEN toolbox
        try
            models{i} = importModel(fullfile(top_model_dir, model_dirs{i}, model_file));
        catch ME
            fprintf('Model could also not be read using RAVEN importModel function.\n')
            disp(ME.message)
        end
        fprintf('Model read successfully using RAVEN import function.\n')
    end
end

save(fullfile(top_model_dir, 'models'), 'models', 'organisms_uniq', 'model_publications', '-v7.3')
clear
function cmf = getCarbonMolarFraction(model, umolFlag, scaleFlag)
%% cmf = getCarbonMolarFraction(model)
% calculate mmol(C) / gDW based on composition of biomass reaction
% INPUT
%   struct model:           metabolic model
%   logical umolFlag:       (false) whether biomass coefficients are given
%                           in umol/gDW and not mmol/gDW
%   logical scaleFlag:      (false) whether to divide the cmf value by the
%                           sum of weights of biomass components
% OUTPUT
%   double cmf:             molar fraction of carbon in biomass

if nargin < 2 || isempty(umolFlag)
    umolFlag = false;
end

if nargin < 3 || isempty(scaleFlag)
    scaleFlag = false;
end
if sum(model.c) > 1
    error('More than one objective')
elseif sum(model.c) == 0
    error('No objective set')
elseif any(model.c < 0)
    error('Negative objective detected')
end
    
% biomass reaction index
bio_idx = find(model.c==1);
% metabolite indices of biomass reaction substrates
met_idx = find(model.S(:, bio_idx) < 0);
% metabolite formulas
met_form = model.metFormulas(met_idx);
empty_idx = cellfun(@isempty, met_form);
met_form(empty_idx) = {'X'};
% molecular masses [g/umol]
if umolFlag
    met_mw = getMolecularMass(met_form) / 1e6;
else
    met_mw = getMolecularMass(met_form) / 1e3;
end
% coefficients
met_coeff = abs(model.S(met_idx, bio_idx));
% calculate w/w of carbon atoms per molecule
pct_carbon = zeros(size(met_idx));
for i = 1:numel(met_idx)
    if ismember('C', met_form{i})
        [Ematrix, elements] = getElementalComposition(met_form{i});
        element_mw = getMolecularMass(elements);
        element_mass = Ematrix' .* element_mw;
        pct_carbon(i) = element_mass(ismember(elements, 'C')) / sum(element_mass);
    end
end
% calculate w/w of metabolites [g/gDW]
met_weights = met_coeff .* met_mw;
% calculate carbon weight fraction of metabolites
carbon_w_w = met_weights .* pct_carbon;
% calculate mmol carbon per gram dry weight
cmf = 1000 * sum(carbon_w_w) / getMolecularMass('C');

if scaleFlag
    % divide cmf by the sum of weights of biomass components (should be 1 g/gDW
    % but this does not always add up)
    cmf = cmf / sum(met_weights);
end
end


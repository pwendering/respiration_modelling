function model = improveAraCoreMetNames(model)
model.metNames(ismember(model.metNames,'CO2, carbon dioxide')) = ...
    {'CO2'};
model.metNames(ismember(model.metNames,'O2, oxygen')) = ...
    {'O2'};
model.metNames(ismember(model.metNames,'Ribulose 1,5-bisphosphate')) = ...
    {'D-ribulose 1,5-bisphosphate'};
model.metNames(ismember(model.metNames,'Histidinol phosphate')) = ...
    {'L-histidinol 1-phosphate'};
model.metNames(ismember(model.metNames,'3-Isopropylmalate')) = ...
    {'(2R,3S)-3-isopropylmalate'};
model.metNames(ismember(model.metNames,'Aspartate 4-semialdehyde')) = ...
    {'L-Aspartate 4-semialdehyde'};
model.metNames(ismember(model.metNames,'Fructose 1,6-bisphosphate')) = ...
    {'D-Fructose 1,6-bisphosphate'};
model.metNames(ismember(model.metNames,'3-Phosphoglycerate')) = ...
    {'(R)-3-phosphoglycerate'};
model.metNames(ismember(model.metNames,'Glucose 6-phosphate')) = ...
    {'D-Glucose 6-phosphate'};
model.metNames(ismember(model.metNames,'Glucose 1-phosphate')) = ...
    {'D-Glucose 1-phosphate'};
model.metNames(ismember(model.metNames,'Glucose')) = ...
    {'D-Glucose'};
model.metNames(ismember(model.metNames,'Fructose')) = ...
    {'D-Fructose'};
model.metNames(ismember(model.metNames,'Fructose 6-phosphate')) = ...
    {'D-Fructose 6-phosphate'};
model.metNames(ismember(model.metNames,'Fructose 2,6-bisphosphate')) = ...
    {'D-Fructose 2,6-bisphosphate'};
model.metNames(ismember(model.metNames,'Nicotinamide adenine dinucleotide')) = ...
    {'NAD+'};
model.metNames(ismember(model.metNames,'Nicotinamide adenine dinucleotide - reduced')) = ...
    {'NADH'};
model.metNames(ismember(model.metNames,'Nicotinamide adenine dinucleotide phosphate')) = ...
    {'NADP+'};
model.metNames(ismember(model.metNames,'Nicotinamide adenine dinucleotide phosphate - reduced')) = ...
    {'NADPH'};
model.metNames(ismember(model.metNames,'Ribose')) = ...
    {'D-Ribose'};
model.metNames(ismember(model.metNames,'Serine')) = ...
    {'L-Serine'};
model.metNames(ismember(model.metNames,'Glutamate')) = ...
    {'L-Glutamate'};
model.metNames(ismember(model.metNames,'Alanine')) = ...
    {'L-Alanine'};
model.metNames(ismember(model.metNames,'Glutamine')) = ...
    {'L-Glutamine'};
model.metNames(ismember(model.metNames,'Arginine')) = ...
    {'L-Arginine'};
model.metNames(ismember(model.metNames,'Aspartate')) = ...
    {'L-Aspartate'};
model.metNames(ismember(model.metNames,'Asparagine')) = ...
    {'L-Asparagine'};
model.metNames(ismember(model.metNames,'Cysteine')) = ...
    {'L-Cysteine'};
model.metNames(ismember(model.metNames,'Threonine')) = ...
    {'L-Threonine'};
model.metNames(ismember(model.metNames,'Histidine')) = ...
    {'L-Histidine'};
model.metNames(ismember(model.metNames,'Isoleucine')) = ...
    {'L-Isoleucine'};
model.metNames(ismember(model.metNames,'Leucine')) = ...
    {'L-Leucine'};
model.metNames(ismember(model.metNames,'Lysine')) = ...
    {'L-Lysine'};
model.metNames(ismember(model.metNames,'Homocysteine')) = ...
    {'L-Homocysteine'};
model.metNames(ismember(model.metNames,'Methionine')) = ...
    {'L-Methionine'};
model.metNames(ismember(model.metNames,'Arogenate')) = ...
    {'L-Arogenate'};
model.metNames(ismember(model.metNames,'Ornithine')) = ...
    {'L-Ornithine'};
model.metNames(ismember(model.metNames,'Proline')) = ...
    {'L-Proline'};
model.metNames(ismember(model.metNames,'Tryptophan')) = ...
    {'L-Tryptophan'};
model.metNames(ismember(model.metNames,'Valine')) = ...
    {'L-Valine'};
model.metNames(ismember(model.metNames,'Tyrosine')) = ...
    {'L-Tyrosine'};
model.metNames(ismember(model.metNames,'Inosine 5-phosphate')) = ...
    {'Inosine 5''phosphate'};
model.metNames(ismember(model.metNames,'Histidine')) = ...
    {'L-Histidine'};
model.metNames(ismember(model.metNames,'Phenylalanine')) = ...
    {'L-Phenylalanine'};
end
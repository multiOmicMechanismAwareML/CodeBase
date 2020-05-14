
%initCobraToolbox
load('iSce926.mat');
load('reaction_expression.mat'); % Reaction expressions
load('pos_genes_in_react_expr.mat'); % Position of the genes in the reaction expression var
load('idxs_genes_sorted_by_length.mat');

dataset = 'main'; % Select dataset ('main' or 'ITS' for the experimentally-independent test set)

% Load gene expression and growth data
if strcmp(dataset, 'main')
    load('MSB_data_expression.mat'); 
    load('MSB_data_growth_rates.mat');
    expr_data = msbdata;
    growth_rates = growthRatesMSB;
elseif strcmp(dataset, 'ITS')
    load('independent_test_set_expression.mat');
    load('independent_test_set_growth_rates.mat');
    expr_data = delmutantslimmaSameithetal2015;
    growth_rates = growthRatesSameithetal2015;
end

growthReactName = 'r_4041'; % yeast 8 biomass pseudoreaction
model = changeObjective(model, growthReactName); % Set the objetive to be the biomass

genes = model.genes;
genes_in_dataset = expr_data.geneName;
GeneExpressionArray = ones(numel(genes),1);
gamma = 1;
lowerPercentBounded = 1; % We initially explored different low bound values and settled on 1
bounds = -50;
model = setMediaConditions(model, bounds);
fluxes = model.rxns;

for t=3:(width(expr_data))
    expr_profile = table2array(expr_data(:,t));
    if contains(expr_data.Properties.VariableNames{t}, '_') % Double KO case
        deletedGene = strsplit(expr_data.Properties.VariableNames{t}, '_');
    else % Single KO case
        deletedGene = expr_data.Properties.VariableNames{t};
    end
    deletedGeneIndex = find(strcmp(expr_data.commonName,deletedGene),1);
    deletedGeneExpressionIndexInArray = inf;
    pos_genes_in_dataset = zeros(numel(genes),1);
    for i=1:numel(genes)
        position = find(strcmp(genes{i},genes_in_dataset));
        if ~isempty(position)
            pos_genes_in_dataset(i) = find(strcmp(genes{i},genes_in_dataset));
            if pos_genes_in_dataset(i) == deletedGeneIndex
                deletedGeneExpressionIndexInArray = i;
            end
            GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));
        end
    end

    NumberOfGenes = numel(genes);

    [v_out, f_out] = evaluate_objective(gamma,GeneExpressionArray,model,genes,reaction_expression,pos_genes_in_react_expr,idxs_genes_sorted_by_length,...
        lowerPercentBounded,deletedGeneExpressionIndexInArray);
    f_out
    
    if f_out > 0 
        fluxes = [fluxes, array2table(v_out,'VariableNames',{expr_data.Properties.VariableNames{t}})];
    end
end

% Organise the data and plot the correlated results 
justGrowth =  fluxes(strcmp(model.rxns,growthReactName),2:end);
justGrowthColumn = table2array(justGrowth);
justGrowthColumn = array2table(justGrowthColumn.');
justGrowthColumn = [justGrowth.Properties.VariableNames', justGrowthColumn];
justGrowthColumn.Properties.VariableNames{1} = 'commonName';
growth_rates.commonName = cellstr(growth_rates.commonName);
joined = innerjoin(justGrowthColumn, growth_rates);
joinedForComp = table2array(joined(:,2:end));
scatter(joinedForComp(:,1),joinedForComp(:,2));
disp(['correlation between FBA growth rate and doubling time fold change: ', num2str(corr(joinedForComp(:,1),joinedForComp(:,2)))])
xlabel('FBA growth rate');
ylabel ('Doubling time fold change log2(strain/WT)');
p = polyfit(joinedForComp(:,1),joinedForComp(:,2),1);
f = polyval(p,joinedForComp(:,1));
hold on
plot(joinedForComp(:,1),f,'--r', 'LineWidth', 1)

% Merge data and save
growth_rates_array = table2array(growth_rates(:, 2));
growth_rates_t = array2table(growth_rates_array');
growth_rates.commonName = cellfun(@strrep, growth_rates.commonName, repmat({'-'}, length(growth_rates.commonName), 1), ...
    repmat({''}, length(growth_rates.commonName), 1), 'UniformOutput', false);
growth_rates.commonName = cellfun(@strrep, growth_rates.commonName, repmat({','}, length(growth_rates.commonName), 1), ...
    repmat({''}, length(growth_rates.commonName), 1), 'UniformOutput', false);
growth_rates.commonName = cellfun(@strrep, growth_rates.commonName, repmat({'('}, length(growth_rates.commonName), 1), ...
    repmat({''}, length(growth_rates.commonName), 1), 'UniformOutput', false);
growth_rates.commonName = cellfun(@strrep, growth_rates.commonName, repmat({')'}, length(growth_rates.commonName), 1), ...
    repmat({''}, length(growth_rates.commonName), 1), 'UniformOutput', false);
growth_rates_t.Properties.VariableNames = cellstr(growth_rates.commonName);
growth_rates_t.geneName = {'log2relT'};
fluxes.Properties.VariableNames{1} = 'geneName';

[a, b] = ismember(growth_rates_t.Properties.VariableNames, fluxes.Properties.VariableNames);
t = outerjoin(growth_rates_t(:,a), fluxes(:,b(a)), 'MergeKeys', true);
[a, b] = ismember(t.Properties.VariableNames, expr_data.Properties.VariableNames);
t = outerjoin(t(:,a), expr_data(:,b(a)), 'MergeKeys', true);
t = sortrows(t, 'geneName');
t.geneName = strrep(t.geneName, '-', '_');

if strcmp(dataset, 'main')
    t_array = table2array(t(:, 1:end-1));
    completeDataset = array2table(t_array');
    completeDataset.Properties.VariableNames = cellstr(t.geneName);
    completeDataset.Properties.RowNames = t.Properties.VariableNames(1:end-1);
    rxn_idxs = contains(completeDataset.Properties.VariableNames, 'r_');
    target_idx = strcmp(completeDataset.Properties.VariableNames, 'log2relT');
    expr_idxs = and(not(rxn_idxs), not(target_idx));
    completeDataset = completeDataset(:, [find(target_idx), find(rxn_idxs), find(expr_idxs)]);
    writetable(completeDataset, 'completeDataset.csv', 'WriteRowNames', true);
elseif strcmp(dataset, 'ITS')
    t_array = table2array(t(:, 1:end-1));
    completeDataset_ITS = array2table(t_array');
    completeDataset_ITS.Properties.VariableNames = cellstr(t.geneName);
    completeDataset_ITS.Properties.RowNames = t.Properties.VariableNames(1:end-1);
    rxn_idxs = contains(completeDataset_ITS.Properties.VariableNames, 'r_');
    target_idx = strcmp(completeDataset_ITS.Properties.VariableNames, 'log2relT');
    expr_idxs = and(not(rxn_idxs), not(target_idx));
    completeDataset_ITS = completeDataset_ITS(:, [find(target_idx), find(rxn_idxs), find(expr_idxs)]);
    writetable(completeDataset_ITS, 'completeDataset_ITS.csv', 'WriteRowNames', true);
end



%%
% Setting the initial media conditions, adjusting the uptake rates
function model = setMediaConditions(model,bound)
exchangeReactions = {'ammonium exchange' ...
'sulphate exchange' ...
'biotin exchange' ...
'(R)-pantothenate exchange' ...
'folic acid exchange' ...
'myo-inositol exchange' ...
'nicotinate exchange' ...
'4-aminobenzoate exchange' ...
'pyridoxine exchange' ...
'H+ exchange' ...
'riboflavin exchange' ...
'thiamine(1+) exchange' ...
'sulphate exchange' ...
'potassium exchange' ...
'phosphate exchange' ...
'sulphate exchange' ...
'sodium exchange' ...
'L-alanine exchange' ...
'L-arginine exchange' ...
'L-asparagine exchange' ...
'L-aspartate exchange' ...
'L-cysteine exchange' ...
'L-glutamate exchange' ...
'L-glutamine exchange' ...
'glycine exchange' ...
'L-histidine exchange' ...
'L-isoleucine exchange' ...
'L-leucine exchange' ...
'L-lysine exchange' ...
'L-methionine exchange' ...
'L-phenylalanine exchange' ...
'L-proline exchange' ...
'L-serine exchange' ...
'L-threonine exchange' ...
'L-tryptophan exchange' ...
'L-tyrosine exchange' ...
'L-valine exchange' ...
'oxygen exchange' ...
'adenine exchange' ...
'uracil exchange'};
for i = 1 : numel(exchangeReactions)
     model.lb(model.rxns{find(ismember(model.rxnNames, exchangeReactions{i}),1)}) = bound;
end
end
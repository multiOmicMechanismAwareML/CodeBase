function [v1, f_out] = evaluate_objective(gammaVal, x,NumberOfObjectives,NumberOfGenes,fbarecon,geni,reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length,lowbound,deletedGeneIndex)

yt=x';

eval_reaction_expression = reaction_expression;


reactionIndexesWithDeletedGene = false;
%indices are sorted by length of their string, so the longest get replaced first. This avoids, for example, that if we have two genes 123.1 and 23.1, the substring '23.1' is replaced in both. If we first replace the longest and finally the shortest, this problem is avoided
for i=ixs_geni_sorted_by_length %loop over the array of the non-1 gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values. All the gene set expressions with only 1s as gene values , will be left empty and at the end of this loop everything empty will be substituted with 1 anyway. This avoids looping over all the genes yt, which is very expensive
    posizioni_gene = pos_genes_in_react_expr{i};
    if i == deletedGeneIndex
        reactionIndexesWithDeletedGene = posizioni_gene;
    end
    for j=1:length(posizioni_gene)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
        eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, geni{i}, num2str(yt(i),'%.15f'));  %Matlab strangely truncates decimal digits when using num2str. Addimg %.12f at least ensures that 12 decimal digits are included in the number converted into string
    end
end
eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally


% Eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions or reactions whose genes have all gene expression 1) with 1, i.e. gene expressed nomally

num_reaction_expression = zeros(1,length(eval_reaction_expression));
gamma = zeros(1,length(reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    
    num_parenthesis = numel(strfind(str,')'));
    while (num_parenthesis > 32) %if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
        substrings_to_replace = regexp(str, to_replace, 'match');
        if isempty(substrings_to_replace)
            num_parenthesis = 0; %if num_parenthesis > 32 and there is nothing that can be replaced with regexp, we force this, in order to avoid an endless loop. Later, eval will catch an exception as it cannot evaluate when num_parenthesis>32
        else
            for j = 1:numel(substrings_to_replace)
                ss_rep = substrings_to_replace{j};
                str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
            end
            num_parenthesis = numel(strfind(str,')'));
        end
    end
    
    str = regexprep(str,'/','');
    
    try
        num_reaction_expression(i) = eval(str);   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
    catch
        %i
        %warning('Problem using function.  Assigning a value of 1.');
        num_reaction_expression(i) = 1;
    end
end

% Try different values for these -> [-3 3]
gamma = gammaVal *  ones(1,length(reaction_expression));


if reactionIndexesWithDeletedGene
    reactionsWithDeletedGene = num_reaction_expression(reactionIndexesWithDeletedGene);
    num_reaction_expression(reactionIndexesWithDeletedGene) = nan;
end
reaction_expressionUpCorrected = prctile(num_reaction_expression,99);
num_reaction_expression(num_reaction_expression > reaction_expressionUpCorrected) = reaction_expressionUpCorrected;
if lowbound > 0
    reaction_expressionLowCorrected = prctile(num_reaction_expression, lowbound);
    num_reaction_expression(num_reaction_expression > reaction_expressionLowCorrected) = reaction_expressionLowCorrected;
end
if reactionIndexesWithDeletedGene
    num_reaction_expression(reactionIndexesWithDeletedGene) = reactionsWithDeletedGene;
end




for i=1:length(num_reaction_expression)   % Loop over the array of the geneset expressions 
    fbarecon.lb(i) = fbarecon.lb(i)*(num_reaction_expression(i)^gamma(i));
    fbarecon.ub(i) = fbarecon.ub(i)*(num_reaction_expression(i)^gamma(i)); 
end


% Min-norm run through
[solution] = optimizeCbModel(fbarecon, 'max',  'one', 1e-6);
f_out(1) = solution.f;
f_out(2) = solution.f;
try
    v1 = solution.v;
catch
    v1 = [0 0];
end



format longG; format compact;
f_out;

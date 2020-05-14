function [v_out, f_out] = evaluate_objective(gammaVal, x, fbarecon, genes, reaction_expression, pos_genes_in_react_expr, idxs_genes_sorted_by_length, lowbound, deletedGeneIndex)

x=x';     
eval_reaction_expression = reaction_expression;

reactionIndexesWithDeletedGene = false;
% indices are sorted by length of their string, so the longest get replaced first. This avoids, for example, that if we have two genes 123.1 and 23.1, the substring '23.1' is replaced in both. If we first replace the longest and finally the shortest, this problem is avoided
for i=idxs_genes_sorted_by_length
    gene_position = pos_genes_in_react_expr{i};
    if i == deletedGeneIndex
        reactionIndexesWithDeletedGene = gene_position;
    end
    for j=1:length(gene_position) % for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
        eval_reaction_expression{gene_position(j)} = strrep(eval_reaction_expression{gene_position(j)}, genes{i}, num2str(x(i),'%.15f'));
    end
end
eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally


num_reaction_expression = zeros(1,length(eval_reaction_expression));
gamma = zeros(1,length(reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};

    num_parenthesis = numel(strfind(str,')'));
    while (num_parenthesis > 32) % if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string
        % search for all the strings of the kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';
        substrings_to_replace = regexp(str, to_replace, 'match');
        if isempty(substrings_to_replace)
            num_parenthesis = 0; % if num_parenthesis > 32 and there is nothing that can be replaced with regexp, we force this, in order to avoid an endless loop. Later, eval will catch an exception as it cannot evaluate when num_parenthesis>32
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
        num_reaction_expression(i) = eval(str); % evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
     catch
         num_reaction_expression(i) = 1;
     end
end

gamma = gammaVal *  ones(1,length(reaction_expression));

% Winsorisation of num_reaction_expression, excluding deleted gene
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

% Alter reaction bounds
for i=1:length(num_reaction_expression)
    fbarecon.lb(i) = fbarecon.lb(i)*(num_reaction_expression(i)^gamma(i));
    fbarecon.ub(i) = fbarecon.ub(i)*(num_reaction_expression(i)^gamma(i));
end

% pFBA
[solution] = optimizeCbModel(fbarecon, 'max', 'one');
f_out = solution.f;
try
    v_out = solution.v;
catch
    v_out = zeros(length(fbarecon.rxns), 1);
end

format longG; format compact;
end

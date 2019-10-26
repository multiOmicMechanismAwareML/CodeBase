
function expFBA(pop,gen,model_name,M,varargin)

if nargin >= 4
    last_gen = varargin{1};
    num_cores = varargin{2};
else
    last_gen = 0;
    num_cores = 1;
end


%last_gen allows to start the evolution from the last population.mat file
%generated, by indicating its number as fourth argument. If no argument is
%indicated, than the function will start the evolution from scratch

% if (num_cores > 1)
%     %parallel computing toolbox matlab
%     if (matlabpool('size') == 0)  %opens only if it is closed 
%         matlabpool('open','local',num_cores)
%     end
% end


geni_aux = load('geni.mat');
geni = geni_aux.geni;
reaction_expression_aux = load('reaction_expression.mat');
reaction_expression = reaction_expression_aux.reaction_expression;


% load model
fbarecon_aux = load([model_name '.mat']);
fbarecon = fbarecon_aux.fbarecon;


%V=fbarecon.nbin;
V = length(geni); %length of the input individuals without ranking, crowding distance and outputs
N=pop;



min_range = zeros(1,V);     %the expression of each gene is >= 0
max_range = 100*ones(1,V);  %the expression of each gene is <= 100


%% Initialize the population

if (last_gen==0)
    chromosome=ones(pop,V+M);   %gene expressions initialized as 1, i.e. all the genes are normally expressed (reference state). A chromosome means, in our case, an array of gene expressions
    
%    chromosome = chromosome + 0.1.*(rand(pop,V+M)-0.5.*ones(pop,V+M));  %adds some initial random noise
%     chromosome = chromosome + 2.*(rand(pop,V+M)-0.5.*ones(pop,V+M));   %adds some initial random noise, random in every cell
     chromosome = chromosome + ( 2.*rand(pop,1)*ones(1,V+M) - ones(pop,V+M) );   %adds some initial random noise, random and different in every row, but the same noise for every cell in the same row
    %[v1, fmax, fmax_max] = flux_balance_trilevel(fbarecon,true);
%     [v1, fmax, fmin] = flux_balance(fbarecon,true);
%     chromosome(i,V+1)= - fbarecon.f' * v1;% biomass
%     chromosome(i,V+2)= fmin;% 2nd obj sum of fluxes of IDH1 and IDH2 reactions (that are irreversible, so lb=0 so we can use fbarecon.g' * v1 to get |v(3065)|+|v(3066)|, because fbarecon.g has 1 in both positions 3065 and 3066
%     chromosome(i,V+3)= sum(abs(chromosome(i,:)-ones(1,numel(chromosome(i,:))))); % 3rd obj: minimize the gene overexpression and underexpression (which is the distance from 1)

    aux = zeros(pop,M);
    parfor i=1:pop
        aux(i,:) = evaluate_objective(chromosome(i,1:V), M, V, fbarecon, geni, reaction_expression, pos_genes_in_react_expr);
    end
    chromosome(:,V+1:V+M) = aux(:,:);
    %% Sort the initialized population
    % Sort the population using non-domination-sort. This returns two columns
    % for each individual which are the rank and the crowding distance
    % corresponding to their position in the front they belong. At this stage
    % the rank and the crowding distance for each chromosome is added to the
    % chromosome vector for easy of computation.
    chromosome = non_domination_sort_mod(chromosome, M, V);
    population=['populations_' model_name '/population' num2str(1)];
    save(population, 'chromosome');
    last_gen = 1;
else
    sol=['populations_' model_name '/population' num2str(last_gen) '.mat'];
    load(sol);
end

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.
%diary('elapsed_time')
tic;
for i = last_gen+1 : gen
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
   
    pool = round(pop/2);
    tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selcted as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    parent_chromosome = tournament_selection(chromosome, pool, tour);

    % Perfrom crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
%     if (num_cores>1)
%         offspring_chromosome = ...
%         genetic_operator_parallel(parent_chromosome, ...
%         M, V, mu, mum, min_range, max_range, fbarecon, geni, reaction_expression);
%     else
        offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range, fbarecon, geni, reaction_expression, pos_genes_in_react_expr);
%     end
    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    % temp is a dummy variable.
    clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
%   chromosome = delete_redundant(chromosome,fbarecon);
    population=['populations_' model_name '/population' num2str(i)]; 
    save(population, 'chromosome');
    if ~mod(i,10)
        clc
        fprintf('%d generations completed\n',i);
    end
end
toc;
%diary off
%matlabpool close
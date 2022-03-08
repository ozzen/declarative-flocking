rows = 3;
cols = 2;
%% Fitness
i = 1;
subplot(rows,2,i);
plotFileData(sprintf([prefix '/gmpcOBST_%d_fitness.txt'], simNum), 'Fitness Value' ,Num);

%%  CBF fitness value 
i = 3;
subplot(rows,2,i)
res = BFbasedFitness( x, y, vx, vy, dmin, maxAcc, LDfname );
plotCurve(res, 'Time', 'LD fitness', 1.5, Num)

%% Cohesion Fitness term
i = 5;
subplot(rows,2,i)


%% LD
i = 2;
subplot(rows,2,i);
LDfname = sprintf([prefix '/gmpcOBST_%d_LD.txt'], simNum);
plotFileData(LDfname, 'LD' ,Num);

%% BF
i = 4;
subplot(rows,2,i);
PlotBf;

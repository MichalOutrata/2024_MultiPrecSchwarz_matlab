clear; clc; close('all')
% addpath(genpath('./RankTruncData'))
% MyData = {'Schwarz methods as preconds','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'theta for dAS', dampingTheta,'# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, ...
%     'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
%     'GMRES prec type (noprec/L/R)', GMRES_PrecType, 'max # GMRES iters', GMRES_nmbiter, 'GMRES conv. tolerance', GMRES_relresacc, ...
%     'GMRES ConvCrvs noprec', GMRES_resnorms_noprec, 'GMRES # iter no prec', GMRES_nmbittoconv_noprec, ...
%     'GMRES ConvCrvs dAS', GMRES_resnorms_dAS, 'GMRES # iter dAS', GMRES_nmbittoconv_dAS, ... 
%     'GMRES ConvCrvs RAS', GMRES_resnorms_RAS, 'GMRES # iter RAS', GMRES_nmbittoconv_RAS, ... 
%     'GMRES ConvCrvs MS', GMRES_resnorms_MS, 'GMRES # iter MS', GMRES_nmbittoconv_MS};

ProblemChoice = 1;
s_LoadString = append(append('SavedData_ClassicmpSM_Stieltjes_AsPrecs_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{8}; list_of_nmb_int_gridcols = LoadedData{14}; GMRES_nmbiter_Prblm1 = LoadedData{18};
GMRES_resnorms_noprec_Prblm1 = LoadedData{22}; GMRES_nmbittoconv_noprec_Prblm1 = LoadedData{24};
GMRES_resnorms_dAS_Prblm1 = LoadedData{26}; GMRES_nmbittoconv_dAS_Prblm1 = LoadedData{28};
GMRES_resnorms_RAS_Prblm1 = LoadedData{30}; GMRES_nmbittoconv_RAS_Prblm1 = LoadedData{32};
GMRES_resnorms_MS_Prblm1 = LoadedData{34}; GMRES_nmbittoconv_MS_Prblm1 = LoadedData{36};


ProblemChoice = 2;
s_LoadString = append(append('SavedData_ClassicmpSM_Stieltjes_AsPrecs_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{8}; list_of_nmb_int_gridcols = LoadedData{14}; GMRES_nmbiter_Prblm2 = LoadedData{18};
GMRES_resnorms_noprec_Prblm2 = LoadedData{22}; GMRES_nmbittoconv_noprec_Prblm2 = LoadedData{24};
GMRES_resnorms_dAS_Prblm2 = LoadedData{26}; GMRES_nmbittoconv_dAS_Prblm2 = LoadedData{28};
GMRES_resnorms_RAS_Prblm2 = LoadedData{30}; GMRES_nmbittoconv_RAS_Prblm2 = LoadedData{32};
GMRES_resnorms_MS_Prblm2 = LoadedData{34}; GMRES_nmbittoconv_MS_Prblm2 = LoadedData{36};






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'o','d','s','^','*','>'}; MyLines = {'--','-.',':','-'};

Cbar_Ticks = []; Cbar_TickLabels = {};
for ind = 1:length(list_of_nmb_int_gridcols)
    if mod(ind-1,3) == 0
    Cbar_Ticks = [Cbar_Ticks,ind];
    Cbar_TickLabels = {Cbar_TickLabels{:}, int2str( list_of_nmb_int_gridcols(ind)^2 ) };
    end
end

nmb_iters_ToPlot = GMRES_nmbiter_Prblm1; iters_mesh = 0:nmb_iters_ToPlot-1; nmb_DigsToPlot = 6;





%%% plot interaction of "d_s " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = 2; nmbColsTiles = 3;
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','loose','Padding','Compact');

ind_MeshSizeToPlot = length(list_of_nmb_int_gridcols); angle1 = 227; angle2 = 26;


% first row: convergence curves of dAS/RAS/MS-prec GMRES for different d_s
%%%%%
nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_dAS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig};
    semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('damped additive Schwarz','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_RAS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('restricted additive Schwarz','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_MS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2);
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('multiplicative Schwarz','FontSize',24, 'interpreter', 'latex');


% legend
plt1 = plot(nan,nan,'Color','black','LineWidth',6); hold on;
LegendHandles(1) = plt1; LegendLabels{1} = 'no precond.';
for ind_nmbdig = nmb_DigsToPlot:-1:1
    plt1 = plot(nan,nan,'Color',[MyColors{ind_nmbdig}],'LineWidth',6); hold on;
    LegendHandles(1+ind_nmbdig) = plt1; LegendLabels{1+ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
    if ind_nmbdig == nmb_DigsToPlot, LegendLabels{1+ind_nmbdig} = append(append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig))),' $\quad$ '); end
end
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','d','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+1} = ' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$'; LegendHandles(1+nmb_DigsToPlot+1) = plt1;
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','^','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+2} = '$\; \mathrm{RAS} \quad$'; LegendHandles(1+nmb_DigsToPlot+2) = plt1;
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','o','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+3} = '$\; \mathrm{MS} \quad$'; LegendHandles(1+nmb_DigsToPlot+3) = plt1;
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 11; Lgnd.Layout.Tile = 'North';


% second row: are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols
%%%%%
nexttile();
ribbon(GMRES_nmbittoconv_dAS_Prblm1'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); 

nexttile();
ribbon(GMRES_nmbittoconv_RAS_Prblm1'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

nexttile();
ribbon(GMRES_nmbittoconv_MS_Prblm1'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

cbar = colorbar(); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
% cbar.Layout.Tile = 'east';



% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontSize',24);
end














%%% plot interaction of "d_s " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = 3; nmbColsTiles = 3;
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

ind_MeshSizeToPlot = length(list_of_nmb_int_gridcols); angle1 = 227; angle2 = 26;


% first row: convergence curves of dAS/RAS/MS-prec GMRES for different d_s - Prblm 1
%%%%%
nexttile(1);
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_dAS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig};
    semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('damped additive Schwarz','FontSize',24,'FontWeight','bold', 'interpreter', 'latex');

nexttile(2);
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_RAS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('restricted additive Schwarz','FontSize',24,'FontWeight','bold','interpreter', 'latex'); %title('GMRES preconditioned with RAS','FontSize',24, 'interpreter', 'latex');

nexttile(3);
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_MS_Prblm1{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec_Prblm1{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData(1:nmb_iters_ToPlot),'Color','k','Marker','o','MarkerSize',12,'LineWidth',2);
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('multiplicative Schwarz','FontSize',24,'FontWeight','bold','interpreter', 'latex'); %title('GMRES preconditioned with MS','FontSize',24, 'interpreter', 'latex');


% legend
plt1 = plot(nan,nan,'Color','black','LineWidth',6); hold on;
LegendHandles(1) = plt1; LegendLabels{1} = 'no precond.';
for ind_nmbdig = nmb_DigsToPlot:-1:1
    plt1 = plot(nan,nan,'Color',[MyColors{ind_nmbdig}],'LineWidth',6); hold on;
    LegendHandles(1+ind_nmbdig) = plt1; LegendLabels{1+ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
    % if ind_nmbdig == nmb_DigsToPlot, LegendLabels{1+ind_nmbdig} = append(append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig))),' $\quad$ '); end
end
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','d','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+1} = ' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$'; LegendHandles(1+nmb_DigsToPlot+1) = plt1;
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','^','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+2} = '$\; \mathrm{RAS} \quad$'; LegendHandles(1+nmb_DigsToPlot+2) = plt1;
% plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','o','MarkerSize',16,'LineWidth',2); hold on;
% LegendLabels{1+nmb_DigsToPlot+3} = '$\; \mathrm{MS} \quad$'; LegendHandles(1+nmb_DigsToPlot+3) = plt1;
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 11; Lgnd.Layout.Tile = 'North';


% third row: are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols - Prblm 3
%%%%%
nexttile(4);
ribbon(GMRES_nmbittoconv_dAS_Prblm1'); view(angle1,angle2);
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({});  zlabel('Problem 4','FontSize',24, 'interpreter', 'latex');

nexttile(5);
ribbon(GMRES_nmbittoconv_RAS_Prblm1'); view(angle1,angle2);
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({}); zticklabels({}); 

nexttile(6);
ribbon(GMRES_nmbittoconv_MS_Prblm1'); view(angle1,angle2);
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({}); zticklabels({}); 

% cbar = colorbar(); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
% cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
% cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
% cbar.Layout.Tile = 'east';


% third row: are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols - Prblm 4
%%%%%
nexttile(7);
ribbon(GMRES_nmbittoconv_dAS_Prblm2'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0); zlabel('Problem 5','FontSize',24, 'interpreter', 'latex');

nexttile(8);
ribbon(GMRES_nmbittoconv_RAS_Prblm2'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0); zticklabels({});

nexttile(9);
ribbon(GMRES_nmbittoconv_MS_Prblm2'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0); zticklabels({});

cbar = colorbar(); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
cbar.Layout.Tile = 'east';



% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontSize',33);
end



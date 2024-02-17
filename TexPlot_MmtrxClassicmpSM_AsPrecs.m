clear; clc; close('all')
% addpath(genpath('./RankTruncData'))
% MyData = {'Schwarz methods as preconds','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'theta for dAS', dampingTheta,'# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, ...
%     'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
%     'GMRES prec type (noprec/L/R)', GMRES_PrecType, 'max # GMRES iters', GMRES_nmbiter, 'GMRES conv. tolerance', GMRES_relresacc, ...
%     'GMRES ConvCrvs noprec', GMRES_resnorms_noprec, 'GMRES # iter no prec', GMRES_nmbittoconv_noprec, ...
%     'GMRES ConvCrvs dAS', GMRES_resnorms_dAS, 'GMRES # iter dAS', GMRES_nmbittoconv_dAS, ... 
%     'GMRES ConvCrvs RAS', GMRES_resnorms_RAS, 'GMRES # iter RAS', GMRES_nmbittoconv_RAS, ... 
%     'GMRES ConvCrvs MS', GMRES_resnorms_MS, 'GMRES # iter MS', GMRES_nmbittoconv_MS};

ProblemChoice = 3;
s_LoadString = append(append('SavedData_ClassicmpSM_Mmtrx_AsPrecs_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{8}; list_of_nmb_int_gridcols = LoadedData{14}; GMRES_nmbiter = LoadedData{18};
GMRES_resnorms_noprec = LoadedData{22}; GMRES_nmbittoconv_noprec = LoadedData{24};
GMRES_resnorms_dAS = LoadedData{26}; GMRES_nmbittoconv_dAS = LoadedData{28};
GMRES_resnorms_RAS = LoadedData{30}; GMRES_nmbittoconv_RAS = LoadedData{32};
GMRES_resnorms_MS = LoadedData{34}; GMRES_nmbittoconv_MS = LoadedData{36};































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



%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
ind_MeshSizeToPlot = length(list_of_nmb_int_gridcols); angle1 = 227; angle2 = 26;

% left subplot are convergence curves of MS-prec GMRES for different d_s
nmb_iters_ToPlot = GMRES_nmbiter; iters_mesh = 0:nmb_iters_ToPlot-1;
nmb_DigsToPlot = 6;

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_MS{ind_MeshSizeToPlot,ind_nmbdig};
    plt1 = semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
    LegendHandles(ind_nmbdig+1) = plt1; LegendLabels{ind_nmbdig+1} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
end
PlotData = GMRES_resnorms_noprec{ind_MeshSizeToPlot};
plt1 = semilogy(iters_mesh, PlotData,'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); hold on;
LegendHandles(1) = plt1; LegendLabels{1} = 'no precond.';

xlabel('iteration','FontSize',24,'interpreter', 'latex'); ylabel('rel. res.', 'FontSize',24,'interpreter', 'latex'); %title('convergence curves for MS','FontSize',24, 'interpreter', 'latex');
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles(1:end),LegendLabels{1:end},'Location','east','Orientation','Horizontal'); Lgnd.NumColumns = 1;
fontsize(Lgnd,25,'points'); 


% right subplot are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols
nexttile();
ribbon(GMRES_nmbittoconv_RAS'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});
cbar = colorbar; cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';


% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19);
end







%%% plot interaction of "d_s " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = 2; nmbColsTiles = 3;
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','loose','Padding','Compact');

ind_MeshSizeToPlot = length(list_of_nmb_int_gridcols); 


% first row: convergence curves of dAS/RAS/MS-prec GMRES for different d_s
%%%%%
nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_dAS{ind_MeshSizeToPlot,ind_nmbdig};
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData,'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('GMRES preconditioned with dAS','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_RAS{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData,'Color','k','Marker','o','MarkerSize',12,'LineWidth',2); 
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('GMRES preconditioned with RAS','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = GMRES_resnorms_MS{ind_MeshSizeToPlot,ind_nmbdig}; semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
PlotData = GMRES_resnorms_noprec{ind_MeshSizeToPlot}; semilogy(iters_mesh, PlotData,'Color','k','Marker','o','MarkerSize',12,'LineWidth',2);
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('GMRES preconditioned with MS','FontSize',24, 'interpreter', 'latex');

% legend
plt1 = plot(nan,nan,'Color','black','LineWidth',6); hold on;
LegendHandles(1) = plt1; LegendLabels{1} = 'no precond.';
for ind_nmbdig = nmb_DigsToPlot:-1:1
    plt1 = plot(nan,nan,'Color',[MyColors{ind_nmbdig}],'LineWidth',6); hold on;
    LegendHandles(1+ind_nmbdig) = plt1; LegendLabels{1+ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
end
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','d','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{1+nmb_DigsToPlot+1} = ' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$'; LegendHandles(1+nmb_DigsToPlot+1) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','^','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{1+nmb_DigsToPlot+2} = '$\; \mathrm{RAS} \quad$'; LegendHandles(1+nmb_DigsToPlot+2) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','o','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{1+nmb_DigsToPlot+3} = '$\; \mathrm{MS} \quad$'; LegendHandles(1+nmb_DigsToPlot+3) = plt1;
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 11; Lgnd.Layout.Tile = 'North';


% second row: are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols
%%%%%
nexttile();
ribbon(GMRES_nmbittoconv_dAS'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

nexttile();
ribbon(GMRES_nmbittoconv_RAS'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

nexttile();
ribbon(GMRES_nmbittoconv_MS'); view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

cbar = colorbar(); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
% cbar.Layout.Tile = 'east';



% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19);
end



% addpath(genpath('./RankTruncData'))
% MyData = {'Schwarz methods','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'SM #iter', SM_nmbiter,'SM relres accuracy stop', SM_relresacc, 'theta for dAS', dampingTheta, ...
%    '# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, 'Err Mtrx conditions calculated', CalcErrMtrx, ...
%     'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
%     'exact solution (only for the largest problem)', u_ExactSol, ...
%     'ConvCurvs dAS', ConvCurves_dAS, 'approx ConvFacts dAS', ConvFactApprox_dAS, 'Error plots dAS (only for the largest problem)', PlotData_ErrPlot_dAS, ...
%     'ConvCurvs RAS', ConvCurves_RAS, 'approx ConvFacts RAS', ConvFactApprox_RAS, 'Error plots RAS (only for the largest problem)', PlotData_ErrPlot_RAS, ...
%     'ConvCurvs MS', ConvCurves_MS, 'approx ConvFacts MS', ConvFactApprox_MS, 'Error plots MS (only for the largest problem)', PlotData_ErrPlot_MS, ...
%     'RM_nmbdigs_list indices for which we plot error', indsDigs_PlotErr, 'SM ityeration indices for which we plot error', indsIter_PlotErr, 'RM_nmbdigs_list inds that satisfy the conv cond', inds_DigsThatStsfyConvCond_dAS};

clear; clc; close('all')
ProblemChoice = 3;
s_LoadString = append(append('SavedData_ClassicmpSM_Mmtrx_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{12}; CalcErrMtrx  = LoadedData{16}; list_of_nmb_int_gridcols = LoadedData{20};
u_ExactSol = LoadedData{22};
ConvCurves_dAS = LoadedData{24}; ConvFactApprox_dAS = LoadedData{26}; PlotData_ErrPlot_dAS = LoadedData{28};
ConvCurves_RAS = LoadedData{30}; ConvFactApprox_RAS = LoadedData{32}; PlotData_ErrPlot_RAS = LoadedData{34};
ConvCurves_MS = LoadedData{36}; ConvFactApprox_MS = LoadedData{38}; PlotData_ErrPlot_MS = LoadedData{40};
indsDigs_PlotErr = LoadedData{42}; indsIter_PlotErr = LoadedData{44}; 


% %%% place holder before we re-run the experiments
% inds_DigsThatStsfyConvCond_dAS = LoadedData{46};
% inds_FrstIndDigThatStsfyConvCond = ones(length(list_of_nmb_int_gridcols),1)*inds_DigsThatStsfyConvCond_dAS;
%%% after we re-run, replace with
inds_FrstIndDigThatStsfyConvCond = LoadedData{46};









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'o','d','s','^','*','>'}; MyLines = {'--','-.',':','-'};



%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(1,2,'TileSpacing','loose','Padding','Compact');
ind_MeshSizeToPlot = 1;


% left subplot are convergence curves of MS for different d_s
nmb_iters_ToPlot = 60; iters_mesh = 0:nmb_iters_ToPlot;
nmb_DigsToPlot = 7;

ax1 = nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_MS(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    plt1 = semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
    LegendHandles(ind_nmbdig) = plt1; LegendLabels{ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
end
xlabel('iteration','FontSize',24,'interpreter', 'latex'); ylabel('error 2-norm', 'FontSize',24,'interpreter', 'latex'); %title('convergence curves for MS','FontSize',24, 'interpreter', 'latex');
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles(1:end),LegendLabels{1:end},'Location','southwest','Orientation','Horizontal'); Lgnd.NumColumns = 1;
fontsize(Lgnd,25,'points'); set(ax1,'FontWeight','bold','FontSize',19,'Color','white');


% right subplot are approximate convergence factors of MS, dAS,RAS for different d_s
LegendHandles = []; LegendLabels = {};
ax2 = nexttile();
plot(RM_nmbdigits_list,ConvFactApprox_dAS(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(1) = plt; LegendLabels{1} = 'dAS, $\theta = \frac{1}{3}$';
if CalcErrMtrx, inds_to_color = inds_FrstIndDigThatStsfyConvCond(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_dAS(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_RAS(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(2) = plt; LegendLabels{2} = 'RAS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_RAS(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_MS(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(3) = plt; LegendLabels{3} = 'MS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_MS(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on; end

xlabel('$d_s$','FontSize',24, 'interpreter', 'latex'); ylabel('approx. $\rho_{\mathrm{conv}}$','FontSize',24,'interpreter', 'latex'); %title('convergence curves for MS', 'interpreter', 'latex');

% legend
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','d','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{1} = ' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$'; LegendHandles(1) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','^','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{2} = '$\; \mathrm{RAS} \quad$'; LegendHandles(2) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','o','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{3} = '$\; \mathrm{MS} \quad$'; LegendHandles(3) = plt1;
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 3; Lgnd.Layout.Tile = 'North';
pos=Lgnd.Position;        % retrieve existing position
pos(3)=15*pos(3);       % increase width value 50% in position 4-vector
Lgnd.Position=pos;        % set new position
fontsize(Lgnd,25,'points'); set(ax2,'FontWeight','bold','FontSize',19,'Color','white');





% %%% plot the errors -- old
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nmbRowsTiles = length(indsDigs_PlotErr); nmbColsTiles = length(indsIter_PlotErr)+1;
% figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
% 
% ind_MeshSizeToPlot = 1;%length(list_of_nmb_int_gridcols); 
% nmb_int_gridcols = list_of_nmb_int_gridcols(ind_MeshSizeToPlot); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;
% 
% for ind_dig = 1:length(indsDigs_PlotErr)
% 
%     %%% first column is the initial error, i.e., the exact solution if we took initial guess as all zeros
%     nexttile()
%     PlotData = u_plot;
%     PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_ExactSol,nmb_int_gridcols,nmb_int_gridcols)'; 
%     mesh(x_mesh,x_mesh,PlotData); xlabel('$x_1$','FontSize',24,'interpreter', 'latex');ylabel('$x_2$','FontSize',24,'interpreter', 'latex'); view(angle1,angle2);
%     u_plot = zeros(nmb_int_gridcols+2);
%     if ind_dig == 1; myT = 'solution'; title(myT, 'interpreter', 'latex'); end
%     myT = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); zlabel(myT, 'interpreter', 'latex');
% 
% 
%     for ind_iter = 1:length(indsIter_PlotErr)
% 
%         nexttile()
%         u_err = reshape( PlotData_ErrPlot_MS{indsDigs_PlotErr(ind_dig),ind_iter}, 1,nmb_int_gridcols^2);
%         PlotData = u_plot;
%         PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
%         mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); xlabel('$x_1$','FontSize',24,'interpreter', 'latex');ylabel('$x_2$','FontSize',24,'interpreter', 'latex');
%         if ind_dig == 1; myT = append( append('error after ',num2str(indsIter_PlotErr(ind_iter))), ' it.'); title(myT, 'interpreter', 'latex'); end
%         u_plot = zeros(nmb_int_gridcols+2);
% 
%     end
% end
% 
% % font sizes
% for ind=1:length(t.Children)
%     curr_ax = t.Children(ind);
%     set(curr_ax,'FontWeight','bold','FontSize',19,'Color','white');
% end

%%% plot the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = length(indsIter_PlotErr); nmbColsTiles = length(indsDigs_PlotErr);
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

ind_MeshSizeToPlot = 1;%length(list_of_nmb_int_gridcols); 
nmb_int_gridcols = list_of_nmb_int_gridcols(ind_MeshSizeToPlot); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;

for ind_iter = 1:length(indsIter_PlotErr)

    % %%% first column is the initial error, i.e., the exact solution if we took initial guess as all zeros
    % nexttile()
    % if ind_dig == 1; myT = 'solution'; title(myT, 'interpreter', 'latex'); end
    % myT = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); zlabel(myT, 'interpreter', 'latex');


    
    for ind_dig = 1:length(indsDigs_PlotErr)
        nexttile()
        u_err = reshape( PlotData_ErrPlot_MS{indsDigs_PlotErr(ind_dig),ind_iter}, 1,nmb_int_gridcols^2);
        PlotData = u_plot;
        PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
        mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); xlabel('$x_1$','FontSize',24,'interpreter', 'latex');ylabel('$x_2$','FontSize',24,'interpreter', 'latex');
        if ind_dig == 1 
            mylabel = append( append('error after ',num2str(indsIter_PlotErr(ind_iter))), ' it.'); zlabel(mylabel, 'interpreter', 'latex'); 
            if ind_iter == 3, zl=get(gca,'zlabel'); pzl = get(zl,'position'); pzl(1) = 1.001*pzl(1); set(zl,'position',pzl); end
        end
        if ind_iter == 1, mylabel = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); title(mylabel, 'interpreter', 'latex'); end
        u_plot = zeros(nmb_int_gridcols+2);

    end
end

% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19,'Color','white');
end






%%% plot interaction of "d_s " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = 2; nmbColsTiles = 3; 
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Loose','Padding','Compact');

ind_MeshSizeToPlot = 1; angle1 = 227; angle2 = 26;
Cbar_Ticks = []; Cbar_TickLabels = {};
for ind = 1:length(list_of_nmb_int_gridcols)
    if mod(ind-1,3) == 0
    Cbar_Ticks = [Cbar_Ticks,ind];
    Cbar_TickLabels = {Cbar_TickLabels{:}, int2str( list_of_nmb_int_gridcols(ind)^2 ) };
    end
end


% first row: convergence curves of dAS/RAS/MS-prec GMRES for different d_s
%%%%%
nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_dAS(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','d','MarkerSize',12,'LineWidth',2); hold on;
end
xlabel('iteration','FontSize',24,'interpreter', 'latex'); ylabel('error 2-norm', 'FontSize',24,'interpreter', 'latex'); title('damped additive Schwarz','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_RAS(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','^','MarkerSize',12,'LineWidth',2); hold on;
end
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('restricted additive Schwarz','FontSize',24, 'interpreter', 'latex');

nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_MS(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
xlabel('iteration','FontSize',24,'interpreter', 'latex'); title('multiplicative Schwarz','FontSize',24, 'interpreter', 'latex');

% legend for the first row plots
for ind_nmbdig = nmb_DigsToPlot:-1:1
    plt1 = semilogy(nan,nan,'Color',[MyColors{ind_nmbdig}],'Marker','none','MarkerSize',12,'LineWidth',2); hold on;
    LegendHandles(ind_nmbdig) = plt1; LegendLabels{ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
end
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 7; Lgnd.Layout.Tile = 'North';
fontsize(Lgnd,25,'points'); 



% second row: are # prec GMRES iter to converge for diff d_s and nmb_int_gridcols
%%%%%
nexttile();
ribbon(ConvFactApprox_dAS'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond(ind); z = ConvFactApprox_dAS(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});
zlabel('approx. $\rho_{\mathrm{conv}}$','FontSize',24,'interpreter', 'latex');

nexttile();
ribbon(ConvFactApprox_RAS'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond(ind); z = ConvFactApprox_RAS(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});

nexttile(); LegendHandles = []; LegendLabels = {};
ribbon(ConvFactApprox_MS');  hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond(ind); z = ConvFactApprox_MS(x,y);
    plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
    LegendHandles(1) = plt; LegendLabels{1} = 'first $d_s$ such that $\| A_i^{-1} \mathtt{E}_i \| < 1$';
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});
cbar = colorbar; cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';

set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','southoutside','Orientation','Horizontal'); Lgnd.NumColumns = 7; Lgnd.Layout.Tile = 'South';
fontsize(Lgnd,25,'points'); 
% ah1=axes('position',get(gca,'position'),'visible','off');
% leg2=legend(ah1,p1,magnitude{ind},'Location','NorthWest');set(leg2,'FontSize',9);


% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19);%,'Color','white');
end


% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19);%,'Color','white');
end


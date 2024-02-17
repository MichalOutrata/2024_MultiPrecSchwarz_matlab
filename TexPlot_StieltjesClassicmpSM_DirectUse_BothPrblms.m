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


ProblemChoice = 1;
s_LoadString = append(append('SavedData_ClassicmpSM_Stieltjes_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{12}; CalcErrMtrx  = LoadedData{16}; list_of_nmb_int_gridcols = LoadedData{20};
u_ExactSol_Prblm1 = LoadedData{22};
ConvCurves_dAS_Prblm1 = LoadedData{24}; ConvFactApprox_dAS_Prblm1 = LoadedData{26}; PlotData_ErrPlot_dAS_Prblm1 = LoadedData{28};
ConvCurves_RAS_Prblm1 = LoadedData{30}; ConvFactApprox_RAS_Prblm1 = LoadedData{32}; PlotData_ErrPlot_RAS_Prblm1 = LoadedData{34};
ConvCurves_MS_Prblm1 = LoadedData{36}; ConvFactApprox_MS_Prblm1 = LoadedData{38}; PlotData_ErrPlot_MS_Prblm1 = LoadedData{40};
indsDigs_PlotErr_Prblm1 = LoadedData{42}; indsIter_PlotErr_Prblm1 = LoadedData{44}; 

% %%% place holder before we re-run the experiments
% inds_DigsThatStsfyConvCond_dAS_Prblm3 = LoadedData{46};
% inds_FrstIndDigThatStsfyConvCond_Prblm3 = ones(length(list_of_nmb_int_gridcols),1)*inds_DigsThatStsfyConvCond_dAS_Prblm3;
%%% after we re-run, replace with
inds_FrstIndDigThatStsfyConvCond_Prblm1 = LoadedData{46};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ProblemChoice = 2;
s_LoadString = append(append('SavedData_ClassicmpSM_Stieltjes_Prblm',num2str(ProblemChoice)),'.mat'); LoadedData = load(s_LoadString).MyData;
RM_nmbdigits_list = LoadedData{12}; CalcErrMtrx  = LoadedData{16}; list_of_nmb_int_gridcols = LoadedData{20};
u_ExactSol_Prblm2 = LoadedData{22};
ConvCurves_dAS_Prblm2 = LoadedData{24}; ConvFactApprox_dAS_Prblm2 = LoadedData{26}; PlotData_ErrPlot_dAS_Prblm2 = LoadedData{28};
ConvCurves_RAS_Prblm2 = LoadedData{30}; ConvFactApprox_RAS_Prblm2 = LoadedData{32}; PlotData_ErrPlot_RAS_Prblm2 = LoadedData{34};
ConvCurves_MS_Prblm2 = LoadedData{36}; ConvFactApprox_MS_Prblm2 = LoadedData{38}; PlotData_ErrPlot_MS_Prblm2 = LoadedData{40};
indsDigs_PlotErr_Prblm2 = LoadedData{42}; indsIter_PlotErr_Prblm2 = LoadedData{44}; 

% %%% place holder before we re-run the experiments
% inds_DigsThatStsfyConvCond_dAS_Prblm4 = LoadedData{46};
% inds_FrstIndDigThatStsfyConvCond_Prblm4 = ones(length(list_of_nmb_int_gridcols),1)*inds_DigsThatStsfyConvCond_dAS_Prblm3;
%%% after we re-run, replace with
inds_FrstIndDigThatStsfyConvCond_Prblm2 = LoadedData{46};












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Clean data     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
nmb_iters_ToPlot = 60; iters_mesh = 0:nmb_iters_ToPlot; Max_ds = 4;

%%% Prblm 1, dAS -> mesh sizes 1,2,3,4,5,10 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; t(1) = tiledlayout(2,3,'TileSpacing','compact','Padding','Compact');
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(1,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(2,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(3,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% 
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(4,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(5,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm1(10,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
%%% hence:
ConvFactApprox_dAS_Prblm1(2:end,1:2) = nan(27,2);


%%% Prblm 1, RAS -> mesh sizes 10 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; t(1) = tiledlayout(1,1,'TileSpacing','compact','Padding','Compact');
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_RAS_Prblm1(10,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
%%% hence:
ConvFactApprox_RAS_Prblm1(:,1:2) = nan(28,2); 
ConvCrv_mshsiz10 = reshape(ConvCurves_RAS_Prblm1(10,3,:), 1,length(ConvCurves_RAS_Prblm1(10,3,:))); ConvCrv_mshsiz10_cut = ConvCrv_mshsiz10(1:35);
ConsecErrsRatio = ConvCrv_mshsiz10_cut(2:end)./ConvCrv_mshsiz10_cut(1:end-1);
ConvFactApproxSequence = ( ConsecErrsRatio(1:end-2) + ConsecErrsRatio(2:end-1) ) / 2;
% figure(); plot(ConvFactApproxSequence);
ConvFactApprox_RAS_Prblm1(10,3) = ConvFactApproxSequence(18);

%%% Prblm 1, MS -> mesh sizes 1,10 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(1,3,'TileSpacing','compact','Padding','Compact');
nexttile()
for ind_nmbdig = Max_ds:-1:1
    PlotData = reshape(ConvCurves_MS_Prblm1(1,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
legend();
nexttile()
for ind_nmbdig = Max_ds:-1:1
    PlotData = reshape(ConvCurves_MS_Prblm1(10,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
%%% hence:
ConvFactApprox_MS_Prblm1(:,1:2) = nan(28,2);
ConvCrv_mshsiz10 = reshape(ConvCurves_MS_Prblm1(10,3,:), 1,length(ConvCurves_MS_Prblm1(10,3,:))); ConvCrv_mshsiz10_cut = ConvCrv_mshsiz10(1:16);
ConsecErrsRatio = ConvCrv_mshsiz10_cut(2:end)./ConvCrv_mshsiz10_cut(1:end-1);
ConvFactApproxSequence = ( ConsecErrsRatio(1:end-2) + ConsecErrsRatio(2:end-1) ) / 2;
%figure(); plot(ConvFactApproxSequence);
ConvFactApprox_MS_Prblm1(10,3) = ConvFactApproxSequence(10);



%%% Prblm 2, dAS -> mesh sizes 1,2,3,10 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; t(1) = tiledlayout(1,4,'TileSpacing','compact','Padding','Compact');
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm2(1,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm2(2,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm2(3,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_dAS_Prblm2(10,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
%%% hence:
ConvFactApprox_dAS_Prblm2(2:end,2) = nan(27,1); ConvFactApprox_dAS_Prblm2(4:end,1) = nan(25,1);
ConvCrv_mshsiz10 = reshape(ConvCurves_dAS_Prblm2(10,3,:), 1,length(ConvCurves_dAS_Prblm2(10,3,:))); ConvCrv_mshsiz10_cut = ConvCrv_mshsiz10(1:23);
ConsecErrsRatio = ConvCrv_mshsiz10_cut(2:end)./ConvCrv_mshsiz10_cut(1:end-1);
ConvFactApproxSequence = ( ConsecErrsRatio(1:end-2) + ConsecErrsRatio(2:end-1) ) / 2;
%figure(); plot(ConvFactApproxSequence);
ConvFactApprox_dAS_Prblm2(10,3) = ConvFactApproxSequence(18);



%%% Prblm 2, RAS -> mesh sizes 1,2,10 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; t(1) = tiledlayout(1,4,'TileSpacing','compact','Padding','Compact');
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_RAS_Prblm2(1,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_RAS_Prblm2(2,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_RAS_Prblm2(10,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
%%% hence:
ConvFactApprox_RAS_Prblm2(2:end,1:2) = nan(27,2);


%%% Prblm 1, MS -> mesh sizes 1,2 (checked manually the conv. curves) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; t(1) = tiledlayout(1,2,'TileSpacing','compact','Padding','Compact');
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_MS_Prblm2(1,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
% legend();
% nexttile()
% for ind_nmbdig = Max_ds:-1:1
%     PlotData = reshape(ConvCurves_MS_Prblm2(2,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
%     semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
% end
%%% hence:
ConvFactApprox_MS_Prblm2(2:end,1:2) = nan(27,2);

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'o','d','s','^','*','>'}; MyLines = {'--','-.',':','-'};



%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(2,2,'TileSpacing','compact','Padding','Compact');
ind_MeshSizeToPlot = 1;
nmb_iters_ToPlot = 60; iters_mesh = 0:nmb_iters_ToPlot;
nmb_DigsToPlot = 6;

% top row: convergence curves of MS for different d_s for Problem 1
nexttile(1); 
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_MS_Prblm1(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
title('Problem 4', 'FontSize',24,'FontWeight','bold','interpreter', 'latex'); xlabel('iteration','FontSize',24,'interpreter', 'latex'); %ylabel('2-norm of the error','FontSize',24, 'interpreter', 'latex');
nexttile(2); 
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves_MS_Prblm2(ind_MeshSizeToPlot,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
    semilogy(iters_mesh, PlotData,'Color',[MyColors{ind_nmbdig}],'Marker','o','MarkerSize',12,'LineWidth',2); hold on;
end
ylim([1e-13,1])
title('Problem 5', 'FontSize',24,'FontWeight','bold','interpreter', 'latex'); xlabel('iteration','FontSize',24,'interpreter', 'latex');


% bottom row: observed convergence factors of MS,dAS,RAS for different d_s for Problem 1
nexttile(3);
plot(RM_nmbdigits_list,ConvFactApprox_dAS_Prblm1(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(1) = plt; LegendLabels{1} = 'dAS, $\theta = \frac{1}{3}$';
if CalcErrMtrx, inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_dAS_Prblm1(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_RAS_Prblm1(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(2) = plt; LegendLabels{2} = 'RAS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_RAS_Prblm1(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_MS_Prblm1(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(3) = plt; LegendLabels{3} = 'MS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_MS_Prblm1(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on; end
xlabel('$d_s$','FontSize',24, 'interpreter', 'latex'); %title('observed $\rho_{\mathrm{conv}}$','FontSize',24,'interpreter', 'latex');

nexttile(4);
plot(RM_nmbdigits_list,ConvFactApprox_dAS_Prblm2(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(1) = plt; LegendLabels{1} = 'dAS, $\theta = \frac{1}{3}$';
if CalcErrMtrx, inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_dAS_Prblm2(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_RAS_Prblm2(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(2) = plt; LegendLabels{2} = 'RAS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_RAS_Prblm2(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on; end

plot(RM_nmbdigits_list,ConvFactApprox_MS_Prblm2(ind_MeshSizeToPlot,:),'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(3) = plt; LegendLabels{3} = 'MS';
if CalcErrMtrx,  inds_to_color = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind_MeshSizeToPlot):length(RM_nmbdigits_list);
    plot(inds_to_color,ConvFactApprox_MS_Prblm2(ind_MeshSizeToPlot,inds_to_color),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on; end
xlabel('$d_s$','FontSize',24, 'interpreter', 'latex');



% legend
for ind_nmbdig = nmb_DigsToPlot:-1:1
    plt1 = plot(nan,nan,'Color',[MyColors{ind_nmbdig}],'LineWidth',6); hold on;
    LegendHandles(ind_nmbdig) = plt1; LegendLabels{ind_nmbdig} = append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig)));
    if ind_nmbdig == nmb_DigsToPlot, LegendLabels{ind_nmbdig} = append(append('$d_s$ = ',num2str(RM_nmbdigits_list(ind_nmbdig))),' $\quad$ '); end
end
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','d','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{nmb_DigsToPlot+1} = ' $\; \mathrm{dAS \; with} \; \theta = \frac{1}{3} \quad$'; LegendHandles(nmb_DigsToPlot+1) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','^','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{nmb_DigsToPlot+2} = '$\; \mathrm{RAS} \quad$'; LegendHandles(nmb_DigsToPlot+2) = plt1;
plt1 = plot(nan, nan, 'color', '[.7 .7 .7]','LineStyle', '-','Marker','o','MarkerSize',16,'LineWidth',2); hold on;
LegendLabels{nmb_DigsToPlot+3} = '$\; \mathrm{MS} \quad$'; LegendHandles(nmb_DigsToPlot+3) = plt1;
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 10; Lgnd.Layout.Tile = 'North';

% font sizes
fontsize(Lgnd,35,'points');
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontSize',30);%,'Color','white');
end










%%% plot the errors - Problem 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = length(indsIter_PlotErr_Prblm1); nmbColsTiles = length(indsDigs_PlotErr_Prblm1)-1;
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

ind_MeshSizeToPlot = 1;%length(list_of_nmb_int_gridcols); 
nmb_int_gridcols = list_of_nmb_int_gridcols(ind_MeshSizeToPlot); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;

for ind_iter = 1:length(indsIter_PlotErr_Prblm1)

    % %%% first column is the initial error, i.e., the exact solution if we took initial guess as all zeros
    % nexttile()
    % if ind_dig == 1; myT = 'solution'; title(myT, 'interpreter', 'latex'); end
    % myT = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); zlabel(myT, 'interpreter', 'latex');



    for ind_dig = 1:length(indsDigs_PlotErr_Prblm1)-1
        nexttile()
        u_err = reshape( PlotData_ErrPlot_MS_Prblm1{indsDigs_PlotErr_Prblm1(ind_dig),ind_iter}, 1,nmb_int_gridcols^2);
        PlotData = u_plot;
        PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
        mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); 
        if ind_iter == length(indsIter_PlotErr_Prblm1), xlabel('$x_1$','interpreter', 'latex');ylabel('$x_2$','FontSize',26,'interpreter', 'latex'); end
        max_err = max(max(abs(u_err))); max_err_exp = log10(max_err);
        ax = gca; ax.ZAxis.Exponent = floor(max_err_exp);
        % if ind_dig == 1 
        %     mylabel = append( append('error after ',num2str(indsIter_PlotErr_Prblm1(ind_iter))), ' it.'); zlabel(mylabel, 'interpreter', 'latex'); 
        %     if ind_iter == 3, zl=get(gca,'zlabel'); pzl = get(zl,'position'); pzl(1) = 1.001*pzl(1); set(zl,'position',pzl); end
        % end
        if ind_iter == 1, mylabel = append('$d_s$ = ',num2str(indsDigs_PlotErr_Prblm1(ind_dig))); title(mylabel, 'interpreter', 'latex'); end
        u_plot = zeros(nmb_int_gridcols+2);

    end
end

% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontSize',35);
end















%%% plot interaction of "d_s " and "h"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = 2; nmbColsTiles = 3; 
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

ind_MeshSizeToPlot = 1; angle1 = 227; angle2 = 26;
Cbar_Ticks = []; Cbar_TickLabels = {};
for ind = 1:length(list_of_nmb_int_gridcols)
    if mod(ind-1,3) == 0
    Cbar_Ticks = [Cbar_Ticks,ind];
    Cbar_TickLabels = {Cbar_TickLabels{:}, int2str( list_of_nmb_int_gridcols(ind)^2 ) };
    end
end

% first row: observed convergence factor for diff d_s and nmb_int_gridcols for Problem 3
%%%%% 
nexttile();
ribbon(ConvFactApprox_dAS_Prblm1'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind); z = ConvFactApprox_dAS_Prblm1(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2); title('damped additive Schwarz', 'interpreter', 'latex');
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({});
zlabel('Problem 4','FontSize',24,'interpreter', 'latex');

nexttile();
ribbon(ConvFactApprox_RAS_Prblm1'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind); z = ConvFactApprox_RAS_Prblm1(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2); title('restricted additive Schwarz', 'interpreter', 'latex');
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({});

nexttile(); LegendHandles = []; LegendLabels = {};
ribbon(ConvFactApprox_MS_Prblm1');  hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm1(ind); z = ConvFactApprox_MS_Prblm1(x,y);
    plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
    LegendHandles(1) = plt; LegendLabels{1} = 'first $d_s$ such that $\| A_i^{-1} \mathtt{E}_i \| < 1$';
end
view(angle1,angle2);  title('multiplicative Schwarz', 'interpreter', 'latex');
xticks([]); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({});


% first row: observed convergence factor for diff d_s and nmb_int_gridcols for Problem 4
%%%%%
nexttile();
ribbon(ConvFactApprox_dAS_Prblm2'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind); z = ConvFactApprox_dAS_Prblm2(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2); 
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0);
zlabel('Problem 5','FontSize',24,'interpreter', 'latex');

nexttile();
ribbon(ConvFactApprox_RAS_Prblm2'); hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind); z = ConvFactApprox_RAS_Prblm2(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0);

nexttile(); LegendHandles = []; LegendLabels = {};
ribbon(ConvFactApprox_MS_Prblm2');  hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond_Prblm2(ind); z = ConvFactApprox_MS_Prblm2(x,y);
    plt = plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
    LegendHandles(1) = plt; LegendLabels{1} = 'first $d_s$ such that $\| A_i^{-1} \mathtt{E}_i \| < 1$';
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'}); ytickangle(0);
cbar = colorbar('Location','eastoutside'); cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';
cbar.Layout.Tile = 'east';

% Legend
set(groot,'defaultLegendInterpreter','latex');
Lgnd = legend(LegendHandles,LegendLabels,'Location','northoutside','Orientation','Horizontal'); Lgnd.NumColumns = 7; Lgnd.Layout.Tile = 'North';
fontsize(Lgnd,35,'points'); 


% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontSize',33);%,'Color','white');
end


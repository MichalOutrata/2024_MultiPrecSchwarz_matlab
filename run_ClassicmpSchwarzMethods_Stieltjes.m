clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
SM_nmbiter = 61; % number of RAS iterations
SM_relresacc = 1e-12; % number of RAS iterations
RM_nmbdigits_list = 1:16; % number of digits to keep for the subdomain solves
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
Advanpix = true; CalcErrMtrx = true; DoPlotting = false; SavePlotData = true;

%%% plotting
indsIter_PlotErr = [10,20,60];
indsDigs_PlotErr = [2,3,4,5];
indsMesh_PlotErr = 1;
inds_FrstIndDigThatStsfyConvCond = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 50:10:320; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 22;



for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
    nmb_int_gridcols = list_of_nmb_int_gridcols(ind_meshsize); h=1/(nmb_int_gridcols+1); % mesh size
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(append('Problem size = ',num2str(nmb_int_gridcols^2)))

    if ProblemChoice == 0 %%% do negative laplacian
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) zeros(size(x)); a = @(x,y) ones(size(x)); b1 = @(x,y) zeros(size(x)); b2 = b1;
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1); rhs = 1/h^2*rhs; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
        u_init = zeros(size(A,1),1); u_plot = zeros(nmb_int_gridcols+2); %u_plot(2:nmb_int_gridcols+1,1) = 1; 
        angle1 = 52.5; angle2 = 30; x_mesh = 0:h:1;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';
    
    elseif ProblemChoice == 1 %%% do symmetric AdvecDiff based on SzyldGander eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; alpha =  @(x,y) (x+y).^2.*exp(x-y); nu = @(x,y) 0*(1.*x+1.*y); mu = @(x,y) 0*(1.*x+1.*y);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,nu,mu});
        rhs = ones(size(A,1),1)*1/sqrt(nmb_int_gridcols^2)*1/h^2; u_init = randn(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';
    
    elseif ProblemChoice == 2 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2); 
        eta = @(x,y) 500.*x+1.*y; alpha = @(x,y) 1+9.*(x+y); beta = @(x,y) 1+9.*(x+y); nu = @(x,y) 0.*x+0.*y; mu = @(x,y) 0.*x+0.*y;
        A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
        rhs = ones(size(A,1),1)*1/sqrt(nmb_int_gridcols^2)*1/h^2; u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';

    elseif ProblemChoice == 22 %%% based on Problem 4, only we added a region of much higher diffusion to artificially bump up condition number
        G=numgrid('S',nmb_int_gridcols+2);
        AdvectStrngth = 0;
        addpath(genpath('./CoeffFuncs_ReacAdvecDiff'));
        eta_handle = @(x,y) eta(x,y); a_handle = @(x,y) alpha(x,y); b1_handle = @(x,y) b1(x,y,AdvectStrngth); b2_handle = @(x,y) b2(x,y,AdvectStrngth);
        A=ReactAdvDiff_Sqr_FD('SzyldGander_SclrFeval',G,{eta_handle,a_handle,b1_handle,b2_handle});
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';

    % %%% optional more advanced rhs building
    % [m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h); IntriorInds = find(G);
    % rhsfun = @(x,y) x.*y.*(1-x).*(1-y); rhsfun_discr = feval(rhsfun,X,Y); rhs = A*rhsfun_discr(IntriorInds);
    end

    if ind_meshsize == indsMesh_PlotErr, u_ExactSol = A\rhs; end

    debug = 1; if DoPlotting, MyPlot = {x_mesh,u_plot,0.05,angle1,angle2}; else, MyPlot = {}; end % MyPlot = {x,u_plot,waitime,angle1,angle2};
    u_PlotErr = zeros(length(indsDigs_PlotErr),length(indsIter_PlotErr),length(rhs));


    % run the dAS
    SchwarzMethod = {'dAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   dAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);

        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves_dAS(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        if ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) > 1 
            ConvFactApprox_dAS(ind_meshsize,ind_nmbdig) = nan; ConvCurves_dAS(ind_meshsize,ind_nmbdig,:) = nan(size(error_nrms));
        else
            ConvFactApprox_dAS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);
        end

        % get the errors when we want them
        if ind_meshsize == indsMesh_PlotErr
            indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
            if indDig_CurrOrdrInIndsToPlot ~= 0 
                for ind = 1:length(indsIter_PlotErr)
                    u_solseq = SMoutput{1}; 
                    u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                    u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
                    % save the data if we wanna
                    if SavePlotData, PlotData_ErrPlot_dAS{ind_nmbdig, ind} = u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:);  end
                end
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            if ~(any(SMoutput{4})), inds_FrstIndDigThatStsfyConvCond(ind_meshsize) = ind_nmbdig; CalcErrMtrx_rplc = false; end
            InvErrMtrx_normest(ind_meshsize,ind_nmbdig,:) = SMoutput{5}; 
            if ind_nmbdig == 1, CondEst_Ai(ind_meshsize,:) = SMoutput{6}; end
        end
    end


    % run the RAS
    SchwarzMethod = {'RAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; CalcErrMtrx_rplc = false; %if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   RAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);


        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves_RAS(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        if ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) > 1 
            ConvFactApprox_RAS(ind_meshsize,ind_nmbdig) = nan; ConvCurves_RAS(ind_meshsize,ind_nmbdig,:) = nan(size(error_nrms));
        else
            ConvFactApprox_RAS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);
        end

        % get the errors when we want them
        if ind_meshsize == indsMesh_PlotErr
            indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
            if indDig_CurrOrdrInIndsToPlot ~= 0 
                for ind = 1:length(indsIter_PlotErr)
                    u_solseq = SMoutput{1}; 
                    u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                    u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
                    % save the data if we wanna
                    if SavePlotData, PlotData_ErrPlot_RAS{ind_nmbdig, ind} = u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:);  end
                end
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            % if ~(any(SMoutput{4})), inds_FrstIndDigThatStsfyConvCond(ind_meshsize) = ind_nmbdig; CalcErrMtrx_rplc = false; end
        end
    end


    % run the MS
    SchwarzMethod = {'MS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; CalcErrMtrx_rplc = false; %if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   MS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);

        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves_MS(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        if ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) > 1 
            ConvFactApprox_MS(ind_meshsize,ind_nmbdig) = nan; ConvCurves_MS(ind_meshsize,ind_nmbdig,:) = nan(size(error_nrms));
        else
            ConvFactApprox_MS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);
        end

        % get the errors when we want them
        if ind_meshsize == indsMesh_PlotErr
            indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
            if indDig_CurrOrdrInIndsToPlot ~= 0 
                for ind = 1:length(indsIter_PlotErr)
                    u_solseq = SMoutput{1}; 
                    u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                    u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
                    % save the data if we wanna
                    if SavePlotData, PlotData_ErrPlot_MS{ind_nmbdig, ind} = u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:);
                    else, PlotData_ErrPlot_MS{ind_nmbdig, ind} = u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:); disp('hi, we plotting MS'); end
                end
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            % if ~(any(SMoutput{4})), inds_FrstIndDigThatStsfyConvCond(ind_meshsize) = ind_nmbdig; CalcErrMtrx_rplc = false; end
        end
    end

end





if SavePlotData
    MyData = {'Schwarz methods','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'SM #iter', SM_nmbiter,'SM relres accuracy stop', SM_relresacc, 'theta for dAS', dampingTheta, ...
       '# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, 'Err Mtrx conditions calculated', CalcErrMtrx, ...
        'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
        'exact solution (only for the largest problem)', u_ExactSol, ...
        'ConvCurvs dAS', ConvCurves_dAS, 'approx ConvFacts dAS', ConvFactApprox_dAS, 'Error plots dAS (only for the largest problem)', PlotData_ErrPlot_dAS, ...
        'ConvCurvs RAS', ConvCurves_RAS, 'approx ConvFacts RAS', ConvFactApprox_RAS, 'Error plots RAS (only for the largest problem)', PlotData_ErrPlot_RAS, ...
        'ConvCurvs MS', ConvCurves_MS, 'approx ConvFacts MS', ConvFactApprox_MS, 'Error plots MS (only for the largest problem)', PlotData_ErrPlot_MS, ...
        'RM_nmbdigs_list indices for which we plot error', indsDigs_PlotErr, 'SM ityeration indices for which we plot error', indsIter_PlotErr, 'RM_nmbdigs_list inds that satisfy the conv cond', inds_FrstIndDigThatStsfyConvCond, ...
        'The saved norms of inv_Ai_ErrMtrx', InvErrMtrx_normest, 'The saved condition number estimates for the subdomain problems', CondEst_Ai };
    s_SaveString = append(append('SavedData_ClassicmpSM_Stieltjes_Prblm',num2str(ProblemChoice)),'.mat');
    save(s_SaveString,'MyData');  
else


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'o','d','s','^','*','>'}; MyLines = {'--','-.',':','-'};



%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(1,2,'TileSpacing','loose','Padding','Compact');
ind_MeshSizeToPlot = length(list_of_nmb_int_gridcols);


% left subplot are convergence curves of MS for different d_s
nmb_iters_ToPlot = SM_nmbiter-1; iters_mesh = 0:nmb_iters_ToPlot;
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




%%% plot the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmbRowsTiles = length(indsDigs_PlotErr); nmbColsTiles = length(indsIter_PlotErr)+1;
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

ind_MeshSizeToPlot = 1;%length(list_of_nmb_int_gridcols); 
nmb_int_gridcols = list_of_nmb_int_gridcols(ind_MeshSizeToPlot); h = 1/(nmb_int_gridcols+1); u_plot = zeros(nmb_int_gridcols+2); x_mesh = 0:h:1; angle1 = 224; angle2 = 35;

for ind_dig = 1:length(indsDigs_PlotErr)

    %%% first column is the initial error, i.e., the exact solution if we took initial guess as all zeros
    nexttile()
    PlotData = u_plot;
    PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_ExactSol,nmb_int_gridcols,nmb_int_gridcols)'; 
    mesh(x_mesh,x_mesh,PlotData); xlabel('$x_1$','FontSize',24,'interpreter', 'latex');ylabel('$x_2$','FontSize',24,'interpreter', 'latex'); view(angle1,angle2);
    u_plot = zeros(nmb_int_gridcols+2);
    if ind_dig == 1; myT = 'solution'; title(myT, 'interpreter', 'latex'); end
    myT = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); zlabel(myT, 'interpreter', 'latex');


    for ind_iter = 1:length(indsIter_PlotErr)

        nexttile()
        u_err = reshape( PlotData_ErrPlot_MS{indsDigs_PlotErr(ind_dig),ind_iter}, 1,nmb_int_gridcols^2);
        PlotData = u_plot;
        PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
        mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); xlabel('$x_1$','FontSize',24,'interpreter', 'latex');ylabel('$x_2$','FontSize',24,'interpreter', 'latex');
        if ind_dig == 1; myT = append( append('error after ',num2str(indsIter_PlotErr(ind_iter))), ' it.'); title(myT, 'interpreter', 'latex'); end
        u_plot = zeros(nmb_int_gridcols+2);

    end
end

% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19,'Color','white');
end








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
        if ind_dig == 1; mylabel = append( append('error after ',num2str(indsIter_PlotErr(ind_iter))), ' it.'); zlabel(mylabel, 'interpreter', 'latex'); end
        if ind_iter == 1; mylabel = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); title(mylabel, 'interpreter', 'latex'); end
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

nexttile();
ribbon(ConvFactApprox_MS');  hold on; 
for ind = 1:length(list_of_nmb_int_gridcols)
    x = ind; y = inds_FrstIndDigThatStsfyConvCond(ind); z = ConvFactApprox_MS(x,y);
    plot3(x,y,z,'Marker','o','Color','black','MarkerSize',15,'MarkerFaceColor','[.7 .7 .7]'); hold on;
end
view(angle1,angle2);
xticks([]); ylabel('$d_s$','FontSize',24,'interpreter', 'latex'); ylim([0,17]); yticks([1,4,8,12,16]); yticklabels({'1','4','8','12','16'});
cbar = colorbar; cbar.Label.String= '$N$'; cbar.Label.Color= 'black'; 
cbar.Ticks = Cbar_Ticks; cbar.TickLabels = Cbar_TickLabels; 
cbar.Label.Interpreter = 'latex'; cbar.TickLabelInterpreter = 'latex';



% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19);%,'Color','white');
end

end

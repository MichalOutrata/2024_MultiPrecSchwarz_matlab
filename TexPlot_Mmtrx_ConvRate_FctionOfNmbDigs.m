clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_type = 'MS'; % type of Schwarz method
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
SM_nmbiter = 101; % number of RAS iterations
SM_relresacc = 1e-12; % number of RAS iterations
RM_nmbdigits_list = 1:14; % number of digits to keep for the subdomain solves
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
Advanpix = true; CalcErrMtrx = true; DoPlotting = false; SaveAllData = true;

%%% plotting
indsIter_PlotErr = [10,20,50];
indsDigs_PlotErr = [2,4,6];
inds_DigsThatStsfyConvCond_dAS = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 60; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 3;



for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
    nmb_int_gridcols = list_of_nmb_int_gridcols(ind_meshsize); h=1/(nmb_int_gridcols+1); % mesh size
    
    if ProblemChoice == 3 %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; a =  @(x,y) (x+y).^2.*exp(x-y); b1 =  @(x,y) (y-0.5); b2 =  @(x,y) -(x-0.5);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1); u_ExactSol = A\rhs;
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; RM_type = 'Mmtrx';
    
    elseif ProblemChoice == 4 %%% do non-symmetric AdvecDiff based SzyldFrommer eqn (15) and p.660 -> changed the advection strength to ~40 so that we have an Mmtrx
        G=numgrid('S',nmb_int_gridcols+2); AdvectStrngth = 100;
        eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1 + 0.*x+0.*y; beta = @(x,y) 1 + 0.*x+0.*y; nu = @(x,y) AdvectStrngth.*(1.*x.*(x-1).*(1-2.*y)); mu = @(x,y) AdvectStrngth.*(-1.*y.*(y-1).*(1-2.*x));
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,beta,nu,mu}); %A_szyld = A*h^2+speye(size(A,1)); disp(isMmtrx(A)); disp(isMmtrx(A_szyld)); disp(issymmetric(A))
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1); u_ExactSol = A\rhs;
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; x_mesh = 0:h:1; RM_type = 'Mmtrx';
    
    % %%% optional more advanced rhs building
    % [m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h); IntriorInds = find(G);
    % rhsfun = @(x,y) x.*y.*(1-x).*(1-y); rhsfun_discr = feval(rhsfun,X,Y); rhs = A*rhsfun_discr(IntriorInds);
    
    end
    
    
    debug = 1; if DoPlotting, MyPlot = {x_mesh,u_plot,0.05,angle1,angle2}; else, MyPlot = {}; end % MyPlot = {x,u_plot,waitime,angle1,angle2};
    u_PlotErr = zeros(length(indsDigs_PlotErr),length(indsIter_PlotErr),length(rhs));


    % run the dAS
    SchwarzMethod = {'dAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);
    
        % save the data if we wanna
        if SaveAllData, AllData_dAS{ind_meshsize,ind_nmbdig} = SMoutput; end
        
        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        ConvFactApprox_dAS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);

        % get the errors when we want them
        indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
        if indDig_CurrOrdrInIndsToPlot ~= 0 
            for ind = 1:length(indsIter_PlotErr)
                u_solseq = SMoutput{1}; 
                u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            if ~SMoutput{4}, inds_DigsThatStsfyConvCond_dAS = ind_nmbdig:length(RM_nmbdigits_list); CalcErrMtrx_rplc = false; end
        end


    end

    % run the RAS
    SchwarzMethod = {'RAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; CalcErrMtrx_rplc = false; %if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);
    
        % save the data if we wanna
        if SaveAllData, AllData_RAS{ind_meshsize,ind_nmbdig} = SMoutput; end
    
        
        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        ConvFactApprox_RAS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);

        % get the errors when we want them
        indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
        if indDig_CurrOrdrInIndsToPlot ~= 0 
            for ind = 1:length(indsIter_PlotErr)
                u_solseq = SMoutput{1}; 
                u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            % if ~SMoutput{4}, inds_DigsThatStsfyConvCond_RAS = ind_nmbdig:length(RM_nmbdigits_list); CalcErrMtrx_rplc = false; end
            inds_DigsThatStsfyConvCond_RAS = inds_DigsThatStsfyConvCond_dAS;
        end
    end

    % run the MS
    SchwarzMethod = {'MS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; CalcErrMtrx_rplc = false; %if CalcErrMtrx, CalcErrMtrx_rplc = true; else, CalcErrMtrx_rplc = false;  end
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx_rplc};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);
    
        % save the data if we wanna
        if SaveAllData, AllData_MS{ind_meshsize,ind_nmbdig} = SMoutput; end
    
        
        % get the approximate convergence factors
        error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ConvCurves(ind_meshsize,ind_nmbdig,:) = error_nrms; ind_aux = 1;
        while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
        ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
        ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
        ConvFactApprox_MS(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);

        % get the errors when we want them
        indDig_CurrOrdrInIndsToPlot = find(~(indsDigs_PlotErr-RM_nmbdigits_list(ind_nmbdig)));
        if indDig_CurrOrdrInIndsToPlot ~= 0 
            for ind = 1:length(indsIter_PlotErr)
                u_solseq = SMoutput{1}; 
                u_sol = u_solseq(indsIter_PlotErr(ind),:); 
                u_PlotErr(indDig_CurrOrdrInIndsToPlot,ind,:) = reshape(SMoutput{2}'-u_sol, size(u_PlotErr(1,1,:)) );
            end
        end

        % if we did the CalcErrMtrx, then check if we satisfied the condition
        if CalcErrMtrx
            % if ~SMoutput{4}, inds_DigsThatStsfyConvCond_MS = ind_nmbdig:length(RM_nmbdigits_list); CalcErrMtrx_rplc = false; end
            inds_DigsThatStsfyConvCond_MS = inds_DigsThatStsfyConvCond_dAS;
        end
    end

end





if SaveAllData
    MyData = {'Schwarz methods','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'SM #iter', SM_nmbiter,'SM relres accuracy stop', SM_relresacc, 'theta for dAS', dampingTheta, ...
       '# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, 'Err Mtrx conditions calculated', CalcErrMtrx, ...
        'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
        'data dAS', AllData_dAS, 'data RAS', AllData_RAS, 'data MS', AllData_MS, ... 
        'inds that satisfy the conv cond', inds_DigsThatStsfyConvCond_dAS};
    s_SaveString = append(append(append(append('SavedData_mpSM_NonSymm_ConstSiz_Digs_',num2str(RM_nmbdigits_list(1))),'_'),num2str(RM_nmbdigits_list(end))),'.mat');
    save(s_SaveString,'MyData');  
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyColors = {"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F","#000000","#FF0000","#00FF00","#0000FF","#FFFF00","#FF00FF","#00FFFF"};
MyMarkers = {'o','d','s','^','*','>'}; MyLines = {'--','-.',':','-'};



%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; t(1) = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

% left subplot are convergence curves of MS for different d_s
nmb_iters_ToPlot = 100; iters_mesh = 0:nmb_iters_ToPlot;
nmb_DigsToPlot = 6;

ax1 = nexttile();
for ind_nmbdig = nmb_DigsToPlot:-1:1
    PlotData = reshape(ConvCurves(ind_meshsize,ind_nmbdig,1:nmb_iters_ToPlot+1) ,nmb_iters_ToPlot+1,1);
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
plt = plot(RM_nmbdigits_list,ConvFactApprox_dAS(ind_meshsize,:),'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(1) = plt; LegendLabels{1} = 'dAS, $\theta = \frac{1}{3}$';
if CalcErrMtrx, plot(inds_DigsThatStsfyConvCond_dAS,ConvFactApprox_dAS(ind_meshsize,inds_DigsThatStsfyConvCond_dAS),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','d','MarkerSize',14,'LineWidth',2); hold on; end

plt = plot(RM_nmbdigits_list,ConvFactApprox_RAS(ind_meshsize,:),'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(2) = plt; LegendLabels{2} = 'RAS';
if CalcErrMtrx, plot(inds_DigsThatStsfyConvCond_RAS,ConvFactApprox_RAS(ind_meshsize,inds_DigsThatStsfyConvCond_RAS),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','^','MarkerSize',14,'LineWidth',2); hold on; end

plt = plot(RM_nmbdigits_list,ConvFactApprox_MS(ind_meshsize,:),'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on;
% LegendHandles(3) = plt; LegendLabels{3} = 'MS';
if CalcErrMtrx, plot(inds_DigsThatStsfyConvCond_MS,ConvFactApprox_MS(ind_meshsize,inds_DigsThatStsfyConvCond_MS),'MarkerEdgeColor',[MyColors{7}], ...
    'MarkerFaceColor',[MyColors{7}],'Color',[MyColors{7}],'Marker','o','MarkerSize',14,'LineWidth',2); hold on; end

xlabel('$d_s$','FontSize',24, 'interpreter', 'latex'); ylabel('approx. conv. factor', 'FontSize',24,'interpreter', 'latex'); %title('convergence curves for MS', 'interpreter', 'latex');

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


for ind_dig = 1:length(indsDigs_PlotErr)

    %%% first column is the initial error, i.e., the exact solution if we took initial guess as all zeros
    nexttile()
    PlotData = u_plot;
    PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_ExactSol,nmb_int_gridcols,nmb_int_gridcols)'; 
    mesh(x_mesh,x_mesh,PlotData); xlabel('x');ylabel('y'); view(angle1,angle2);
    u_plot = zeros(nmb_int_gridcols+2);
    if ind_dig == 1; myT = 'solution'; title(myT, 'interpreter', 'latex'); end
    myT = append('$d_s$ = ',num2str(indsDigs_PlotErr(ind_dig))); zlabel(myT, 'interpreter', 'latex');


    for ind_iter = 1:length(indsIter_PlotErr)

        nexttile()
        u_err = reshape( u_PlotErr(ind_dig,ind_iter,:), 1,length(rhs));
        PlotData = u_plot;
        PlotData(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = reshape(u_err,nmb_int_gridcols,nmb_int_gridcols)'; 
        mesh(x_mesh,x_mesh,PlotData); view(angle1,angle2); %xlabel('x');ylabel('y');
        if ind_dig == 1; myT = append( append('error after ',num2str(indsIter_PlotErr(ind_iter))), ' it.'); title(myT, 'interpreter', 'latex'); end
        u_plot = zeros(nmb_int_gridcols+2);

    end
end

% font sizes
for ind=1:length(t.Children)
    curr_ax = t.Children(ind);
    set(curr_ax,'FontWeight','bold','FontSize',19,'Color','white');
end


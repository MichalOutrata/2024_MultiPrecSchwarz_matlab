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
indsDigs_PlotErr = [3,4,5,6];
indsMesh_PlotErr = 1;
inds_FrstIndDigThatStsfyConvCond = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 50:10:320; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 1;



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
    
    elseif ProblemChoice == 1 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
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































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
indsDigs_PlotErr = [3,4,5,6];
indsMesh_PlotErr = 1;
inds_FrstIndDigThatStsfyConvCond = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 50:10:320; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 2;



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
    
    elseif ProblemChoice == 1 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
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


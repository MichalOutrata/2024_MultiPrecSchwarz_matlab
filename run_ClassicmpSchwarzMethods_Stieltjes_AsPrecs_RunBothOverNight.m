clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
RM_nmbdigits_list = 1:16; % number of digits to keep for the subdomain solves
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
Advanpix = true; debug = 1; SavePlotData = true;
SM_nmbiter = 42; SM_relresacc = 42.42;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 50:10:100; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 1;


%%% choose GMRES set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GMRES_PrecType = 'L'; % type of preconditioner - none/left/right
GMRES_nmbiter = 100; % number of GMRES iterations
GMRES_relresacc = 1e-12; % GMRES accuracy
Prec = @mpSchwarzMethodsPreconds;





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
    GMRES_initguess = zeros(size(rhs)); 


    % run non-preconditioned GMRES
    PrecInfo = {'noprec'};
    [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, 500, GMRES_relresacc );
    if AddOutputs{4} >= 1 
        conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
        for ind = 1:GMRES_nmbiter
            if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
        end
    end
    GMRES_resnorms_noprec{ind_meshsize} = conv_curves_prettier; GMRES_nmbittoconv_noprec(ind_meshsize) = AddOutputs{4};


    % run the dAS
    SchwarzMethod = {'dAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta};
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   dAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_dAS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_dAS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

    % run the RAS
    SchwarzMethod = {'RAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc};
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   RAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_RAS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_RAS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

    % run the MS
    SchwarzMethod = {'MS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; 
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   MS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_MS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_MS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

end





if SavePlotData
    MyData = {'Schwarz methods as preconds','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'theta for dAS', dampingTheta,'# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, ...
        'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
        'GMRES prec type (noprec/L/R)', GMRES_PrecType, 'max # GMRES iters', GMRES_nmbiter, 'GMRES conv. tolerance', GMRES_relresacc, ...
        'GMRES ConvCrvs noprec', GMRES_resnorms_noprec, 'GMRES # iter no prec', GMRES_nmbittoconv_noprec, ...
        'GMRES ConvCrvs dAS', GMRES_resnorms_dAS, 'GMRES # iter dAS', GMRES_nmbittoconv_dAS, ... 
        'GMRES ConvCrvs RAS', GMRES_resnorms_RAS, 'GMRES # iter RAS', GMRES_nmbittoconv_RAS, ... 
        'GMRES ConvCrvs MS', GMRES_resnorms_MS, 'GMRES # iter MS', GMRES_nmbittoconv_MS};
    s_SaveString = append(append('SavedData_ClassicmpSM_Stieltjes_AsPrecs_Prblm',num2str(ProblemChoice)),'.mat');
    save(s_SaveString,'MyData');

else
    %%% plot the ConvFact as a function of number of digits
    ind_meshsize = 2; % choose
    ClsstIntSqrt = floor(sqrt(length(RM_nmbdigits_list))); nmbColsTiles = floor( length(RM_nmbdigits_list) / ClsstIntSqrt );
    if length(RM_nmbdigits_list) <= 3, nmbColsTiles = length(RM_nmbdigits_list); nmbRowsTiles = 1; else
        if floor( length(RM_nmbdigits_list) / nmbColsTiles ) == length(RM_nmbdigits_list) / nmbColsTiles, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ); 
        else, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ) + 1; end
    end
    figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        ax(ind_nmbdig) = nexttile();
        semilogy(GMRES_resnorms_MS{ind_meshsize,ind_nmbdig},'ko-');
    end
    title(t, append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));
    
    figure(97)
    plot(RM_nmbdigits_list,GMRES_nmbittoconv_MS(ind_meshsize,:),'ro-'); xlabel('# digits');ylabel('# GMRES it to conv'); title( append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));
    
    %%% plot the ConvFact as a function of mesh size
    ind_nmbdig = 3; % choose
    ClsstIntSqrt = floor(sqrt(length(list_of_nmb_int_gridcols))); nmbColsTiles = floor( length(list_of_nmb_int_gridcols) / ClsstIntSqrt );
    if length(list_of_nmb_int_gridcols) <= 3, nmbColsTiles = length(list_of_nmb_int_gridcols); nmbRowsTiles = 1; else
        if floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ) == length(list_of_nmb_int_gridcols) / nmbColsTiles, nmbRowsTiles = floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ); 
        else, nmbRowsTiles = floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ) + 1; end
    end
    figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');
    for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
        ax(ind_meshsize) = nexttile(); 
        semilogy(GMRES_resnorms_MS{ind_meshsize,ind_nmbdig},'ko-');
    end
    title(t, append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));
    
    figure(98)
    plot(list_of_nmb_int_gridcols,GMRES_nmbittoconv_MS(:,ind_nmbdig),'ro-'); xlabel('# grid points on one side');ylabel('# GMRES it to conv'); title( append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));
    
    %%% plot the ConvFact heat map
    figure(99)
    imagesc(GMRES_nmbittoconv_MS); colorbar;
end




































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
RM_nmbdigits_list = 1:16; % number of digits to keep for the subdomain solves
dampingTheta = 1/3; % generally one should use 1/(2^SM_nmbsubdoms_PwrOfTwo+1) but since we have sausage-like domain decomposition, we can do with alternating coloring no matter how many subdomains
Advanpix = true; debug = 1; SavePlotData = true;
SM_nmbiter = 42; SM_relresacc = 42.42;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 50:10:100; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 2;


%%% choose GMRES set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GMRES_PrecType = 'L'; % type of preconditioner - none/left/right
GMRES_nmbiter = 100; % number of GMRES iterations
GMRES_relresacc = 1e-12; % GMRES accuracy
Prec = @mpSchwarzMethodsPreconds;





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
    GMRES_initguess = zeros(size(rhs)); 


    % run non-preconditioned GMRES
    PrecInfo = {'noprec'};
    [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, 500, GMRES_relresacc );
    if AddOutputs{4} >= 1 
        conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
        for ind = 1:GMRES_nmbiter
            if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
        end
    end
    GMRES_resnorms_noprec{ind_meshsize} = conv_curves_prettier; GMRES_nmbittoconv_noprec(ind_meshsize) = AddOutputs{4};


    % run the dAS
    SchwarzMethod = {'dAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta};
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   dAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_dAS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_dAS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

    % run the RAS
    SchwarzMethod = {'RAS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc};
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   RAS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_RAS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_RAS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

    % run the MS
    SchwarzMethod = {'MS', 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; 
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        disp('   -------------------')
        disp(append('   MS with d_s = ',num2str(RM_nmbdigits_list(ind_nmbdig))))

        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, false};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug); PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );

        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = conv_curves;
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) < GMRES_relresacc, conv_curves_prettier(ind:end) = nan; break; end
            end
        end
        GMRES_resnorms_MS{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv_MS(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    end

end





MyData = {'Schwarz methods as preconds','dAS,RAS,MS', '# subdoms = 2^{...}',SM_nmbsubdoms_PwrOfTwo,'theta for dAS', dampingTheta,'# digits kept', RM_nmbdigits_list, 'advanpix used', Advanpix, ...
    'problem index', ProblemChoice, '# interior grid points on the side of the square',list_of_nmb_int_gridcols, ...
    'GMRES prec type (noprec/L/R)', GMRES_PrecType, 'max # GMRES iters', GMRES_nmbiter, 'GMRES conv. tolerance', GMRES_relresacc, ...
    'GMRES ConvCrvs noprec', GMRES_resnorms_noprec, 'GMRES # iter no prec', GMRES_nmbittoconv_noprec, ...
    'GMRES ConvCrvs dAS', GMRES_resnorms_dAS, 'GMRES # iter dAS', GMRES_nmbittoconv_dAS, ... 
    'GMRES ConvCrvs RAS', GMRES_resnorms_RAS, 'GMRES # iter RAS', GMRES_nmbittoconv_RAS, ... 
    'GMRES ConvCrvs MS', GMRES_resnorms_MS, 'GMRES # iter MS', GMRES_nmbittoconv_MS};
s_SaveString = append(append('SavedData_ClassicmpSM_Stieltjes_AsPrecs_Prblm',num2str(ProblemChoice)),'.mat');
save(s_SaveString,'MyData');


clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_type = 'MS'; % type of Schwarz method
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
SM_nmbiter = 50; % number of RAS iterations
SM_relresacc = 1e-12; % number of RAS iterations
RM_nmbdigits_list = 2:8; % number of digits to keep for the subdomain solves
dampingTheta = 1/(SM_nmbsubdoms_PwrOfTwo+1);
Advanpix = false; CalcErrMtrx = true; DoPlotting = true;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 170:5:180; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 2;



for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
    nmb_int_gridcols = list_of_nmb_int_gridcols(ind_meshsize); h=1/(nmb_int_gridcols+1); % mesh size
    
    if ProblemChoice == 0 %%% do negative laplacian
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) zeros(size(x)); a = @(x,y) ones(size(x)); b1 = @(x,y) zeros(size(x)); b2 = b1;
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1); rhs = 1/h^2*rhs; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
        u_init = zeros(size(A,1),1); u_plot = zeros(nmb_int_gridcols+2); %u_plot(2:nmb_int_gridcols+1,1) = 1; 
        angle1 = 52.5; angle2 = 30;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';
    
    elseif ProblemChoice == 1 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; alpha =  @(x,y) (x+y).^2.*exp(x-y); nu = @(x,y) 0*(1.*x+1.*y); mu = @(x,y) 0*(1.*x+1.*y);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,nu,mu});
        rhs = ones(size(A,1),1)*1/sqrt(nmb_int_gridcols^2)*1/h^2; u_init = randn(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';
    
    elseif ProblemChoice == 2 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2); 
        eta = @(x,y) 500.*x+1.*y; alpha = @(x,y) 1+9.*(x+y); beta = @(x,y) 1+9.*(x+y); nu = @(x,y) 0.*x+0.*y; mu = @(x,y) 0.*x+0.*y;
        A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
        rhs = ones(size(A,1),1)*1/sqrt(nmb_int_gridcols^2)*1/h^2; u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;
        RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts';
    
    elseif ProblemChoice == 3 %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) x.^2.*cos(x+y).^2'; a =  @(x,y) (x+y).^2.*exp(x-y); b1 =  @(x,y) (y-0.5); b2 =  @(x,y) -(x-0.5);
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; RM_type = 'Mmtrx';
    
    elseif ProblemChoice == 4 %%% do non-symmetric AdvecDiff based SzyldFrommer eqn (15) and p.660 -> changed the advection strength to ~40 so that we have an Mmtrx
        G=numgrid('S',nmb_int_gridcols+2); AdvectStrngth = 100;
        eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1 + 0.*x+0.*y; beta = @(x,y) 1 + 0.*x+0.*y; nu = @(x,y) AdvectStrngth.*(1.*x.*(x-1).*(1-2.*y)); mu = @(x,y) AdvectStrngth.*(-1.*y.*(y-1).*(1-2.*x));
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,alpha,beta,nu,mu}); %A_szyld = A*h^2+speye(size(A,1)); disp(isMmtrx(A)); disp(isMmtrx(A_szyld)); disp(issymmetric(A))
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; RM_type = 'Mmtrx';
    
    % %%% optional more advanced rhs building
    % [m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h); IntriorInds = find(G);
    % rhsfun = @(x,y) x.*y.*(1-x).*(1-y); rhsfun_discr = feval(rhsfun,X,Y); rhs = A*rhsfun_discr(IntriorInds);
    
    end
    
    
    
    if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; end
    debug = 1; if DoPlotting, MyPlot = {0:h:1,u_plot,0.05,angle1,angle2}; else, MyPlot = {}; end % MyPlot = {x,u_plot,waitime,angle1,angle2};
    
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx};
        [SMoutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, MyPlot, debug);
    
        if debug ~= 0
            error_nrms = SMoutput{3}; error_nrms = error_nrms / error_nrms(1); ind_aux = 1;
            while ind_aux < length(error_nrms) && error_nrms(ind_aux+1) > 1e-15, ind_aux = ind_aux + 1; end; ind_EndOfConv(ind_meshsize,ind_nmbdig) = ind_aux; %if ind_aux > 20, ind_aux = ind_aux - 5; end
            ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) = error_nrms(2:ind_EndOfConv(ind_meshsize,ind_nmbdig))./error_nrms(1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1);
            ConvFactApproxSequence(ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) = ( ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-2) + ConsecErrsRatio(ind_meshsize,ind_nmbdig,2:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1) ) / 2;
            ConvFactApprox(ind_meshsize,ind_nmbdig) = ConvFactApproxSequence(ind_nmbdig,ind_EndOfConv(ind_meshsize,ind_nmbdig)-2);
        end
    end
    %%% we assume that the overlap is corresponding to one matrix block,
    %%% i.e., to one grid-column, i.e., the overlap is 2*h
    %ConvFactContinuous(ind_meshsize)
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_meshsize = 1; % choose
ClsstIntSqrt = floor(sqrt(length(RM_nmbdigits_list))); nmbColsTiles = floor( length(RM_nmbdigits_list) / ClsstIntSqrt );
if length(RM_nmbdigits_list) <= 3, nmbColsTiles = length(RM_nmbdigits_list); nmbRowsTiles = 1; else
    if floor( length(RM_nmbdigits_list) / nmbColsTiles ) == length(RM_nmbdigits_list) / nmbColsTiles, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ); 
    else, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ) + 1; end
end
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

for ind_nmbdig = 1:length(RM_nmbdigits_list)
    ax(ind_nmbdig) = nexttile();
    PlotData = reshape(ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1), 1,ind_EndOfConv(ind_meshsize,ind_nmbdig)-1 );
    plot(PlotData,'ko-'); hold on; plot(ones(ind_EndOfConv(ind_meshsize,ind_nmbdig)-1,1)*ConvFactApprox(ind_meshsize,ind_nmbdig),'r-'); hold off;
end
title(t, append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));

figure(97)
plot(RM_nmbdigits_list,ConvFactApprox(ind_meshsize,:),'ro-'); xlabel('# digits');ylabel('observed ConvFact'); title( append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));




%%% plot the ConvFact as a function of mesh size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_nmbdig = 3; % choose
ClsstIntSqrt = floor(sqrt(length(list_of_nmb_int_gridcols))); nmbColsTiles = floor( length(list_of_nmb_int_gridcols) / ClsstIntSqrt );
if length(list_of_nmb_int_gridcols) <= 3, nmbColsTiles = length(list_of_nmb_int_gridcols); nmbRowsTiles = 1; else
    if floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ) == length(list_of_nmb_int_gridcols) / nmbColsTiles, nmbRowsTiles = floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ); 
    else, nmbRowsTiles = floor( length(list_of_nmb_int_gridcols) / nmbColsTiles ) + 1; end
end
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');

for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
    ax(ind_meshsize) = nexttile(); 
    PlotData = reshape(ConsecErrsRatio(ind_meshsize,ind_nmbdig,1:ind_EndOfConv(ind_meshsize,ind_nmbdig)-1), 1,ind_EndOfConv(ind_meshsize,ind_nmbdig)-1 );
    plot(PlotData,'ko-'); hold on; plot(ones(ind_EndOfConv(ind_meshsize,ind_nmbdig)-1,1)*ConvFactApprox(ind_meshsize,ind_nmbdig),'r-'); hold off;
end
title(t, append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));

figure(98)
plot(list_of_nmb_int_gridcols,ConvFactApprox(:,ind_nmbdig),'ro-'); xlabel('# grid points on one side');ylabel('observed ConvFact'); title( append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));






%%% plot the ConvFact heat map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(99)
imagesc(ConvFactApprox); colorbar;




clear; clc; close('all');


%%% choose Schwarz Method set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_type = 'RAS'; % type of Schwarz method
SM_nmbsubdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
RM_nmbdigits_list = 1:12; % number of digits to keep for the subdomain solves
dampingTheta = 1/(SM_nmbsubdoms_PwrOfTwo+1);
Advanpix = false; CalcErrMtrx = false; SandwichScaling = false; DoPlotting = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_of_nmb_int_gridcols = 10:10:160; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
ProblemChoice = 4;


%%% choose GMRES set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GMRES_PrecType = 'L'; % type of preconditioner - left/right
GMRES_nmbiter = 50; % number of GMRES iterations
GMRES_relresacc = 1e-15; % GMRES accuracy
GMRES_resnorms = cell(length(list_of_nmb_int_gridcols),length(RM_nmbdigits_list)); GMRES_nmbittoconv = nan(length(list_of_nmb_int_gridcols),length(RM_nmbdigits_list));
Prec = @mpSchwarzMethodsPreconds;



for ind_meshsize = 1:length(list_of_nmb_int_gridcols)
    nmb_int_gridcols = list_of_nmb_int_gridcols(ind_meshsize); h=1/(nmb_int_gridcols+1); % mesh size
    
    if ProblemChoice == 1 %%% do negative laplacian
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) zeros(size(x)); a = @(x,y) ones(size(x)); b1 = @(x,y) zeros(size(x)); b2 = b1;
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1); rhs = 1/h^2*rhs; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
        u_init = zeros(size(A,1),1); u_plot = zeros(nmb_int_gridcols+2); %u_plot(2:nmb_int_gridcols+1,1) = 1; 
        angle1 = 52.5; angle2 = 30;
        if SandwichScaling, RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts'; else, RM_type = 'Stieltjess_RoundBi_Facts'; end
    
    elseif ProblemChoice == 2 %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
        G=numgrid('S',nmb_int_gridcols+2);
        eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1+9.*(x+y); beta = @(x,y) 1+9.*(x+y); nu = @(x,y) 0.*x+0.*y; mu = @(x,y) 0.*x+0.*y;
        A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;
        if SandwichScaling, RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts'; else, RM_type = 'Stieltjess_RoundBi_Facts'; end
    
    elseif ProblemChoice == 3 %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
        G=numgrid('S',nmb_int_gridcols+2);
        eta=inline('x.^2.*cos(x+y).^2','x','y'); a=inline('(x+y).^2.*exp(x-y)','x','y'); b1=inline('(y-0.5)','x','y'); b2=inline('-(x-0.5)','x','y');
        A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; RM_type = 'Mmtrx';
    
    elseif ProblemChoice == 4 %%% do non-symmetric AdvecDiff based SzyldFrommer eqn (15) and p.660 -> changed the advection strength to ~40 so that we have an Mmtrx
        G=numgrid('S',nmb_int_gridcols+2); AdvectStrngth = 185;
        eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1 + 0.*x+0.*y; beta = @(x,y) 1 + 0.*x+0.*y; nu = @(x,y) AdvectStrngth.*(1.*x.*(x-1).*(1-2.*y)); mu = @(x,y) AdvectStrngth.*(-1.*y.*(y-1).*(1-2.*x));
        A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu}); %A_szyld = A*h^2+speye(size(A,1)); disp(isMmtrx(A)); disp(isMmtrx(A_szyld)); disp(issymmetric(A))
        rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
        u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; RM_type = 'Mmtrx';
    
    % %%% optional more advanced rhs building
    % [m,n] = size(G); h=1/(n-1); [X,Y]=meshgrid(0:h:1,0:h:(m-1)*h); IntriorInds = find(G);
    % rhsfun = @(x,y) x.*y.*(1-x).*(1-y); rhsfun_discr = feval(rhsfun,X,Y); rhs = A*rhsfun_discr(IntriorInds);
    
    end


    if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo}; end; debug = 0;
    GMRES_initguess = zeros(size(rhs)); 
    
    for ind_nmbdig = 1:length(RM_nmbdigits_list)
    
        RoundingMethod = {RM_type, RM_nmbdigits_list(ind_nmbdig), Advanpix, CalcErrMtrx};
        [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug);
        PrecInfo = {GMRES_PrecType,'Handle',{A,MyPrecondPackage}};
    
        [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, rhs, PrecInfo, Prec, GMRES_initguess, GMRES_nmbiter, GMRES_relresacc );
        if AddOutputs{8} >= 1 
            conv_curves = res_norms/res_norms(1); conv_curves_prettier = 1e-16*ones(GMRES_nmbiter,1);
            for ind = 1:GMRES_nmbiter
                if ind <= length(conv_curves) && conv_curves(ind) > 1e-16, conv_curves_prettier(ind) = conv_curves(ind); end
            end
        end
        GMRES_resnorms{ind_meshsize,ind_nmbdig} = conv_curves_prettier; GMRES_nmbittoconv(ind_meshsize,ind_nmbdig) = AddOutputs{8};
    
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plotting     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% plot the ConvFact as a function of number of digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_meshsize = 2; % choose
ClsstIntSqrt = floor(sqrt(length(RM_nmbdigits_list))); nmbColsTiles = floor( length(RM_nmbdigits_list) / ClsstIntSqrt );
if length(RM_nmbdigits_list) <= 3, nmbColsTiles = length(RM_nmbdigits_list); nmbRowsTiles = 1; else
    if floor( length(RM_nmbdigits_list) / nmbColsTiles ) == length(RM_nmbdigits_list) / nmbColsTiles, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ); 
    else, nmbRowsTiles = floor( length(RM_nmbdigits_list) / nmbColsTiles ) + 1; end
end
figure; t(1) = tiledlayout(nmbRowsTiles,nmbColsTiles,'TileSpacing','Compact','Padding','Compact');


for ind_nmbdig = 1:length(RM_nmbdigits_list)
    ax(ind_nmbdig) = nexttile();
    semilogy(GMRES_resnorms{ind_meshsize,ind_nmbdig},'ko-');
end
title(t, append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));

figure(97)
plot(RM_nmbdigits_list,GMRES_nmbittoconv(ind_meshsize,:),'ro-'); xlabel('# digits');ylabel('# GMRES it to conv'); title( append('# grid points on one side = ',num2str(list_of_nmb_int_gridcols(ind_meshsize))));




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
    semilogy(GMRES_resnorms{ind_meshsize,ind_nmbdig},'ko-');
end
title(t, append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));

figure(98)
plot(list_of_nmb_int_gridcols,GMRES_nmbittoconv(:,ind_nmbdig),'ro-'); xlabel('# grid points on one side');ylabel('# GMRES it to conv'); title( append('# digits = ',num2str(RM_nmbdigits_list(ind_nmbdig))));






%%% plot the ConvFact heat map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(99)
imagesc(GMRES_nmbittoconv); colorbar;









% 
% for ind_nmbdig = 1:length(RM_nmbdigits_list)
%     ax(ind_nmbdig) = nexttile(); 
%     semilogy(1:GMRES_nmbiter,GMRES_resnorms{ind_nmbdig},'ko-'); %hold on; plot(ones(ind_EndOfConv(ind_nmbdig)-1,1)*ConvFactApprox(ind_nmbdig),'r-'); hold off;
% end
% 
% % figure(99)
% % plot(RM_nmbdigits_list,ConvFactApprox,'ro-'); xlabel('# digits');ylabel('observed ConvFact');



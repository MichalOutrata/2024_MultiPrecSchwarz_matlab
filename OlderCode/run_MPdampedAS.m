clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb_subdoms_PwrOfTwo = 2; % number of subdomains for the Schwarz method
nmb_subdoms = 2^nmb_subdoms_PwrOfTwo; % number of subdomains for the Schwarz method
d_s = 8; % number of digits to keep for the subdomain solves
nmb_iter = 100; % number of RAS iterations
waitime = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb_int_gridcols = 20; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
damping_theta = 1/nmb_subdoms; % damping parameter for damped additive Schwarz
h=1/(nmb_int_gridcols+1); % mesh size
N = nmb_int_gridcols^2; % total nmb of unknowns

% %%% do negative laplacian
% G=numgrid('S',nmb_int_gridcols+2);
% eta = @(x,y) zeros(size(x)); a = @(x,y) ones(size(x));
% b1 = @(x,y) zeros(size(x)); b2 = b1;
% A=delsqnonsymmetric_Gander(eta,a,b1,b2,G);
% RoundType = 'Stieltjess';
% rhs = zeros(size(A,1),1); rhs(1:nmb_int_gridcols:end)=-(-1); rhs = 1/h^2*rhs; % right-hand side, one minus for shifting to the RHS, one minus for having negative laplacian
% u_init = zeros(size(A,1),1); u_plot = zeros(nmb_int_gridcols+2); u_plot(2:nmb_int_gridcols+1,1) = 1; angle1 = 52.5; angle2 = 30;


% 
% %%% do non-symmetric AdvecDiff
% G=numgrid('S',nmb_int_gridcols+2);
% eta = @(x,y) x.^2.*cos(x+y).^2; a = @(x,y) (x+y).^2.*exp(x-y);
% b1 = @(x,y) (y-0.5); b2 = @(x,y) -(x-0.5);
% A=delsqnonsymmetric_Gander(eta,a,b1,b2,G);
% RoundType = 'Mmtrx';
% rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
% u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;


if size(A,1) ~= N, disp('error, if the matrix was FD'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,nmb_subdoms_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(nmb_subdoms,1); last_rows = NaN(nmb_subdoms,1); 
frst_cols = NaN(nmb_subdoms,1); last_cols = NaN(nmb_subdoms,1);
for ind_subdom = 1:nmb_subdoms
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end
    
    if ind_subdom ~= nmb_subdoms
        last_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +3 ); last_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +3 );
    else
        last_rows(nmb_subdoms) = N; last_cols(nmb_subdoms) = N;
    end
end
% plotPartition_Gander( A, ni, nj ); % check for indices


%%% do the LowPrecision blocks
if strcmp(RoundType,'Stieltjess')
    MySubdomFactors_LP = cell(5,nmb_subdoms); % need to keep the block diagonal as well as L,U,P,Q
elseif strcmp(RoundType,'Mmtrx')
    MySubdomFactors_LP = cell(4,nmb_subdoms); % need to keep the L,U,P,Q
end
 
for ind_subdom = 1:nmb_subdoms
    if strcmp(RoundType,'Stieltjess')
        [curr_L,curr_U,curr_Pvec,curr_Qvec_inv, Ai_diag] = MtrxFactorsLP_Stieltjess( A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)), d_s );
        MySubdomFactors_LP{1,ind_subdom} = curr_L; MySubdomFactors_LP{2,ind_subdom} = curr_U; MySubdomFactors_LP{3,ind_subdom} = curr_Pvec; MySubdomFactors_LP{4,ind_subdom} = curr_Qvec_inv; MySubdomFactors_LP{5,ind_subdom} = Ai_diag;
    elseif strcmp(RoundType,'Mmtrx')
        [curr_L,curr_U,curr_Pvec,curr_Qvec_inv] = MtrxFactorsLP_Mmtrx( A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)), d_s );
        MySubdomFactors_LP{1,ind_subdom} = curr_L; MySubdomFactors_LP{2,ind_subdom} = curr_U; MySubdomFactors_LP{3,ind_subdom} = curr_Pvec; MySubdomFactors_LP{4,ind_subdom} = curr_Qvec_inv;
    end
end








% %%% run damped AS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 u_exact = A\rhs; u_dd = u_init; res_prev = rhs-A*u_dd; 

 u_err = zeros(nmb_iter+1,N); err_norm = zeros(nmb_iter+1,1);
 u_err(1,:) = u_exact - u_init; err_norm(1) = norm(u_err(1,:));

% check convergence
u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)';
x=0:h:1; mesh(x,x,u_plot)
xlabel('x');ylabel('y'); view(angle1,angle2); pause(2);


% % check for indices
% figure(2); plotspecs = {'bo','ro','go','ko','yo','mo','co','bo','ro','go'};


for ind_iter = 1:nmb_iter
    for ind_subdom = 1:nmb_subdoms

        curr_L = MySubdomFactors_LP{1,ind_subdom}; curr_U = MySubdomFactors_LP{2,ind_subdom}; curr_Pvec = MySubdomFactors_LP{3,ind_subdom}; curr_Qvec_inv = MySubdomFactors_LP{4,ind_subdom}; 
        curr_rhs_bfrperm = res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));
        if strcmp(RoundType,'Stieltjess')
            Ai_diag = MySubdomFactors_LP{5,ind_subdom}; curr_rhs_bfrperm = 1 ./ Ai_diag .* curr_rhs_bfrperm;
        end
        curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_U \ ( curr_L \ curr_rhs); curr_u = sol_perm(curr_Qvec_inv);
        % %%%% check
        % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); curr_u_check = Ai \ res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));

        if ind_subdom == 1
            u_dd = u_dd + damping_theta * [curr_u; zeros(N-last_rows(ind_subdom),1)];
            % check convergence
            u_plot_new=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = u_plot_new; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);

        elseif ind_subdom == nmb_subdoms
            u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u]; 
            % check convergence
            u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);

        else
            u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u; zeros(N-last_rows(ind_subdom),1)];
            % check convergence
            u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);
        end
    end
    res_prev = rhs-A*u_dd; u_err(ind_iter+1,:) = u_exact - u_dd; err_norm(ind_iter+1) = norm(u_err(ind_iter+1,:));
end

semilogy([0,1:nmbiter], err_norm, 'ro-');










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[L,U,p,q_inv] = MtrxFactorsLP_Mmtrx(A,d_s)
    if ~issparse(A), disp('error, input not sparse'); end
    [i,j,v] = find(A); [v_sorted,perm] = sort(v); perm_inv(perm) = 1:numel(perm);
    
    [~,ind_AbsValMin] = min(abs(v_sorted));
    if v_sorted(ind_AbsValMin) > 0 % if we have a positive value being smallest in magn. -> the entry bfr the first occurence (=index "ind_AbsValMin") is the smallest negative entry in abs. val.
        ind_AbsValMin = ind_AbsValMin - 1;
    else % if we have a negative value being smallest in magn. -> if it is repeated, the frst occurence (=index "ind_AbsValMin") is not what we want, we want the last one (so that the next entry is positive)
        while v_sorted(ind_AbsValMin+1) < 0
            ind_AbsValMin = ind_AbsValMin + 1;
        end
    end
    
    % round positive entries
    pwrs = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; pwrs_pos = max(pwrs);
    v_PosRounded = ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos))./10.^(d_s-pwrs_pos);
    
    % round negative entries
    pwrs = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; pwrs_neg = max(pwrs);
    v_NegRounded = ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    
    v_Rounded = [v_NegRounded;v_PosRounded]; A_rounded_aux = sparse(i,j,v_Rounded(perm_inv)); A_rounded = mp(A_rounded_aux);
    [L,U,p,q] = lu(A_rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    % err = A - A_rounded; if ~isempty(err(err<0)), disp('error'); end
end
function[L,U,p,q_inv,A_diag] = MtrxFactorsLP_Stieltjess(A,d_s)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A);
    % bcs A was Stieltjess, we know that the diagonal is nonzero
    
    [i,j,v] = find(A_DiagRplcZeros); 
    
    % round off-diag entries (they're all engative)
    pwrs = floor(log10(-v)) +1; pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    
    A_rounded = sparse(i,j,v_Rounded); A_rounded_OnesDiag_aux = spdiags(ones(size(A,1),1),0, diag(A_diag)\A_rounded);
    A_rounded_OnesDiag = mp(A_rounded_OnesDiag_aux);
    [L,U,p,q] = lu(A_rounded_OnesDiag,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse
end

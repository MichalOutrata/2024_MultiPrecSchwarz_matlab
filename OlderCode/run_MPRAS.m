clear; clc; close('all');

%%% choose set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb_subdoms_PwrOfTwo = 1; % number of subdomains for the Schwarz method
nmb_subdoms = 2^nmb_subdoms_PwrOfTwo; % number of subdomains for the Schwarz method
d_s = 4; % number of digits to keep for the subdomain solves
nmb_iter = 100; % number of RAS iterations
waitime = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmb_int_gridcols = 20; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
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



%%% do non-symmetric AdvecDiff
G=numgrid('S',nmb_int_gridcols+2);
eta = @(x,y) x.^2.*cos(x+y).^2; a = @(x,y) (x+y).^2.*exp(x-y);
b1 = @(x,y) (y-0.5); b2 = @(x,y) -(x-0.5);
A=delsqnonsymmetric_Gander(eta,a,b1,b2,G);
RoundType = 'Mmtrx';
rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;


if size(A,1) ~= N, disp('error, if the matrix was FD'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,nmb_subdoms_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(nmb_subdoms,1); mid_rows = NaN(nmb_subdoms,1); last_rows = NaN(nmb_subdoms,1); 
frst_cols = NaN(nmb_subdoms,1); mid_cols = NaN(nmb_subdoms,1); last_cols = NaN(nmb_subdoms,1);
for ind_subdom = 1:nmb_subdoms
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end

    if ind_subdom ~= nmb_subdoms
        mid_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +2 ); mid_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +2 );
    end
    
    if ind_subdom ~= nmb_subdoms
        last_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +3 ); last_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +3 );
    else
        last_rows(nmb_subdoms) = N; last_cols(nmb_subdoms) = N;
    end
end
frst_rows_RAS = NaN(nmb_subdoms,1); last_rows_RAS = NaN(nmb_subdoms,1);
for ind_subdom = 1:nmb_subdoms
    if ind_subdom == 1
        frst_rows_RAS(ind_subdom) = 1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
    elseif ind_subdom == nmb_subdoms
        frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = N;
    else
        frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
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








% %%% run RAS
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
                u_dd = u_dd + [curr_u(1:last_rows_RAS(ind_subdom)); zeros(N-last_rows_RAS(ind_subdom),1)];
                
                % % check for indices
                % disp(append('updated unknowns: ',append(num2str(1),append('-',num2str(last_rows_RAS(ind_subdom))))))
                % spy(u_dd,plotspecs{ind_subdom}); hold on;
    
    
                % check convergence
                u_plot_new=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = u_plot_new; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);
    
            elseif ind_subdom == nmb_subdoms
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1;
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:end)]; 
                
                % % check for indices
                % disp(append('updated unknowns: ',append(num2str(frst_rows(ind_subdom) + frst_row_inside_subdom-1),append('-',num2str(N)))))
                % spy([zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:end)],plotspecs{ind_subdom}); hold on;
    
    
                % check convergence
                u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);
    
            else
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_in_left_overlap = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                last_row_inside_subdom = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_to_right_overlap = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)];
    
                % % check for indices
                % disp(append('updated unknowns: ',append(num2str(frst_rows(ind_subdom) + frst_row_inside_subdom-1),append('-',num2str(frst_rows(ind_subdom) + last_row_inside_subdom-1)))))
                % spy([zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)],plotspecs{ind_subdom}); hold on; 
    
    
                % check convergence
                u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);
            end
    end
    res_prev = rhs-A*u_dd; u_err(ind_iter+1,:) = u_exact - u_dd; err_norm(ind_iter+1) = norm(u_err(ind_iter+1,:));
end

semilogy([0,1:nmbiter], err_norm, 'ro-');
















% %%% run Schwarz method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% check AS
% R1=[speye(last_row_A1) , sparse(last_row_A1, dim_A-last_row_A1 )];
% R2=[sparse(dim_A-frst_row_A2+1, frst_row_A2-1 ) , speye(dim_A-frst_row_A2+1) ];
% A1_check=R1*A*R1'; err_1 = A1_check - A1; check_1 = isempty( err_1(err_1~=0) );
% A2_check=R2*A*R2'; err_2 = A2_check - A2; check_1 = isempty( err_2(err_2~=0) );
% 
% %%% check RAS
% R1_tilde = R1; R1_tilde(ind_MidInA2:end,:) = 0;
% R2_tilde = R2; R2_tilde(1:ind_MidInA2-frst_row_A2,:) = 0;
% 
% 
% 
% %%% calculate exact soltuion
% u = A\rhs;
% 
% u_dd=zeros(size(rhs)); 
% 
% %%% plotting
% N_int_gridlength = fix(sqrt(size(A,1)));
% u_plot=zeros(N_int_gridlength+2); u_plot(1:N_int_gridlength+2,1)=1;
% 
% 
% 
% for i=1:10
%     %%% plotting 
%     u_plot(2:N_int_gridlength+1,2:N_int_gridlength+1)=reshape(u_dd,N_int_gridlength,N_int_gridlength)';
%     h=1/(N_int_gridlength+1); x=0:h:1;
%     mesh(x,x,u_plot)
%     xlabel('x');ylabel('y');
%     view(52.5,30)
%     %pause
% 
%     %%% calculating eror
%     err(i) = max(max(abs(u-u_dd)));
% 
%     % %%% run RAS
%     % res_prev = rhs-A*u_dd; 
%     % u_A1sol = A1\res_prev(1:last_row_A1); u_A2sol = A2\res_prev(frst_row_A2:end);
%     % u_dd = u_dd + [u_A1sol(1:ind_MidInA1); zeros(dim_A-ind_MidInA1,1)] + [zeros(dim_A-ind_MidInA1,1) ; u_A2sol(1+ ind_MidInA2-frst_row_A2:end)];
%     % u_dd_check = u_dd + R1_tilde'*(A1\(R1*res_prev)) + R2_tilde'*(A2\(R2*res_prev));
% 
%     %%% run AS
%     res_prev = rhs-A*u_dd;
%     u_dd = u_dd + [A1\res_prev(1:last_row_A1); zeros(dim_A-last_row_A1,1)] + [zeros(dim_A-last_row_A1,1) ; A2\res_prev(frst_row_A2:end)];
%     u_dd_check = u_dd + R1'*(A1\(R1*res_prev)) + R2'*(A2\(R2*res_prev));
% 
%     %%% run MS
%     res_prev = rhs-A*u_dd;
%     u_dd_half = u_dd + [A1\res_prev(1:last_row_A1); zeros(dim_A-last_row_A1,1)]; res_half = rhs-A*u_dd;
%     u_dd = u_dd_half + [zeros(dim_A-last_row_A1,1) ; A2\res_half(frst_row_A2:end)];
%     u_dd_check_half = u_dd + R1'*(A1\(R1*res_prev)); res_check_half = rhs-A*u_dd_check_half;
%     u_dd_check = u_dd + R2'*(A2\(R2*res_check_half));
% end













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

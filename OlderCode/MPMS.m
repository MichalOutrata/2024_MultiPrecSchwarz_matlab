function[u_seq] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, plot, debug)


SM_type = SchwarzMethod{1}; SM_nmbsubdom_input = SchwarzMethod{2}; SM_nmbiter = SchwarzMethod{3};
RM_type = RoundingMethod{1}; RM_nmbdigits = RoundingMethod{2};

N = length(rhs);
SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo;
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end


if ~isempty(plot)
    waitime = 0.05;
end


%%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,SM_nmbsubdom_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(SM_nmbsubdom,1); last_rows = NaN(SM_nmbsubdom,1);  frst_cols = NaN(SM_nmbsubdom,1); last_cols = NaN(SM_nmbsubdom,1);
for ind_subdom = 1:SM_nmbsubdom
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end
 
    if ind_subdom ~= SM_nmbsubdom
        last_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +3 ); last_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +3 );
    else
        last_rows(nmb_subdoms) = N; last_cols(nmb_subdoms) = N;
    end
end
if debug == 1, plotPartition_Gander( A, ni, nj ); end % check for indices

%%% do the LowPrecision blocks
MySubdomFactors_LP = cell(nmb_subdoms,1); % need to keep the L,U,P,Q as well as other things (eg, A_diag for Stieltjess-type rounding)

for ind_subdom = 1:SM_nmbsubdom
    SubdomPrblmFacts = SubdomProbRounding( A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)), RM_nmbdigits, RM_type);
    MySubdomFactors_LP{ind_subdom} = SubdomPrblmFacts;
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



%     res_prev = rhs-A*u_dd;
%     u_dd_half = u_dd + [A1\res_prev(1:last_row_A1); zeros(dim_A-last_row_A1,1)]; res_half = rhs-A*u_dd;
%     u_dd = u_dd_half + [zeros(dim_A-last_row_A1,1) ; A2\res_half(frst_row_A2:end)];
%     u_dd_check_half = u_dd + R1'*(A1\(R1*res_prev)); res_check_half = rhs-A*u_dd_check_half;
%     u_dd_check = u_dd + R2'*(A2\(R2*res_check_half));



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
            u_dd = u_dd + [curr_u; zeros(N-last_rows(ind_subdom),1)];
            % check convergence
            u_plot_new=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1) = u_plot_new; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);

        elseif ind_subdom == nmb_subdoms
            u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u]; 
            % check convergence
            u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);

        else
            u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u; zeros(N-last_rows(ind_subdom),1)];
            % check convergence
            u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x,x,u_plot); view(angle1,angle2); pause(waitime);
        end
        res_prev = rhs-A*u_dd; 
    end
    u_err(ind_iter+1,:) = u_exact - u_dd; err_norm(ind_iter+1) = norm(u_err(ind_iter+1,:));
end

semilogy(err_norm,'ro-')



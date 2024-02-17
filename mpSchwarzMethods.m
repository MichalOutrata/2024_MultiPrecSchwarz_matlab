function[MyOutput] = mpSchwarzMethods(A,rhs,u_init, SchwarzMethod,RoundingMethod, plot,debug)

    function v_out = MatVecFun(v_in,tflag, Msolve,Mmult)
    if strcmp(tflag,'notransp'), v_out_aux = Mmult * v_in; v_out = Msolve \ v_out_aux;
    else, v_out_aux = Msolve' \ v_in; v_out = Mmult' * v_out_aux; end
    end

SM_type = SchwarzMethod{1}; SM_nmbsubdom_input = SchwarzMethod{2}; SM_nmbiter = SchwarzMethod{3}; SM_relresacc = SchwarzMethod{4}; if strcmp(SM_type,'dAS'), damping_theta = SchwarzMethod{end}; end
RM_type = RoundingMethod{1}; RM_nmbdigits = RoundingMethod{2}; RM_Advanpix = RoundingMethod{3}; RM_CalcErrMtrx = RoundingMethod{4};
if isempty(plot)
    plot_convergence = 0;
else    
    plot_convergence = 1; x_mesh = plot{1}; u_plot = plot{2}; waittime = plot{3}; angle1 = plot{4}; angle2 = plot{5};
end

N = length(rhs); u_seq = NaN(SM_nmbiter,N);
SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo;
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end
CalcErrMtrxWithMoreDigs = true(SM_nmbsubdom,1);


%%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% do the partitioning
[ni,nj,~,~] = FindkPartition_Gander(A,SM_nmbsubdom_PwrOfTwo); rows_cumsum = cumsum(ni); cols_cumsum = cumsum(nj);
frst_rows = NaN(SM_nmbsubdom,1); last_rows = NaN(SM_nmbsubdom,1);  frst_cols = NaN(SM_nmbsubdom,1); last_cols = NaN(SM_nmbsubdom,1);
if strcmp(SM_type,'RAS'), mid_rows = NaN(SM_nmbsubdom,1); mid_cols = NaN(SM_nmbsubdom,1); end
for ind_subdom = 1:SM_nmbsubdom
    if ind_subdom == 1
        frst_rows(1) = 1; frst_cols(1) = 1;
    else
        frst_rows(ind_subdom) = rows_cumsum( (ind_subdom-1-1)*3 +1 ) + 1; frst_cols(ind_subdom) = cols_cumsum( (ind_subdom-1-1)*3 +1 ) + 1;
    end

    if strcmp(SM_type,'RAS') % if we do RAS -> get the midpoints of the subdomain overlaps
        if ind_subdom ~= SM_nmbsubdom
            mid_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +2 ); mid_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +2 );
        end
    end
 
    if ind_subdom ~= SM_nmbsubdom
        last_rows(ind_subdom) = rows_cumsum( (ind_subdom-1)*3 +3 ); last_cols(ind_subdom) = cols_cumsum( (ind_subdom-1)*3 +3 );
    else
        last_rows(ind_subdom) = N; last_cols(ind_subdom) = N;
    end
end
if debug == 2, plotPartition_Gander( A, ni, nj ); end % check for indices

if strcmp(SM_type,'RAS') % if we do RAS -> get the restriction indices
    frst_rows_RAS = NaN(SM_nmbsubdom,1); last_rows_RAS = NaN(SM_nmbsubdom,1);
    for ind_subdom = 1:SM_nmbsubdom
        if ind_subdom == 1
            frst_rows_RAS(ind_subdom) = 1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
        elseif ind_subdom == SM_nmbsubdom
            frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = N;
        else
            frst_rows_RAS(ind_subdom) = mid_rows(ind_subdom-1)+1; last_rows_RAS(ind_subdom) = mid_rows(ind_subdom);
        end
    end
end

%%% do the LowPrecision blocks
MySubdomFactors_LP = cell(SM_nmbsubdom,1); % need to keep the L,U,P,Q as well as other things (eg, A_diag for Stieltjess-type rounding)

for ind_subdom = 1:SM_nmbsubdom
    SubdomPrblmFacts = SubdomProbRounding( A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)), RM_nmbdigits, RM_type, RM_Advanpix, RM_CalcErrMtrx);
    MySubdomFactors_LP{ind_subdom} = SubdomPrblmFacts;
end






% %%% run Schwarz method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_dd = u_init; res_prev = rhs-A*u_dd; 


if debug ~= 0 % track convergence
    u_exact = A\rhs; u_err = NaN(SM_nmbiter+1,N); 
    err_norm = NaN(SM_nmbiter+1,1); u_err(1,:) = u_exact - u_init; err_norm(1) = norm(u_err(1,:));
end

if plot_convergence % plot convergence
    %%% plot the approximate solution
    nmb_int_gridcols = length(x_mesh)-2; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_exact,nmb_int_gridcols,nmb_int_gridcols)';
    mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Exact Solution'); pause(2);

    %%% plot the error
    nmb_int_gridcols = length(x_mesh)-2; u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_err(1,:),nmb_int_gridcols,nmb_int_gridcols)';
    mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Error'); pause(.1);
end


for ind_iter = 1:SM_nmbiter
    for ind_subdom = 1:SM_nmbsubdom
    
        %%% get my subdomain solutions  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curr_rhs_bfrperm = res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));

        if strcmp(RM_type,'Mmtrx')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); 
            if RM_Advanpix
                mp.Digits(RM_nmbdigits); curr_rhs_mp = mp(curr_rhs);
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); curr_u  = double(sol_perm(curr_Qvec_inv));
            else
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u  = double(sol_perm(curr_Qvec_inv));
            end
            % %%% check
            % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom));
            % curr_Qvec(curr_Qvec_inv) = 1:length(curr_Qvec_inv); %q_inv(q) = 1:length(q);
            % norm(  Ai(curr_Pvec,curr_Qvec) - curr_Lfact*curr_Ufact,'fro') / norm(Ai,'fro')
            % curr_u_check = Ai \ curr_rhs_bfrperm; norm(curr_u-curr_u_check)/norm(curr_u_check) * 1/condest(Ai)
        elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm_scaled = 1 ./ sqrt(Ai_diag) .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm_scaled(curr_Pvec); 
            if RM_Advanpix
                warning('off','all')
                mp.Digits(RM_nmbdigits); curr_rhs_mp = mp(curr_rhs); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); 
                warning('on','all')
            else
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs);
            end
             curr_u_unscaled = double(sol_perm(curr_Qvec_inv)); curr_u = 1 ./ sqrt(Ai_diag) .* curr_u_unscaled;
            %%% check
            % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom));
            % curr_Qvec(curr_Qvec_inv) = 1:length(curr_Qvec_inv); %q_inv(q) = 1:length(q);
            % Ai_sandwichScaled = diag(1./sqrt(Ai_diag)) * Ai * diag(1./sqrt(Ai_diag));
            % norm(  Ai_sandwichScaled(curr_Pvec,curr_Qvec) - curr_Lfact*curr_Ufact,'fro')
            % curr_u_check = Ai \ curr_rhs_bfrperm; norm(curr_u-curr_u_check)
        elseif strcmp(RM_type,'Stieltjess_RoundBi_Facts')
            ... % need to fill in
        elseif strcmp(RM_type,'Stieltjess_RoundBi_NoFacts') % use GMRES for subdomain solves
            ... % need to fill in
        end
        % % check
        % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); curr_u_check = Ai \ res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));


        %%% some checks to align with theory
        if ind_iter == 1

            if RM_CalcErrMtrx && strcmp(RM_type,'Mmtrx')
                Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{5};
                % Ai_inv_ErrMtrx = Ai\curr_ErrMtrx; check_Ai_inv_ErrMtrx_normest = pnorm(Ai_inv_ErrMtrx,2,1e-4,false);
                Ai_inv_ErrMtrx_Matvec = @(v,tflag) MatVecFun(v,tflag,Ai,curr_ErrMtrx); 
                Ai_inv_ErrMtrx_normest(ind_subdom) = svds(Ai_inv_ErrMtrx_Matvec,size(Ai),1);
                disp( append(append(append(append(append(append('For ',num2str(RM_nmbdigits)),' digits, we get ||inv(A'),num2str(ind_subdom)),')*Ei'),num2str(ind_subdom)),'||_2 ~ ',num2str(Ai_inv_ErrMtrx_normest(ind_subdom),'%.2e')) )
                if Ai_inv_ErrMtrx_normest(ind_subdom) < 1, CalcErrMtrxWithMoreDigs(ind_subdom) = false; end
                if RM_nmbdigits == 1 && strcmp(SM_type,'dAS') % i wanted to check subdomain blck conditioning
                    CondEst_Ai(ind_subdom) = condest(Ai);
                    disp( append(append(append(append(append(append('--> For ',num2str(size(Ai,1))),' subdom size, we get cond(A'),num2str(ind_subdom)),')'),' ~ '),num2str(CondEst_Ai(ind_subdom),'%.2e')) )
                else
                    CondEst_Ai(ind_subdom) = nan;
                end
            
            elseif RM_CalcErrMtrx && strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
                curr_ErrMtrx = MySubdomFactors_LP{ind_subdom}{6};
                Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); Ai_diag = spdiags(A,0); Ai_DiagRplcZeros = spdiags(zeros(size(Ai,1),1),0,Ai); 
                Ai_DiagSqrtInv = spdiags(1./sqrt(Ai_diag),0,size(Ai,1),size(Ai,2)); Ci = Ai_DiagSqrtInv * Ai_DiagRplcZeros * Ai_DiagSqrtInv;
                I_p_Ci = spdiags(ones(size(Ci,1),1),0,Ci);
                I_p_Ci_inv_ErrMtrx_Matvec = @(v,tflag) MatVecFun(v,tflag,I_p_Ci,curr_ErrMtrx); 
                I_p_Ci_inv_ErrMtrx_normest(ind_subdom) = svds(I_p_Ci_inv_ErrMtrx_Matvec,size(Ai),1);
                disp( append(append(append(append(append(append('For ',num2str(RM_nmbdigits)),' digits, we get ||inv(A'),num2str(ind_subdom)),')*Fi'),num2str(ind_subdom)),'||_2 ~ ',num2str(I_p_Ci_inv_ErrMtrx_normest(ind_subdom),'%.2e')) )
                %%%
                CondEst_Ai(ind_subdom) = condest(Ai);
                disp( append(append(append(append(append(append('For ',num2str(RM_nmbdigits)),' digits, we get cond(A'),num2str(ind_subdom)),')*u_s'),num2str(ind_subdom)),' ~ ',num2str(CondEst_Ai(ind_subdom)*10^(-RM_nmbdigits),'%.2e')) )
                if I_p_Ci_inv_ErrMtrx_normest(ind_subdom) < 1 &&  CondEst_Ai(ind_subdom)*10^(-RM_nmbdigits) < 1
                    CalcErrMtrxWithMoreDigs(ind_subdom) = false; end
            end
        end


        %%% the subdomain updating strategies, depending on Schwarz method type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ind_subdom == 1 %%% update my global solution for subdomain 1
            if strcmp(SM_type,'RAS') 
                u_dd = u_dd + [curr_u(1:last_rows_RAS(ind_subdom)); zeros(N-last_rows_RAS(ind_subdom),1)];
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [curr_u; zeros(N-last_rows(ind_subdom),1)];
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [curr_u; zeros(N-last_rows(ind_subdom),1)];
            end

        elseif ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
            if strcmp(SM_type,'RAS') 
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1;
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:end)]; 
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u]; 
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u]; 
            end

        else %%% update my global solution for middle subdomains
            if strcmp(SM_type,'RAS') 
                frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_in_left_overlap = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                last_row_inside_subdom = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_to_right_overlap = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
                u_dd = u_dd + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)];
            elseif strcmp(SM_type,'dAS')
                u_dd = u_dd + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u; zeros(N-last_rows(ind_subdom),1)];
            elseif strcmp(SM_type,'MS')
                u_dd = u_dd + [zeros(frst_rows(ind_subdom)-1,1); curr_u; zeros(N-last_rows(ind_subdom),1)];
            end
        end



        if strcmp(SM_type,'MS') % if multiplicative Schwarz -> update the residual after each subdomain solve
            res_prev = rhs-A*u_dd; 
        end
    end
    
    if strcmp(SM_type,'RAS') || strcmp(SM_type,'dAS') % if (restrcited) additive Schwarz -> update the residual after all subdomain solve
        res_prev = rhs-A*u_dd; 
    end
    u_seq(ind_iter,:) = u_dd; res_norm = norm(res_prev); 
    if res_norm/norm(rhs) < SM_relresacc, break; end

    if debug == 1 % track convergence 
        u_err(ind_iter+1,:) = u_exact - u_dd; err_norm(ind_iter+1) = norm(u_err(ind_iter+1,:)); 
        % if err_norm(ind_iter+1) < 1e-15 && err_norm(ind_iter) < 1e-15
        %     break
        % end
    end
    if plot_convergence % plot convergence
        % %%% plot the approximate solution
        % u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_dd,nmb_int_gridcols,nmb_int_gridcols)'; mesh(x_mesh,x_mesh,u_plot); view(angle1,angle2); pause(waittime);
        %%% plot the error
        u_plot(2:nmb_int_gridcols+1,2:nmb_int_gridcols+1)=reshape(u_err(ind_iter+1,:),nmb_int_gridcols,nmb_int_gridcols)'; mesh(x_mesh,x_mesh,u_plot); xlabel('x');ylabel('y'); view(angle1,angle2); title('Error'); pause(waittime);
    end
end


MyOutput = {u_seq};
if debug == 1 % track convergence 
    if plot_convergence, semilogy([0,1:SM_nmbiter], err_norm/err_norm(1), 'ro-'); pause(1); end
    if RM_CalcErrMtrx, if strcmp(RM_type,'Mmtrx'), InvErrMtrx_normest =  Ai_inv_ErrMtrx_normest; elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts'), InvErrMtrx_normest = I_p_Ci_inv_ErrMtrx_normest; end; else, InvErrMtrx_normest = nan; CondEst_Ai = nan; end
    MyOutput = {u_seq,u_exact,err_norm,CalcErrMtrxWithMoreDigs,InvErrMtrx_normest,CondEst_Ai};
end

end
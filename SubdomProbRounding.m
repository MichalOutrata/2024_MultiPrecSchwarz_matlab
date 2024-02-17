function[output] = SubdomProbRounding(A,d_s,TypeOfRounding,Advanpix,CalculateErrMtrx)
    if Advanpix
        mp.Digits(d_s);
        if strcmp(TypeOfRounding,'Mmtrx')
            output = MtrxFactorsLP_Mmtrx(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundBi_NoFacts')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('old version, not kept updated - wait for button press'); waitforbuttonpress;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            output = MtrcsLP_Stieltjess_DiagZerosUnscaled(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundBi_Facts')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('old version, not kept updated - wait for button press'); waitforbuttonpress;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            output = MtrxFactorsLP_Stieltjess(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundSandwichScaledBi_Facts')
            output = MtrxFactorsLP_Stieltjess_SandwichScaling(A,d_s,CalculateErrMtrx);
        end
    else
        if strcmp(TypeOfRounding,'Mmtrx')
            output = MtrxFactorsLP_Mmtrx_NoAdvanpix(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundBi_NoFacts')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('old version, not kept updated - wait for button press'); waitforbuttonpress;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            output = MtrcsLP_Stieltjess_DiagZerosUnscaled_NoAdvanpix(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundBi_Facts')
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp('old version, not kept updated - wait for button press'); waitforbuttonpress;
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            output = MtrxFactorsLP_Stieltjess_NoAdvanpix(A,d_s,CalculateErrMtrx);
        elseif strcmp(TypeOfRounding,'Stieltjess_RoundSandwichScaledBi_Facts')
            output = MtrxFactorsLP_Stieltjess_SandwichScaling_NoAdvanpix(A,d_s,CalculateErrMtrx);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[MyOutput] = MtrxFactorsLP_Mmtrx(A,d_s,ReturnErrMtrx)
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
    
    %%% old & not functional with advanpix
    % round positive entries
    pwrs_pos = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; %pwrs_pos = max(pwrs);
    v_PosRounded = ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos))./10.^(d_s-pwrs_pos);
    % round negative entries
    pwrs_neg = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; %pwrs_neg = max(pwrs);
    v_NegRounded = ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    v_Rounded = [v_NegRounded;v_PosRounded]; A_rounded_aux = sparse(i,j,v_Rounded(perm_inv)); A_rounded = mp(A_rounded_aux);
    
    %%% new & hopefully better
    % round positive entries
    pwrs_pos = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; %pwrs_pos = max(pwrs);
    v_mpPos = mp( ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos)) )./10.^(d_s-pwrs_pos);
    % round negative entries
    pwrs_neg = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; %pwrs_neg = max(pwrs);
    v_mpNeg = mp( ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg)) )./10.^(d_s-pwrs_neg);
    v_mp = [v_mpNeg;v_mpPos]; A_mp = sparse(i,j,v_mp(perm_inv));

    [L,U,p,q] = lu(A_mp,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    DoCheck = false;
    if ReturnErrMtrx
       if ~DoCheck, ErrMtrx_mp = double(A_mp - A); else
       ErrMtrx = A_rounded_aux - A; ErrMtrx_roundedmp = double(A_rounded - A);  ErrMtrx_mp = double(A_mp - A);  
       disp(append('norm of Ei < 0: ',num2str(norm(ErrMtrx(ErrMtrx<=0),'fro')))); 
       disp(append('norm of Ei_mp_oldway < 0: ',num2str(norm(ErrMtrx_roundedmp(ErrMtrx_roundedmp<=0),'fro')))); 
       disp(append('norm of Ei_mp_newway < 0: ',num2str(norm(ErrMtrx_mp(ErrMtrx_mp<=0),'fro')))); end
    else, ErrMtrx_mp=nan; end

    MyOutput = {L,U,p,q_inv,ErrMtrx_mp};
end

function[MyOutput] = MtrxFactorsLP_Stieltjess(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A);
    % bcs A was Stieltjess, we know that the diagonal is nonzero
    
    [i,j,v] = find(A_DiagRplcZeros);
    
    % round off-diag entries (they're all negative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    
    A_rounded = sparse(i,j,v_Rounded);  A_OnesDiag_BfrAdvanpix = spdiags(ones(size(A,1),1),0, diag(A_diag)\A_rounded);
    A_OnesDiag_Rounded = mp(A_OnesDiag_BfrAdvanpix);
    [L,U,p,q] = lu(A_OnesDiag_Rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    if ReturnErrMtrx, ErrMtrx = A_rounded - A_DiagRplcZeros; else, ErrMtrx=nan; end
    MyOutput = {L,U,p,q_inv,A_diag,ErrMtrx};
end

function[MyOutput] = MtrcsLP_Stieltjess_DiagZerosUnscaled(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A);
    % bcs A was Stieltjess, we know that the diagonal is nonzero
    
    [i,j,v] = find(A_DiagRplcZeros); 
    
    % round off-diag entries (they're all engative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);    
    A_rounded_ZerosDiag = mp( sparse(i,j,v_Rounded) ); 

    if ReturnErrMtrx, ErrMtrx = A_rounded_ZerosDiag - A_DiagRplcZeros; else, ErrMtrx=nan; end
    MyOutput = {A_rounded_ZerosDiag,A_diag,ErrMtrx};
end

function[MyOutput] = MtrxFactorsLP_Stieltjess_SandwichScaling(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A); A_DiagSqrtInv = spdiags(1./sqrt(A_diag),0,size(A,1),size(A,2));
    C_exact = A_DiagSqrtInv * A_DiagRplcZeros * A_DiagSqrtInv;
    
    [i,j,v] = find(C_exact);
    
    % round off-diag entries (they're all negative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_mp = mp( ceil(v.*10.^(d_s-pwrs_neg)) )./10.^(d_s-pwrs_neg);

    C_mp = sparse(i,j,v_mp); 
    I_p_C_mp = mp( spdiags(ones(size(A,1),1),0,C_mp) );
    [L,U,p,q] = lu(I_p_C_mp,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    DoCheck = false;
    if ReturnErrMtrx
        if ~DoCheck, ErrMtrx = double(C_mp - C_exact); else
        ErrMtrx = double(C_mp - C_exact); disp(append('norm of Fi < 0: ',num2str(norm(ErrMtrx(ErrMtrx<0),'fro')))); end
    else, ErrMtrx=nan; end

    MyOutput = {L,U,p,q_inv,A_diag,ErrMtrx};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[MyOutput] = MtrxFactorsLP_Mmtrx_NoAdvanpix(A,d_s,ReturnErrMtrx)
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
    pwrs_pos = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; %pwrs_pos = max(pwrs);
    v_PosRounded = ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos))./10.^(d_s-pwrs_pos);
    
    % round negative entries
    pwrs_neg = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; %pwrs_neg = max(pwrs);
    v_NegRounded = ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    
    v_Rounded = [v_NegRounded;v_PosRounded]; A_rounded = sparse(i,j,v_Rounded(perm_inv));
    [L,U,p,q] = lu(A_rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse
    %%% check
    % err = A_rounded - A; disp( numel(err(err<0)) );
    
    % %%% check 
    if ReturnErrMtrx, ErrMtrx = A_rounded - A; else, ErrMtrx=nan; end
    % err = A - A_rounded; if ~isempty(err(err<0)), disp('error'); end

    MyOutput = {L,U,p,q_inv,ErrMtrx};
end


function[MyOutput] = MtrxFactorsLP_Stieltjess_NoAdvanpix(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A);
    % bcs A was Stieltjess, we know that the diagonal is nonzero
    
    [i,j,v] = find(A_DiagRplcZeros);
    
    % round off-diag entries (they're all negative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
    
    A_rounded = sparse(i,j,v_Rounded);  A_OnesDiag_Rounded = spdiags(ones(size(A,1),1),0, diag(A_diag)\A_rounded);
    [L,U,p,q] = lu(A_OnesDiag_Rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    if ReturnErrMtrx, ErrMtrx = A_rounded - A_DiagRplcZeros; else, ErrMtrx=nan; end
    MyOutput = {L,U,p,q_inv,A_diag,ErrMtrx};
end

function[MyOutput] = MtrcsLP_Stieltjess_DiagZerosUnscaled_NoAdvanpix(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A);
    % bcs A was Stieltjess, we know that the diagonal is nonzero
    
    [i,j,v] = find(A_DiagRplcZeros); 
    
    % round off-diag entries (they're all engative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);    
    A_rounded_ZerosDiag = sparse(i,j,v_Rounded); 

    if ReturnErrMtrx, ErrMtrx = A_rounded_ZerosDiag - A_DiagRplcZeros; else, ErrMtrx=nan; end
    MyOutput = {A_rounded_ZerosDiag,A_diag,ErrMtrx};
end

function[MyOutput] = MtrxFactorsLP_Stieltjess_SandwichScaling_NoAdvanpix(A,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end
    A_diag = spdiags(A,0); A_DiagRplcZeros = spdiags(zeros(size(A,1),1),0,A); A_DiagSqrtInv = spdiags(1./sqrt(A_diag),0,size(A,1),size(A,2));
    C_exact = A_DiagSqrtInv * A_DiagRplcZeros * A_DiagSqrtInv;
    
    [i,j,v] = find(C_exact);
    
    % round off-diag entries (they're all negative)
    pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
    v_Rounded = ceil(v.*10.^(d_s-pwrs_neg)) ./ 10.^(d_s-pwrs_neg);

    C_Rounded = sparse(i,j,v_Rounded); I_p_C_Rounded = spdiags(ones(size(A,1),1),0,C_Rounded);
    [L,U,p,q] = lu(I_p_C_Rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse
    %%% check
    % norm(diag(1./sqrt(diag(A))) * A * diag(1./sqrt(diag(A))) -  B_SandwichScaled_Rounded - eye(size(A)), 'fro')
    % norm(diag(1./sqrt(diag(A))) * A * diag(1./sqrt(diag(A))) -  C_Rounded, 'fro')
    % norm(C_Rounded(p,q) - L*U, 'fro')

    % %%% check 
    DoCheck = false;
    if ReturnErrMtrx
        if ~DoCheck, ErrMtrx = C_Rounded - C_exact; else
        ErrMtrx = C_Rounded - C_exact; disp(append('norm of Fi < 0: ',num2str(norm(ErrMtrx(ErrMtrx<0),'fro')))); end
    else, ErrMtrx=nan; end

    MyOutput = {L,U,p,q_inv,A_diag,ErrMtrx};
end
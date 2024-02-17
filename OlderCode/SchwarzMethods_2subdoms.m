n=12; A=Laplacian(2,n); 

N_int_gridlength = fix(sqrt(size(A,1))); h=1/(N_int_gridlength+1);
f=zeros((n-2)^2,1); f(1:n-2:end)=-1; f = 1/h^2*f;

[~] = Schwarz_AS_RAS_2subdoms(A,f);

function [output] = Schwarz_AS_RAS_2subdoms(A,rhs)

[row_indices,col_indices,dim_A,ind_MidInA1] = Find2Partition_Gander(A); ind_MidInA2 = ind_MidInA1 + 1;

last_row_A1 = sum(row_indices(1:3)); last_col_A1 = sum(col_indices(1:3));
frst_row_A2 = row_indices(1)+1; frst_col_A2 = col_indices(1)+1;

A1 = A(1:last_row_A1,1:last_col_A1); A2 = A(frst_row_A2:end,frst_col_A2:end);

%%% check AS
R1=[speye(last_row_A1) , sparse(last_row_A1, dim_A-last_row_A1 )];
R2=[sparse(dim_A-frst_row_A2+1, frst_row_A2-1 ) , speye(dim_A-frst_row_A2+1) ];
A1_check=R1*A*R1'; err_1 = A1_check - A1; check_1 = isempty( err_1(err_1~=0) );
A2_check=R2*A*R2'; err_2 = A2_check - A2; check_1 = isempty( err_2(err_2~=0) );

%%% check RAS
R1_tilde = R1; R1_tilde(ind_MidInA2:end,:) = 0;
R2_tilde = R2; R2_tilde(1:ind_MidInA2-frst_row_A2,:) = 0;



%%% calculate exact soltuion
u = A\rhs;

u_dd=zeros(size(rhs)); 

%%% plotting
N_int_gridlength = fix(sqrt(size(A,1)));
u_plot=zeros(N_int_gridlength+2); u_plot(1:N_int_gridlength+2,1)=1;



for i=1:10
    %%% plotting 
    u_plot(2:N_int_gridlength+1,2:N_int_gridlength+1)=reshape(u_dd,N_int_gridlength,N_int_gridlength)';
    h=1/(N_int_gridlength+1); x=0:h:1;
    mesh(x,x,u_plot)
    xlabel('x');ylabel('y');
    view(52.5,30); pause(2);

    %%% calculating eror
    err(i) = max(max(abs(u-u_dd)));

    %%% run RAS
    res_prev = rhs-A*u_dd; 
    u_A1sol = A1\res_prev(1:last_row_A1); u_A2sol = A2\res_prev(frst_row_A2:end);
    u_dd = u_dd + [u_A1sol(1:ind_MidInA1); zeros(dim_A-ind_MidInA1,1)] + [zeros(dim_A-ind_MidInA1,1) ; u_A2sol(1+ ind_MidInA2-frst_row_A2:end)];
    u_dd_check = u_dd + R1_tilde'*(A1\(R1*res_prev)) + R2_tilde'*(A2\(R2*res_prev));

    % %%% run AS
    % res_prev = rhs-A*u_dd;
    % u_dd = u_dd + [A1\res_prev(1:last_row_A1); zeros(dim_A-last_row_A1,1)] + [zeros(dim_A-last_row_A1,1) ; A2\res_prev(frst_row_A2:end)];
    % u_dd_check = u_dd + R1'*(A1\(R1*res_prev)) + R2'*(A2\(R2*res_prev));

    % %%% run MS
    % res_prev = rhs-A*u_dd;
    % u_dd_half = u_dd + [A1\res_prev(1:last_row_A1); zeros(dim_A-last_row_A1,1)]; res_half = rhs-A*u_dd;
    % u_dd = u_dd_half + [zeros(dim_A-last_row_A1,1) ; A2\res_half(frst_row_A2:end)];
    % u_dd_check_half = u_dd + R1'*(A1\(R1*res_prev)); res_check_half = rhs-A*u_dd_check_half;
    % u_dd_check = u_dd + R2'*(A2\(R2*res_check_half));
end


output = 0;

end
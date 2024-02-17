clear; clc; close('all'); format longG;

% % %%% test 1
% % R = sprand(3,3,.7)*1e+2; disp(full(R));
% % d_s = 2; mp.Digits(d_s);
% % R_round = mp(R); disp(full(R_round));
% % err = R - R_round; disp(full(err));
% 
% %%% test 2
% % v_sorted = [-0.00178132454400338;-0.00128014399720173; -0.00999080394761361; 0.00171121066356432; 0.00032600820530528; 0.0056119979270966]; ind_AbsValMin = 3; disp(v_sorted);
% v_sorted = [-rand(3,1);rand(3,1)]*10^(0); ind_AbsValMin = 3; disp(v_sorted);
% d_s = 3; mp.Digits(d_s);
% pwrs_pos = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; pwrs_pos = max(pwrs_pos);
% v_PosRounded = ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos))./10.^(d_s-pwrs_pos);
% disp(mp(v_PosRounded))
% pwrs_neg = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; pwrs_neg = max(pwrs_neg);
% v_NegRounded = ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg))./10.^(d_s-pwrs_neg);
% disp(mp(v_NegRounded))
% 
% i = [1;1;2;2;3;3]; j = [1;2;1;2;2;3];
% v_Rounded = [v_NegRounded;v_PosRounded]; A_rounded_aux = sparse(i,j,v_Rounded); 
% %%% check
% A = sparse(i,j,v_sorted); err = A_rounded_aux - A; disp(full(err))
% 
% A_rounded = mp(A_rounded_aux); [L,U,p,q] = lu(A_rounded,'vector'); q_inv(q) = 1:length(q);
% %%% check
% b = A*([1;2;3]);
% rhs_perm = b(p); u_preperm = U\(L\rhs_perm); u_out = u_preperm(q_inv); disp(u_out);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mp.Digits(1);
% x = -0.505625850340136 *1e+3;
% 
% x = -505; disp( append('the number "x": ', num2str(x)) )
% curr_nmb = mp(x); disp( append('the number "x" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(x)-x; disp( append('the number "mp(x) - x" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(x)-x,1); disp( append('the number "round(mp(x) - x,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -500-x; disp( append('the number "MyRoundin(x) - x": ', num2str(curr_nmb) ))
% 
% y = -514.242630385488;
% y = -514; disp( append('the number "y": ', num2str(x)) )
% curr_nmb = mp(y); disp( append('the number "y" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(y)-y; disp( append('the number "mp(y) - y" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(y)-y,1); disp( append('the number "round(mp(y) - y,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -500-y; disp( append('the number "MyRoundin(y) - y": ', num2str(curr_nmb) ))
% 
% %%%
% mp.Digits(1);
% 
% x = -5.05; disp( append('the number "x": ', num2str(x)) )
% curr_nmb = mp(x); disp( append('the number "x" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(x)-x; disp( append('the number "mp(x) - x" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(x)-x,1); disp( append('the number "round(mp(x) - x,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -0.500-x; disp( append('the number "MyRoundin(x) - x": ', num2str(curr_nmb) ))
% 
% y = -5.14; disp( append('the number "y": ', num2str(x)) )
% curr_nmb = mp(y); disp( append('the number "y" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(y)-y; disp( append('the number "mp(y) - y" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(y)-y,1); disp( append('the number "round(mp(y) - y,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -0.500-y; disp( append('the number "MyRoundin(y) - y": ', num2str(curr_nmb) ))
% 
% %%%
% mp.Digits(1);
% 
% x = -0.000505; disp( append('the number "x": ', num2str(x)) )
% curr_nmb = mp(x); disp( append('the number "x" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(x)-x; disp( append('the number "mp(x) - x" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(x)-x,1); disp( append('the number "round(mp(x) - x,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -0.000500-x; disp( append('the number "MyRoundin(x) - x": ', num2str(curr_nmb) ))
% 
% y = -0.000514; disp( append('the number "y": ', num2str(x)) )
% curr_nmb = mp(y); disp( append('the number "y" in mp: ', num2str(curr_nmb)) )
% curr_nmb =  mp(y)-y; disp( append('the number "mp(y) - y" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  round(mp(y)-y,1); disp( append('the number "round(mp(y) - y,1)" in mp: ', num2str(curr_nmb) ))
% curr_nmb =  -0.000500-y; disp( append('the number "MyRoundin(y) - y": ', num2str(curr_nmb) ))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mp.Digits(16);
x = -rand(20,1)*1e+3;
x_mp16 = mp(x);
pwr = floor( log10(-x_mp16) ); 
mp.Digits(1);
x_mp1 = mp(x);
x_mp1_from_mp16 = mp(x_mp16);
x_mp1_rescaled = mp( ceil( x_mp16 .* 10.^(-pwr) ) );
x_mp1_rescaled_indbl = mp( ceil( x .* 10.^(-pwr) ) );
vec1 = double(x_mp1) - x; norm(vec1(vec1<0))
vec2 = double(x_mp1_from_mp16) - x; norm(vec2(vec2<0))
vec3 = double(x_mp1_rescaled).*10.^(pwr) - x; norm( vec3(vec3<0) )
vec4 = double(x_mp1_rescaled_indbl).*10.^(pwr) - x; norm( vec4(vec4<0) )

%%%

mp.Digits(16);
x = [-rand(20,1)*1e+3; -rand(20,1)*1e+1];
x_mp16 = mp(x);
pwr = floor( log10(-x_mp16) ); 
mp.Digits(1);
x_mp1 = mp(x);
x_mp1_from_mp16 = mp(x_mp16);
x_mp1_rescaled = mp( ceil( x_mp16 .* 10.^(-pwr) ) );
x_mp1_rescaled_indbl = mp( ceil( x .* 10.^(-pwr) ) );
vec1 = double(x_mp1) - x; norm(vec1(vec1<0))
vec2 = double(x_mp1_from_mp16) - x; norm(vec2(vec2<0))
vec3 = double(x_mp1_rescaled).*10.^(pwr) - x; norm( vec3(vec3<0) )
vec4 = double(x_mp1_rescaled_indbl).*10.^(pwr) - x; norm( vec4(vec4<0) )

%%%

mp.Digits(16);
x = [-rand(20,1)*1e+3; -rand(20,1)*1e+1; -rand(20,1)*1e-1];
x_mp16 = mp(x);
pwr = floor( log10(-x_mp16) ); 
mp.Digits(1);
x_mp1 = mp(x);
x_mp1_from_mp16 = mp(x_mp16);
x_mp1_rescaled = mp( ceil( x_mp16 .* 10.^(-pwr) ) );
x_mp1_rescaled_indbl = mp( ceil( x .* 10.^(-pwr) ) );
vec1 = double(x_mp1) - x; norm(vec1(vec1<0))
vec2 = double(x_mp1_from_mp16) - x; norm(vec2(vec2<0))
vec3 = double(x_mp1_rescaled.*10.^(pwr)) - x; norm( vec3(vec3<0) )
vec4 = double(x_mp1_rescaled_indbl.*10.^(pwr)) - x; norm( vec4(vec4<0) )


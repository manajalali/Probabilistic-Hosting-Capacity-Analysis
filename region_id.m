function[indic_degen,indic_self,P_coeff,p_bias,M_p,r_p,M_d,r_d]=region_id(theta,A,H,E,b,d,B,f,F,C,C_ineq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function region_id takes a solved QP from MPP_QP and identifies a
% critical region (if possible) and reports the polytopic description as:
%
%                     P_coef*theta <= p_bias.
% This function also calcluates the corresponding affine mapping coeficients
% for primal variables (M_p,r_p) and for dual variables (M_d,r_d).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initializing

P_coeff = 0;
p_bias = 0;

M_p = 0;
r_p = 0;
M_d = 0;
r_d = 0;

active_threshold = 1e-5; % The threshold to declare a constraint active
boundary_accuracy = 1e-4; % The accuracy for region boundaries

indic_degen = 0; % An indicator for degenerate region
indic_self = 1;  % An indicator to show if theta falls into its own region 
%% identifying active and inactive constraints
ind_all = (1:size(A,1))';
ind_active = find(abs(C_ineq)<= active_threshold); 
ind_non_slack = ind_all;
ind_non_slack(ind_active)=[];
ind_non = ind_non_slack;

%% Finding corresponding rows of A,E,b based on active/inactive constraints

A_tilde = A(ind_active,:);
A_bar = A(ind_non,:);
b_tilde = b(ind_active);
b_bar = b(ind_non);
E_tilde = E(ind_active,:);
E_bar = E(ind_non,:);

Y=[ ...
    H A_tilde' B';
    A_tilde zeros(length(ind_active),length(ind_active)) zeros(length(ind_active),size(B,1));
    B zeros(size(B,1),length(ind_active)) zeros(size(B,1),size(B,1))];

degen_indic = 0;

if(rank(Y) < size(Y,1)) %checking if the rgion is degenerate
    indic_degen = 1;
else
    M=inv(Y)*[-C;E_tilde;F];
    r=inv(Y)*[-d;b_tilde;f];
    
    M_p=M(1:size(A,2),:);
    M_d=M(size(A,2)+1:size(A,2)+size(ind_active),:);
    
    r_p=r(1:size(A,2),:);
    r_d=r(size(A,2)+1:size(A,2)+size(ind_active),:);
    

    P_coeff = [A_bar*M_p-E_bar;-M_d];
    p_bias  = [b_bar-A_bar*r_p;r_d];
      
    
    indic = P_coeff*theta - p_bias;
    indic = indic <= boundary_accuracy;
    
    if (sum(indic) == size(indic,1))
        indic_self = 1;
    else
        indic_self = 0;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The master code of running the MPP-PHCA algorithm of the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%% Loading the needed matrices and data to run the MPP-PHCA

%loading the needed matrix and data from the Pre_load() funtion.

%[A,H,Hi,E,b,N,d,B,f,F,C,pc,pg,qc] = Pre_load_123();
[A,H,Hi,E,b,N,d,B,f,F,C,pc,pg,qc] = Pre_load_8500();

scale = 1; %injection scaling
boundary_accuracy = 1e-4; % The accuracy for region boundaries

%%  Building Theta
Theta = Theta_maker(pg,pc,qc,scale);

%% Random shuffling of Theta (Step 3)
rand_ind = randperm(size(Theta,2));
Theta = Theta(:,rand_ind);
%%
x_sol = zeros(size(A,2),size(Theta,2)); %the OPF solutions will be recorded in x

flag_theta=zeros(1,size(Theta,2)); % Flag is the critical region that theta belongs to
                                   % Degeneratethetas are recorded in
                                   % region -1 and thetas with numerical
                                   % issues are recorded in region -2
j=1;
indic_degen = [];

solver_time = [];
yalmip_time = [];
check_time = [];
solve_time = [];

while (find(flag_theta==0))
    
    yalmip('clear')
    numel(find(flag_theta==0))
    ind_theta=find(flag_theta==0);
    theta=Theta(:,min(ind_theta));
    
    tic
    [x,res,C_ineq,yalmiptime,solvertime]=MPP_QP(theta,A,H,E,b,d,B,f,F,C);
    
    solve_time(end+1) = toc;
    solver_time(end+1) = solvertime;
    yalmip_time(end+1) = yalmiptime;
    
    tic
    [indic_degen,indic_self,P_coeff,p_bias,M_p,r_p,M_d,r_d] = ...
        region_id(theta,A,H,E,b,d,B,f,F,C,C_ineq);
    
    if(indic_degen==1)
        
        flag_theta(min(ind_theta))=-1;
        x_sol(:,min(ind_theta)) = x;
        check_time(end+1) = toc;
        
        continue
    end
    
    if(indic_self==0)
        
        flag_theta(min(ind_theta))=-2;
        x_sol(:,min(ind_theta)) = x;
        check_time(end+1) = toc;
        
        continue
    end
    
    poly.P_coeff = P_coeff;
    poly.p_bias = p_bias;
    
    for i=1:length(ind_theta)
        
        if (flag_theta(ind_theta(i))==0)
            
            indic = P_coeff*Theta(:,ind_theta(i)) - p_bias;
            indic = indic <= boundary_accuracy;
            
            if (sum(indic) == size(indic,1))
                x_sol(:,ind_theta(i))=M_p*Theta(:,ind_theta(i))+r_p;
                flag_theta(ind_theta(i))=j;
            end
        end
    end
    check_time(end+1) = toc;
    
    j=j+1;    
end
%% Resorting the solution to comply with the non-shuffled Theta

x_sol_resorted = zeros(size(x_sol));
flag_theta_resorted = zeros(size(flag_theta));

x_sol_resorted(:,rand_ind) = x_sol;
flag_theta_resorted(:,rand_ind) = flag_theta;




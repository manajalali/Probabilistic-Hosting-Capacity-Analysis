function [A,H,Hi,E,b,N,d,B,f,F,C,pc,pg,qc] = Pre_load_123()
% This function requires 2 data files:
% 1) power injections of the feeder
% 2) R, X, I (incidence matrix), N(number of nodes) , resistance and
% reactance of load tap changers (RLDC, XLDC), bus indeces of voltage
% regulators. In case the feeder is divided into several subgraphs, or
% several subgraphs are being used, all the aforementioned parameters must
% be stored in a struct. 
% In this file, the following parameters can be changed :
% (v_min, v_max) (minimum and maximum per unit volatge), vr( voltage regulator setpoints),
% delta(regulator bandwidth), beta (trade-off parameter between power loss
% and voltage regulator), (nu, gamma) penalty coefficients for 's'.

Matrix_generator_123(); % Needed to generate 'Feeder_matrices.mat';

load('power_injection_data.mat');
load('Feeder_matrices.mat');

%% Initial values needed:
beta=0.2;
gamma=1;
nu=20;
vr=[1;1;1;1];
delta=[2;2;1;2]./120;% later make this automatic
v_max=1.03;
v_min=0.97;

% Nt: total number of buses in the feeder.
IN=eye(Nt);

% S_dest, S_origin: selection matrices, selecting the destimation index and
% origin index of the voltage regulators.
S_dest=IN(VR_dest,:);
S_origin=IN(VR_origin,:);

% indeces of zero injection and non-zero injection buses.
ind_zero=find(sum(qg_max_new')==0);
K=[1:1:Nt];
ind_nonzero=setdiff(K,ind_zero);

% Building selection matrices for selecting the zero and non-zero injection
% buses.
S_zero=eye(Nt);
S_nonzero=eye(Nt);
S_nonzero(ind_zero,:)=[];
S_zero(ind_nonzero,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Loss Function %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  H=[ H11  H12 0
%      H12' H22 0
%      0    0   nu];
%  H corresponds to the quadratic term of the loss function: x'*H*x
%  [qg v0 s]'*H*[qg;v0;s]= qg'*H11*qg + 2*qg'*H12*v0+v0*H22*v0+nu*s^2

H11=beta*(S_subgraph'*X_tilda'*X_tilda*S_subgraph) + (1-beta)*(2*S_subgraph'*R_tilda*S_subgraph);

for i=1:subgraph_Number
    NN=N{i};
    X_times_1{i}=X{i}'*ones(NN,1);
end
H12=beta*S_subgraph'*blkdiag(X_times_1{:});

H22=beta*(blkdiag(N{:})+eye(subgraph_Number));

H=2*[H11 H12 zeros(Nt,1);H12' H22 zeros(subgraph_Number,1);zeros(1,Nt) zeros(1,subgraph_Number) nu];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% C= [ C11 C22 0;
%       0  0   0];
% C sorresponds to the linear term in the loss: x'*C*theta
% [qg' v0' s]*C*[p q_max q_min]=
% qg'*X'*R*p +p'*R'*v0

C11=2*beta*S_subgraph'*R_tilda'*X_tilda*S_subgraph;

for i=1:subgraph_Number
    R_times_1{i}=R{i}'*ones(N{i},1);
end
C22=2*beta*S_subgraph'*blkdiag(R_times_1{:});

C=[C11 C22 zeros(Nt,1)
    zeros(2*Nt,Nt+subgraph_Number+1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d corresponds to the linear term in the loss 
% d'*x = qg'*X'*1 + v0*1*1 + gamma*s


N_vector=[];
for i=1:subgraph_Number
   N_vector=[N_vector;N{i}]; 
end
d=[-2*beta*S_subgraph'*X_tilda'*ones(size(X_tilda,1),1); -2*beta*N_vector+ones(subgraph_Number,1);gamma];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Constraints  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints are divided into relaxable and non-relaxable constraints:
% A1*x <= E1*theta + b1
% B1*x == F1*theta + f1

% A2*x <= E2*theta + b2
% B2*x == F2*theta + f2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qg_min <= qg <= qg_max

A1_1=[S_nonzero zeros(length(ind_nonzero),subgraph_Number) zeros(length(ind_nonzero),1);
     -S_nonzero zeros(length(ind_nonzero),subgraph_Number) zeros(length(ind_nonzero),1)];

E1_1=[zeros(length(ind_nonzero),Nt) S_nonzero zeros(length(ind_nonzero),Nt);
    zeros(length(ind_nonzero),Nt) zeros(length(ind_nonzero),Nt) -S_nonzero];

b1_1=zeros(2*length(ind_nonzero),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s>0 
A1_2=[zeros(1,Nt+subgraph_Number) -1];

E1_2=[zeros(1,3*Nt)];

b1_2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q(VR_origin)=sum(q in the feeder connected to VR)

B_0=[S_dest(2:end,:) zeros(subgraph_Number-1,subgraph_Number) zeros(3,1)];

F_0=zeros(subgraph_Number-1,3*Nt);

f_0=zeros(subgraph_Number-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_n - r_LDC * P_mn - x_LDC * Q_mn == vr_mn 

delta_S= S_origin-S_dest;
delta_S(1,:)=S_dest(1,:);

B1_1=-blkdiag(XLDC{:})*delta_S;
B1_1=[B1_1 eye(subgraph_Number) zeros(subgraph_Number,1)];

F1_1=blkdiag(RLDC{:})*delta_S;
F1_1=[F1_1 zeros(subgraph_Number,2*Nt)];

f1_1=vr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q_substation = sum (q_feeder)

B1_2a=-S_origin(2:end,:);

for i=1:subgraph_Number
   ones_v0{i}=ones(1,N{i});
end
B1_2b=blkdiag(ones_v0{:});
B1_2b(1,:)=[];
B1_2b=B1_2b*S_subgraph;

B1_2=B1_2a+B1_2b;
B1_2=[B1_2 zeros(subgraph_Number-1,subgraph_Number) zeros(subgraph_Number-1,1)];

F1_2=zeros(subgraph_Number-1,3*Nt);
f1_2=zeros(subgraph_Number-1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting the injection of the zero-injection buses to zero

B1_3=[S_zero zeros(length(ind_zero),subgraph_Number+1)];

F1_3=[zeros(length(ind_zero),Nt) S_zero zeros(length(ind_zero),Nt)];

f1_3=zeros(length(ind_zero),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_n = vr

B1_4=[zeros(subgraph_Number,Nt) eye(subgraph_Number) zeros(subgraph_Number,1)];

F1_4=zeros(subgraph_Number,3*Nt);

f1_4=vr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatanating the non-relaxable constraints:

A1=[A1_1;A1_2];
E1=[E1_1;E1_2];
b1=[b1_1;b1_2];

B1=[B1_1;B1_2;B1_3];
F1=[F1_1;F1_2;F1_3];
f1=[f1_1;f1_2;f1_3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Relaxable Constraints %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (vr-delta)/1.1 <= vm <= (vr-delta)/0.9

Sm=S_origin*S_subgraph';

XS=[];
for i=1:subgraph_Number
   XS=[XS;X{i}*S{i}']; 
end
A2_1=[Sm*XS  eye(subgraph_Number) -ones(subgraph_Number,1);
     -Sm*XS -eye(subgraph_Number) -ones(subgraph_Number,1)];

RS=[];
for i=1:subgraph_Number
   RS=[RS;R{i}*S{i}']; 
end
E2_1=[-Sm*RS zeros(subgraph_Number,2*Nt);
       Sm*RS zeros(subgraph_Number,2*Nt)];

b2_1=[(vr+delta)/0.9;(delta-vr)/1.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  v_min <= v <= v_max

A_v0=blkdiag(ones_v0{:})';

A2_2=[XS A_v0 -1*ones(Nt,1)];


E2_2=[-RS zeros(Nt,2*Nt)];

b2_2=v_max*ones(Nt,1);

A2_3=[-XS -A_v0 -ones(Nt,1)];

E2_3=-E2_2;

b2_3=-v_min*ones(Nt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concataning the relaxable constraints:

A2=[A2_1;A2_2;A2_3];

E2=[E2_1;E2_2;E2_3];

b2=[b2_1;b2_2;b2_3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All constraints 

A=[A1;A2];
E=[E1;E2];
B=B1;
F=F1;
b=[b1;b2];
f=f1;
H = (H+H')/2;

C = C';

%% Normalizing (optional)
%
C = C/norm(H);
d = d/norm(H);
H = H/norm(H);

E = E/norm(A);
b = b/norm(A);
A = A/norm(A);


F = F/norm(B);
f = f/norm(B);
B = B/norm(B);

%%
Hi = inv(H);
N = Nt;
end

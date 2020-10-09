function [A,H,Hi,E,b,N,d,B,f,F,C,pc,pg,qc,scale] = Pre_load_8500()

% This function requires 2 data files:
% 1) power injections of the feeder
% 2) R, X, I (incidence matrix), N(number of nodes) , bus indeces of voltage
% regulators. In case the feeder is divided into several subgraphs, or
% several subgraphs are being used, all the aforementioned parameters must
% be stored in a struct. 
% In this file, the following parameters can be changed :
% (v_min, v_max) (minimum and maximum per unit volatge), vr( voltage regulator setpoints),
% delta(regulator bandwidth), beta (trade-off parameter between power loss
% and voltage regulator), (nu, gamma) penalty coefficients for 's'.

Matrix_generator_8500();

load data1160.mat % feeder matrices 
load Nodal_data_1160.mat %injection data

data = Nodal_data_1160;
pc = data.pc';
pg = data.pg';
qc = data.qc';


S_base = 1; %MVA


%% Data Pre-processing
V_base = data.base_KV; %KV
S_base = data.base_MVA;

Z_base=(V_base^2)/S_base;

R=R/Z_base;

X=X/Z_base;
%% Initial values needed:
td=0;
scale=1;
gamma=1;
nu=20;
v_max=1.03;
v_min=0.97;

N = size(X,1);
L = size(X,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Loss Function %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  H=[ H11  H12 0
%      H12' H22 0
%      0    0   nu];
%  x=[qg v0 s]
%  H corresponds to the quadratic term of the loss function: x'*H*x
%  [qg v0 s]'*H*[qg;v0;s]= qg'*H11*qg + 2*qg'*H12*v0+v0*H22*v0+nu*s^2

H11=td*(X'*X) + (1-td)*(R);
H12=td*[X'*ones(N,1)];
H22=1;
H=2*[H11 H12 zeros(N,1);H12' H22 0;zeros(1,N) 0 nu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C= [ C11 C22 0;
%       0  0   0];
% C sorresponds to the linear term in the loss: x'*C*theta
% [qg' v0' s]*C*[p q_max q_min]=
% qg'*X'*R*p +p'*R'*v0

C11=2*td*R'*X;
C22=2*td*R'*ones(N,1);
C=[C11 C22 zeros(N,1)
    zeros(2*N,N+2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d corresponds to the linear term in the loss 
% d'*x = qg'*X'*1 + v0*1*1 + gamma*s

d=[-2*td*X'*ones(size(X,1),1); -2*td*N;gamma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Constraints  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints are divided into relaxable and non-relaxable constraints:
% A1*x <= E1*theta + b1
% B1*x == F1*theta + f1

% A2*x <= E2*theta + b2
% B2*x == F2*theta + f2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reactive power generation constraints |q|<sqrt{s^2-pg^2}

scaling =1;

A1_1=[eye(N) zeros(N,1) zeros(N,1);
    -eye(N) zeros(N,1) zeros(N,1)]*scaling;

E1_1=[zeros(N,N) eye(N) zeros(N,N);
    zeros(N,N) zeros(N,N) -eye(N)]*scaling;

b1_1=zeros(2*N,1)*scaling;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s>0
A1_2=[zeros(1,N+1) -1];

E1_2=[zeros(1,3*N)];

b1_2=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatanating the non-relaxable constraints:

A1=[A1_1;A1_2];
E1=[E1_1;E1_2];
b1=[b1_1;b1_2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Relaxable Constraints %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_min<v<v_max

A_v0 = ones(N,1);

A2_2 = X;
A2_2 = [A2_2 A_v0 -1*ones(N,1)];

E2_2 = R;

E2_2 = [-E2_2 zeros(N,2*N)];

b2_2 = v_max*ones(N,1);

A2_3 = X;
A2_3 = [-A2_3 -A_v0 -ones(N,1)];

E2_3=-E2_2;

b2_3 = -v_min*ones(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concataning the relaxable constraints:
A2=[A2_2;A2_3];

E2=[E2_2;E2_3];

b2=[b2_2;b2_3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All constraints 
A=[A1;A2];
E=[E1;E2];

B=zeros(1,N+2);
B(N+1) = 1;
F=zeros(1,3*N);
b=[b1;b2];
f=1;

H = (H+H')/2;

C = C';

%% Normalizing (optional)
%
% C = C/norm(H);
% d = d/norm(H);
% H = H/norm(H);
% 
% E = E/norm(A);
% b = b/norm(A);
% A = A/norm(A);


% F = F/norm(B);
% f = f/norm(B);
% B = B/norm(B);

Hi = inv(H);

end

% Comments:
% For running this file, you will need a matlab 2018 or higher. 
% This code generates a .mat file including the feeder marices used in the master code and loads two files:
% 1) A data struct including the lines, loads, regulators, and breakers
% information. These information include the bus indeces, impedances, etc. This data struct should include the characteristics of a
% single-phase feeder. Here the name of the file including the data struct
% is: 'IEEE_123_data_structure.mat'
% 2) A struct including the configuration of all the lines:'line_configs.mat'
% 
%%%%%%%%%%% NOTE: %%%
% The base_MVA and voltage MVA should be specified:

clear

% loading the required data:
load('IEEE_123_data_structure.mat');
load('line_configs.mat');
base_MVA=1;
base_V=4.160;

% Converting the tables including the feeder data into Matlab arrays,
% otherwise matlab can not process the imported data.
lines=table2array(data.lines);
OLTC=table2array(data.OLTC);
VR_origin=table2array(data.OLTC(:,1));
VR_dest=table2array(data.OLTC(:,2));
breaker=table2array(data.breakers(:,1:2));
buses=table2array(data.buses);

N_buses=[1:1:length(buses)];
buses=[N_buses' buses];
ind_new_buses=find(buses(:,1)~=buses(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subgraph_Number: number of on-line tap changers
subgraph_Number=size(OLTC,1);

% adding the breakers:
for k=1:length(breaker)
    ind_line2=find(lines(:,2)==breaker(k,1));
    ind_line1=find(lines(:,1)==breaker(k,1));
    lines(ind_line2,2)=breaker(k,2);
    lines(ind_line1,1)=breaker(k,2);
end
L=length(lines);

ind_a0=find(lines(:,2)==1);
l0=lines(ind_a0,:);
lines(ind_a0,:)=[];
lines=[l0;lines];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The original numbering of the IEEE feeder includes numbers such as 150
% and 149, the graph-theory matlab commands used here expects to see
% consecutive numbers for buses. Otherwise, graph command assumes
% individual disconnected nodes for the graph.
%% renaming the lines:
ind_x1=find(lines(:,1)>=ind_new_buses(1));
ind_x2=find(lines(:,2)>=ind_new_buses(1));
ind_x3=find(VR_dest>=ind_new_buses(1));
ind_x4=find(VR_origin>=ind_new_buses(1));


high_bus=unique([unique(lines(ind_x1,1));unique(lines(ind_x2,2));VR_dest(ind_x3);VR_origin(ind_x4)]);
lines1=lines(:,1);
lines2=lines(:,2);
for i=1:length(high_bus)
    lines1(lines1==high_bus(i))=ind_new_buses(1)+i;
    lines2(lines2==high_bus(i))=ind_new_buses(1)+i;
    VR_dest(VR_dest==high_bus(i))=ind_new_buses(1)+i;
    VR_origin(VR_origin==high_bus(i))=ind_new_buses(1)+i;
end

%% Renaming the buses

coincidenceMatrix = [lines1 lines2];
nodeIDs = unique(coincidenceMatrix);
numberofNodes = length(nodeIDs);
nodeIDs(: , 2) = [1:numberofNodes]';

ind_VR_dest_old=[];
for i = 1:numberofNodes
    coincidenceMatrix(coincidenceMatrix == nodeIDs(i)) = i;
    ind_VR_dest_old=[ind_VR_dest_old (VR_dest==nodeIDs(i))*i];
%     ind_VR_dest_old(i)=find(VR_dest==nodeIDs(i))
    VR_dest(VR_dest==nodeIDs(i))=i;
    VR_origin(VR_origin==nodeIDs(i))=i;
end

lines(:,1:2)=coincidenceMatrix;

%% Dividing into subgraphs: 
% In this part, the system is divided into subgraphs so that the downstream
% from of each voltage regulator is studied separately. 

% G: The complete graph of the feeder
G=graph(lines(:,1),lines(:,2));

% Removing the edges of the feeder corresponding to voltage regulators
G_divided=rmedge(G,VR_origin(2:end),VR_dest(2:end));


nodenames=[1:1:L+1];

% Incidence matrix of the original feeder
I_G=full(G.incidence);
I_G=I_G';

% Finding the connected components of the feeder after division into
% subgraphs:
[beans,beansizes]=conncomp(G_divided);

% IN=identity matrix of size N
IN=eye(L+1);

% Finding the bus-indeces, graph, selection matrix, and incidence matrix of subgraphs 
for i=1:length(beansizes)
    ind_I{i}= unique([find(nodenames.*(beans==i))],'stable');
    G_sub=subgraph(G,ind_I{i});
    S{i}=IN(:,ind_I{i});
    I{i}=full(G_sub.incidence)';
    I{i}=I{i}(:,2:end);
end

%% Building R and X matrices (single phase):

type=zeros(L,1);
for l=1:L
    type(lines(l,2))=lines(l,4);
    Z_3p{lines(l,2)}=config{type(lines(l,2))}*lines(l,3)*(1/5280);
    Z_1p(l)=sum(diag(Z_3p{lines(l,2)}))/(nnz(diag(Z_3p{lines(l,2)})));
end

r=real(Z_1p);
x=imag(Z_1p);

% Building the structs including R,X,N of the subgraphs:
for i=1:length(beansizes)
    R_inv=I_G(:,ind_I{i}(2:end))'*inv(diag(r))*I_G(:,ind_I{i}(2:end));
    X_inv=I_G(:,ind_I{i}(2:end))'*inv(diag(x))*I_G(:,ind_I{i}(2:end));
    N{i}=length(ind_I{i});
    R{i}=zeros(N{i},N{i});
    X{i}=zeros(N{i},N{i});
    R{i}(2:end,2:end)=inv(R_inv);
    X{i}(2:end,2:end)=inv(X_inv);
    R{i}(1,1)=0.1;
    X{i}(1,1)=0.1;
end

% The on load tap changer resistance and reactance
base_Z=(base_V)^2/base_MVA;
for i=1:subgraph_Number
   RLDC{i}=OLTC(i,4)/base_Z;
   XLDC{i}=OLTC(i,5)/base_Z;
end

% S_subgraph:maps the bus indeces (1,...,N) to the subgraph indeces
% (n1_1,..,n_{N1},n2_1,...,n2_{N2},n3_1,...)
S_subgraph=[];
for i=1:subgraph_Number
    S_subgraph=[S_subgraph;S{i}'];
end

%% Matrices for the MPP code:
% R_tilda=blkdiag(R1,R2,...)
% X_tilda=blkdiag(X1,X2,...)

for i=1:subgraph_Number
   X_tilda=blkdiag(X{:}); 
   R_tilda=blkdiag(R{:}); 
end


save('Feeder_matrices.mat','S','X_tilda','R_tilda','S_subgraph','RLDC','XLDC','L','N','I','subgraph_Number','VR_dest','VR_origin','R','X','S')








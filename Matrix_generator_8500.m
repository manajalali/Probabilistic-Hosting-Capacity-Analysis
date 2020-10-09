%% Reading the 8500-bus feeder:
% This file builds the feeder matrices such as R,X, and incidence matrix
% for the IEEE 8500 subsystems. However, you can modify it to yield all the
% subsystems. 

% This file should be in a folder with the feeder excel files. 
% For running this file, you will need a matlab 2018 or higher. 
% The inputs for this file are:
% 1) The excel files downloaded from the IEEE website. 
% 2) A '.mat' file including the line configurations of the system. 
% 3) In the first 3 lines, you should specify the base voltage and MVA
% ratings. 

% loading the lines, and load files . 
% The bus and line IDs are characters. This code, first changes the format
% of the IDs. Then, we change the IDs from characters to numbers. 
% Each ID is repeated for each pahse, for example, bus SX2767340A,
% SX2767340B and SX2767340C will have the same digital ID. 

clear;

base_kVA = 300;
base_MVA = base_kVA/1000;
base_KV = 7.2;


% Loading the excel files of the IEEE 8500-bus system
lines=readtable('Lines2.xls');
line_length=table2array(lines(:,6));
line_length=convertCharsToStrings(line_length);

line_origin=table2array(lines(:,2));
line_dest=table2array(lines(:,4));
line_conf=table2array(lines(:,10));
line_phases=table2array(lines(:,5));

loads=readtable('Loads.xls');
loads_id=table2array(loads(:,3));
loads_power=table2array(loads(:,9));
loads_id(isnan(loads_power))=[];
loads_power(isnan(loads_power))=[];

% Removing the character files. This way only the digits corresponding to
% each ID remains. 

loads_id = regexp(loads_id,'\d*','Match');
for i=1:length(loads_id)
   if (length(loads_id{i})>1)
       [size1,size2]=cellfun(@size,loads_id{i});
       loads_id{i}=loads_id{i}(find(size2==max(size2)));
   end
end
clear size1 size2
loads_id=string(loads_id);


clear loads 
loads.id=loads_id;
clear loads_id

% Finding multi-phase IDs:
[uniques,idy,idx]=unique(loads.id,'stable');
count=hist(idx,unique(idx));

ind.phase1=find(count==1);
ind.phase2=find(count==2);
ind.phase3=find(count==3);

% Taking the average of powers:
power=zeros(length(idy),1);
for i=1:length(ind.phase1)
   power(ind.phase1(i)) = mean(loads_power(idx==ind.phase1(i)));
end

for i=1:length(ind.phase2)
   power(ind.phase2(i)) = mean(loads_power(idx==ind.phase2(i)));
end

for i=1:length(ind.phase3)
   power(ind.phase3(i)) = mean(loads_power(idx==ind.phase3(i)));
end

unique_loads.power=power;
unique_loads.id=uniques;
unique_loads.phases=count;
clear uniques count idx loads loads_power idy

% Removing the open switches from the lines:
line_dest(2519:2523)=[];
line_origin(2519:2523)=[];
line_conf(2519:2523)=[];
line_phases(2519:2523)=[];

% Removing the substation:
line_dest(1)=[];
line_origin(1)=[];
line_conf(1)=[];
line_phases(1)=[];

line_dest=string(line_dest);
line_origin=string(line_origin);
line_conf=string(line_conf);
%% Merging the connectors:
% In this part we merge both ends of 1-phase and 3-phase connectors and
% change the naming accordingly.

ind_connector=[];
for l=1:length(line_conf)
    if (strcmp(line_conf(l,:),'1PH-Connector') || strcmp(line_conf(l,:),'3PH-Connector'))
        ind_connector=[ind_connector l];
    end
end
line_dest_connector=line_dest(ind_connector);
line_origin_connector=line_origin(ind_connector);
line_origin(ind_connector)=line_dest(ind_connector);

line_dest(ind_connector)=[];
line_origin(ind_connector)=[];
line_conf(ind_connector)=[];


for i=1:length(ind_connector)
    ind_merge_d=find(line_dest==line_origin_connector(i));
    ind_merge_o=find(line_origin==line_dest_connector(i));
    
    if (~isempty(ind_merge_d))
    line_dest(ind_merge_d)=line_dest_connector(i);
    end
   
end
clearvars ind_merge_d ind_merge_o ind_connector line_dest_connector line_origin_connector


% Changing the ID of lines from characters to numbers:
reg_origin = regexp(line_origin,'(regxfmr)','Match');
reg_origin_id=find(~cellfun(@isempty,reg_origin));
line_origin = regexp(line_origin,'\d*','Match');
for i=1:length(line_origin)
   if (length(line_origin{i})>1)
       [size1,size2]=cellfun(@size,line_origin{i});
       line_origin{i}=line_origin{i}(find(size2==max(size2)));
   end
end
clear size1 size2
line_origin=string(line_origin);

reg_dest = regexp(line_dest,'(regxfmr)','Match');
reg_dest_id=find(~cellfun(@isempty,reg_dest));
line_dest = regexp(line_dest,'\d*','Match');
for i=1:length(line_dest)
   if (length(line_dest{i})>1)
       [size1,size2]=cellfun(@size,line_dest{i});
       line_dest{i}=line_dest{i}(find(size2==max(size2)));
   end
end

clear size1 size2
line_dest=string(line_dest);
%%
%  Changing the ID of the buses from characters to numbers:
buses=unique([line_origin;line_dest]);

K=[1:1:length(buses)];

bus_map=containers.Map(buses,K);
k=1;
L = length(line_origin);
origin_ind=zeros(length(line_origin),1);
dest_ind=zeros(length(line_origin),1);


unique_loads.newid=zeros(length(unique_loads.id),1);
for i=1:L
    
    origin_ind(i) = bus_map(line_origin(i));
    dest_ind(i) = bus_map(line_dest(i));
    
    for j=1:length(unique_loads.id)
       if(strcmp(line_dest(i),unique_loads.id(j,:)))
           unique_loads.newid(j)=dest_ind(i);
       end
    end
    
end

line.lines=[origin_ind dest_ind];
line.length=line_length;
line.phase=line_phases;
line.conf=line_conf;
clear line_length line_phases line_conf origin_ind dest_ind lines bus_map

%% making the complete loads:
p_base=zeros(length(buses),1);
p_base(unique_loads.newid)=unique_loads.power;

%% Making the graph of the complete feeder

G=graph(line.lines(:,1),line.lines(:,2));

G=rmedge(G,line.lines(reg_dest_id,1),line.lines(reg_dest_id,2));

[beans,beansizes]=conncomp(G);
%% dividing into subgraphs:
% Here we have chosen the second subgraph. You can change beans==n for
% n=1,2,3,4. 

ind_I=find(K.*(beans==2));
G_sub=subgraph(G,ind_I);

K=[1:1:length(buses)];

% Finding the lines corresponding to the subgraph
subline_ind2=[];
subline_ind=[];
for j=1:length(ind_I)
    subline_ind2=[subline_ind2;find(line.lines(:,1)==K(ind_I(j)))];
    subline_ind=[subline_ind;find(line.lines(:,2)==K(ind_I(j)))];
end
subline=intersect(subline_ind,subline_ind2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The struct "sub" includes all the bus and line information for the chosen
% subgraph. This information includes power, new and old ID of the buses
% and lines, line configuration and lengths.

sub.lines=line.lines(subline,:);
sub.conf=line.conf(subline);
sub.length=line.length(subline);

sub.buses=unique([sub.lines(:,1);sub.lines(:,2)]);

sub.new_id=[1:length(sub.buses)];

sub.p=p_base(ind_I);

new_bus_map=containers.Map(string(sub.buses),sub.new_id);

L_new=length(sub.lines);

for i=1:L_new
    
    origin_ind_new(i) = new_bus_map(string(sub.lines(i,1)));
    dest_ind_new(i) = new_bus_map(string(sub.lines(i,2)));
    
    for j=1:length(ind_I)
        if(sub.lines(i,2)==ind_I(j))
            sub.newbus(j)=dest_ind_new(i);
            sub.newp(j)=p_base(ind_I(j));
        end
    end
end
sub.new_lines=[origin_ind_new' dest_ind_new'];
clear origin_ind_new dest_ind_new

%% Building the R and X matrix for the subgraph
load('8500configs.mat')
type=zeros(L,1);
ind_connector=[];
for l=1:L_new

    R_3p{l}=str2double(sub.length(l))*Rmatrix{str2double(sub.conf(l))};
    R_3p{l}=Rmatrix{str2double(sub.conf(l))};
    R_1p{l}=sum(diag(R_3p{l}))/(nnz(diag(R_3p{l})));
    
    X_3p{l}=str2double(sub.length(l))*Xmatrix{str2double(sub.conf(l))};
    X_3p{l}=Xmatrix{str2double(sub.conf(l))};
    X_1p{l}=sum(diag(X_3p{l}))/(nnz(diag(X_3p{l})));
 end


r=cell2table(R_1p);
r=table2array(r);
x=cell2table(X_1p);
x=x{:,:};

%% Building the incidence matrix for the chosen subgraph
A=zeros(L_new,max(sub.newbus));

for i=1:L_new
   
    A(i,sub.new_lines(i,1))=1;
    A(i,sub.new_lines(i,2))=-1;
    
end

%%% For a reduced incidence matrix, the substation bus is removed. Usually
%%% the substation is bus 1. However, here the substation is the voltage
%%% regulator's bus. To build a reduced incidence matrix, we first find the
%%% ID of the voltage regulators. 
reg_ind=[];
for j=1:length(reg_dest_id)
   reg_name=line.lines(reg_dest_id(j),2);
   reg_ind=[reg_ind find(sub.lines(:,1)==reg_name)];
end
reg_ind=new_bus_map(string(sub.lines(reg_ind,1)))

Ar=A;
Ar(:,reg_ind)=[];

R=inv(Ar)*diag(r)*inv(Ar');
X=inv(Ar)*diag(x)*inv(Ar');

N=max(sub.new_id);
L=L_new;

p_base=sub.p;


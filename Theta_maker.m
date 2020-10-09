function [Theta] = Theta_maker(pg,pc,qc,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This current body of code generates Theta for all injection scalings and
% for 10 solar peneteration values (10-100) and with and without inverter
% oversizing. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta=[];

for i = 1:scale
    
    pg_max_scaled=max(pg');
    pg_max_scaled=repmat(pg_max_scaled',1,size(pg,2));
    sg_max=1.1*pg_max_scaled; % 10% inverter oversizing
    qg_max=(realsqrt((sg_max.^2)-((pg).^2)));
     
    for m_a = 1:10
        m_alpha=m_a/10;
        Theta=[Theta scale.*[m_alpha*pg-pc;m_alpha*qg_max-qc;-m_alpha*qg_max-qc]];
    end
    
    pg_max_scaled=max(pg');
    pg_max_scaled=repmat(pg_max_scaled',1,size(pg,2));
    sg_max=1*pg_max_scaled; % No inverter oversizing
    qg_max=(realsqrt((sg_max.^2)-((pg).^2)));
    
    for m_a = 1:10
        m_alpha=m_a/10;
        Theta=[Theta scale.*[m_alpha*pg-pc;m_alpha*qg_max-qc;-m_alpha*qg_max-qc]];
    end
end
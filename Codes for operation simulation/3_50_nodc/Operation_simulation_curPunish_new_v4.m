load load_data 
load solar_data_by_bus
load wind_data_by_bus
% load('OSdata.mat')
%VX=ones(44,1)';
load HRP-38-50
simplifyModel
%% NOMENCLATURE
% Number of each element 
isExtremeOn=0;
bmva=mpc.baseMVA;
flag=0;
nlc=size(mpc.branchcAC,1);nle=size(mpc.branche,1);
nw=size(mpc.RE,1);
index_wind=find(mpc.RE(:,12)==4);nwind=size(wind_data,3);
index_solar=find(mpc.RE(:,12)==3);nsolar=size(solar_data,3);
ni=size(mpc.gen,1);
ns=1;
nt=24;
nb=size(mpc.bus,1);
% nj=365;
nop=24;% amount of operation condition
dj=365/(nj)*ones(1,nj); 
%--------------Hydro-------------
idx_ScenMonth=spotmonth(1:365);
%--------------Hydro-------------
if nj~=365
   kmeanCluster 
%--------------Hydro-------------
    idx_ScenMonth=spotmonth(idx_c);
%--------------Hydro-------------   
end
%T=mpc.branchc(:,6)*1.5;T(15:16)=[1500,1500];
T=400000*ones(nlc,1);
makeM
c_lshed=3;shedL=1;
% c_reCur=0.00;
Gen_rate=1;
ny=1;% number of planning years
nv=15;
% Sets
% map
mib=[(1:ni)',mpc.gen(:,1)];mwb=[(1:nw)',mpc.RE(:,1)];
mlesb=[(1:nle)',mpc.branche(:,1)];
mlerb=[(1:nle)',mpc.branche(:,2)];
mlcsb=[(1:nlc)',mpc.branchcAC(:,1)];
mlcrb=[(1:nlc)',mpc.branchcAC(:,2)];
bl=intersect(setdiff(1:nb,[mpc.gen(:,1);mpc.RE(:,1)]),find(mpc.bus(:,3)~=0))'; % node connect with load but without generators
bg=intersect(find(mpc.bus(:,3)==0),mpc.gen(:,1))'; % node connect with generators but without load
bnl=find(mpc.bus(:,3)==0);
% Constants
discount_rate=0.08;recovery_year=25;
RY=(((1+discount_rate)^recovery_year*discount_rate)/((1+discount_rate)^recovery_year-1));
c_line_lc=repmat(mpc.branchcAC(:,14),1,ny)*RY;
c_gen_i=repmat(mpc.gencost(:,8),1,ny);
tau=1;
alpha=1;
lg=1;
 
xle=mpc.branche(:,4);xlc=mpc.branchcAC(:,4);
p_gmax_i=mpc.gen(:,9)*Gen_rate;p_gmin_i=mpc.gen(:,10)*0;
r_g=p_gmin_i./p_gmax_i;
r_u_i=mpc.gen(:,20);r_d_i=mpc.gen(:,20);
f_max_le=mpc.branche(:,6);
f_max_lc=mpc.branchcAC(:,6);
p_f_w=zeros(nw,ns,nj,nt,ny);
% w_cap=[20000 50000 9800 20000 70000];

rand('seed',24);

idx_op=1:nj;
for iw = 1 : nw
    for iy = 1:ny
        for is =1:ns
            for ij= 1:nj
                if ismember(iw,index_wind)
                    index_area=mpc.RE(iw,11);
                    p_f_w(iw,is,ij,:,iy)=wind_data( : ,idx_op(ij) ,iw-nsolar)*mpc.RE(iw,9);  
                elseif ismember(iw,index_solar)
                    index_area=mpc.RE(iw,11); 
                    p_f_w(iw,is,ij,:,iy)=solar_data( : ,idx_op(ij) ,iw)*mpc.RE(iw,9);  
                end  
%                    p_f_w(iw,is,ij,:,iy)=wind_data( : ,idx_op(ij) ,mod(mpc.RE(iw,11)-1,(size(wind_data,3)))+1 )*mpc.RE(iw,9);  
            end
        end
    end
end
l_b=zeros(nb,ns,nj,nt,ny);
gr=0.01;
for ib = 1 : nb
    if mpc.bus(ib,3)~=-1
        for is = 1 : ns
            for ij= 1:nj
                for iy= 1 : ny
                    index_area=mpc.bus(ib,7);% **
                    l_b(ib,is,ij,:,iy)=( load_data(:, idx_op(ij) ,index_area)*lg(is)*mpc.bus(ib,3))*(1+gr)^(iy-1);
%                     l_b(ib,is,ij,:,iy)=( load_data(:, idx_op(ij) ,1)*lg(is)*mpc.bus(ib,3))*(1+gr)^(iy-1);
                end
            end
        end
    end
end

% Variables
z=sdpvar(1); % relaxed master objective variable;
Int=binvar(2*nlc*ny,1); % decision of building lines 
%% BUILD MATRIX FOR SUBPROBLEM
Xv=ones(nlc,ny); % initial value
Iv=[ones(nlc,1),zeros(nlc,ny-1)]; % initial value

r_t=(nb+nle+nle*2+nlc*2+nlc*2+nw*2+2*nb+2*ni+1);c_t=(ni+nw+nle+nlc+nb+nb);
r_calu=r_t*nt;c_calu=c_t*nt;
r_scen=r_calu*nj;c_scen=c_calu*nj;
r_year=r_scen*ns;c_year=c_scen*ns;
Row=r_year*ny;Col=c_year*ny;

A_calu=sparse(r_calu,c_calu);B_calu=zeros(r_calu,1);C_calu=zeros(1,c_calu);sig_calu=zeros(r_calu,1);
A_t=sparse(r_t,c_t);B_t=zeros(r_t,1);E_t=sparse(r_t,nv);C_t=zeros(1,c_t);sig_t=zeros(r_t,1);
A_ramp=sparse(2*ni,2*c_t);B_ramp=sparse(ones(2*ni,1));sig_ramp=sparse(ones(2*ni,1));
B_t_lc=sparse(4*nlc,1);

%A_all=sparse(r_scen,c_scen);
Ap_t=sparse(r_t,nlc+nlc);

i_pb=1 : nb;
i_pfle=nb+(1 : nle);
i_cle_1=nb+nle+(1:nle);i_cle_2=nb+nle+nle+(1:nle);
i_pflc_1=nb+nle+2*nle+(1:nlc);i_pflc_2=nb+nle+2*nle+nlc+(1:nlc);
i_clc_1=nb+nle+2*nle+2*nlc+(1:nlc);i_clc_2=nb+nle+2*nle+2*nlc+nlc+(1:nlc);
i_wf_1=nb+nle+2*nle+2*nlc+2*nlc+(1:nw);i_wf_2=nb+nle+2*nle+2*nlc+2*nlc+nw+(1:nw);
i_shed_1=nb+nle+2*nle+2*nlc+2*nlc+2*nw+(1:nb);i_shed_2=nb+nle+2*nle+2*nlc+2*nlc+2*nw+nb+(1:nb);
i_tu_1=nb+nle+2*nle+2*nlc+2*nlc+2*nw+2*nb+(1:ni);i_tu_2=nb+nle+2*nle+2*nlc+2*nlc+2*nw+2*nb+ni+(1:ni);
i_ref=nb+nle+2*nle+2*nlc+2*nlc+2*nw+2*nb+ni*2+1;

x_x=1:nlc;
x_i=nlc+(1:nlc);

i_i=1:ni;
%--------------Hydro-------------
i_hi=mpc.hydro.index;n_hydro=length(i_hi);
%--------------Hydro-------------
i_w=ni+(1:nw);
i_le=ni+nw+(1:nle);
i_lc=ni+nw+nle+(1:nlc);
i_vpa=ni+nw+nle+nlc+(1:nb);
i_lsh=ni+nw+nle+nlc+nb+(1:nb);

ncalu=ny*ns*nj;
Bu=zeros(r_t,ncalu*nt);
Cu=zeros(c_t,ncalu*nt);
Sigu=zeros(r_t,ncalu*nt);
LWS=zeros(nv,ncalu*nt);

% powerbalance
A_t(i_pb,i_i)=sparse(mib(:,2),mib(:,1),ones(ni,1),nb,ni);
A_t(i_pb,i_w)=sparse(mwb(:,2),mwb(:,1),ones(nw,1),nb,nw);
A_t(i_pb,i_le)=sparse([mlerb(:,2);mlesb(:,2)],...
                [mlerb(:,1);mlesb(:,1)],...
                [ones(nle,1);ones(nle,1)*-1],...
                nb,nle);
A_t(i_pb,i_lc)=sparse([mlcrb(:,2);mlcsb(:,2)],...
                [mlcrb(:,1);mlcsb(:,1)],...
                [ones(nlc,1);ones(nlc,1)*-1],...
                nb,nlc);
A_t(i_pb,i_lsh)=sparse(i_pb,i_pb,ones(nb,1),nb,nb);

% dcpf of built lines
A_t(i_pfle,i_le)=sparse(1:nle,1:nle,ones(nle,1),nle,nle);
A_t(i_pfle,i_vpa)=sparse([1:nle,1:nle]',...
                [mpc.branche(:,1);mpc.branche(:,2)],...
                [-1./mpc.branche(:,4)*bmva;1./mpc.branche(:,4)*bmva],...
                nle,nb);
            
% capacity bounds of built lines
A_t(i_cle_1,i_le)=sparse(1:nle,1:nle,ones(nle,1),nle,nle);
A_t(i_cle_2,i_le)=sparse(1:nle,1:nle,-ones(nle,1),nle,nle);

% dcpf of candidate lines
% lower
A_t(i_pflc_1,i_lc)=sparse(1:nlc,1:nlc,ones(nlc,1),nlc,nlc);
A_t(i_pflc_1,i_vpa)=sparse([1:nlc,1:nlc]',...
                [mpc.branchcAC(:,1);mpc.branchcAC(:,2)],...
                [-1./mpc.branchcAC(:,4)*bmva;1./mpc.branchcAC(:,4)*bmva],...
                nlc,nb);
% upper
A_t(i_pflc_2,i_lc)=sparse(1:nlc,1:nlc,-ones(nlc,1),nlc,nlc);
A_t(i_pflc_2,i_vpa)=sparse([1:nlc,1:nlc]',...
                [mpc.branchcAC(:,1);mpc.branchcAC(:,2)],...
                [1./mpc.branchcAC(:,4)*bmva;-1./mpc.branchcAC(:,4)*bmva],...
                nlc,nb);
            
Ap_t(i_pflc_1,x_x)=sparse(1:nlc,1:nlc,-T,nlc,nlc);
Ap_t(i_pflc_2,x_x)=sparse(1:nlc,1:nlc,-T,nlc,nlc);            
            
% capacity bounds of candidate lines
A_t(i_clc_1,i_lc)=sparse(1:nlc,1:nlc,ones(nlc,1),nlc,nlc);
A_t(i_clc_2,i_lc)=sparse(1:nlc,1:nlc,-ones(nlc,1),nlc,nlc);

Ap_t(i_clc_1,x_x)=sparse(1:nlc,1:nlc,f_max_lc,nlc,nlc);
Ap_t(i_clc_2,x_x)=sparse(1:nlc,1:nlc,f_max_lc,nlc,nlc);

% bounds of wind turbines
% lower
A_t(i_wf_1,i_w)=sparse(1:nw,1:nw,ones(nw,1),nw,nw);
% upper
A_t(i_wf_2,i_w)=sparse(1:nw,1:nw,-ones(nw,1),nw,nw);

% bounds of load shedding
A_t(i_shed_1,i_lsh)=sparse(1:nb,1:nb,ones(nb,1),nb,nb);
A_t(i_shed_2,i_lsh)=sparse(1:nb,1:nb,-ones(nb,1),nb,nb);

% bounds of thermal units
A_t(i_tu_1,i_i)=sparse(1:ni,1:ni,ones(ni,1),ni,ni);
A_t(i_tu_2,i_i)=sparse(1:ni,1:ni,-ones(ni,1),ni,ni);

% reference node vpa
A_t(i_ref,i_vpa(1)-1+13)=1;


for iy = 1 : ny
    for is = 1 : ns
        for ij = 1 : nj
            for it = 1 : nt 
                index_calu=(iy-1)*ns*nj*nt+(is-1)*nj*nt+(ij-1)*nt+it;
                LWS(:,index_calu)=[l_b(mpc.region_index,is,ij,it,iy)./mpc.bus(mpc.region_index,3);...
                    p_f_w(mpc.theta_index,is,ij,it,iy)./mpc.RE(mpc.theta_index,9)];
                %************************CalUnit**************************************************
                %***************************** B_t & sig_t *********************************************
                % powerbalance
                B_t(1:nb,1)=l_b(:,is,ij,it,iy);
                sig_t(1:nb,1)=ones(nb,1)*2;
                
                % dcpf of built lines
                sig_t(nb+(1:nle),1)=ones(nle,1)*2;
                
                % capacity bounds of built lines
                % lower bound               
                B_t(i_cle_1,1)=-f_max_le;
                sig_t(i_cle_1,1)=ones(nle,1)*1;
                % upper bound             
                B_t(i_cle_2,1)=-f_max_le;
                sig_t(i_cle_2,1)=ones(nle,1)*1;

                % dcpf of candidate lines
                % lower bound
                B_t(i_pflc_1,1)=-T*(1);
                sig_t(i_pflc_1,1)=ones(nlc,1)*4;
                % upper bound 
                B_t(i_pflc_2,1)=-T*(1);
                sig_t(i_pflc_2,1)=ones(nlc,1)*4;
                
                % capacity bounds of candidate lines
                % lower bound             
                B_t(i_clc_1,1)=-f_max_lc*0;
                sig_t(i_clc_1,1)=ones(nlc,1)*5;
                % upper bound               
                B_t(i_clc_2,1)=-f_max_lc*0;
                sig_t(i_clc_2,1)=ones(nlc,1)*5;
                
                % bounds of wind turbines
                % lower bound  
                B_t(i_wf_1,1)=p_f_w(:,is,ij,it,iy)*0;
                sig_t(i_wf_1,1)=ones(nw,1)*1;  
                % upper bound             
                B_t(i_wf_2,1)=-p_f_w(:,is,ij,it,iy);
                sig_t(i_wf_2,1)=ones(nw,1)*1; 
                
                % bounds of load shedding
                % lower bound              
                sig_t(i_shed_1,1)=ones(nb,1)*1;  
                % upper bound             
                B_t(i_shed_2,1)=-l_b(:,is,ij,it,iy)*shedL;%***********************************************
                sig_t(i_shed_2,1)=ones(nb,1)*1;     
                
                % bounds of thermal units
                % lower bound      
                B_t(i_tu_1,1)=p_gmin_i(:);
                sig_t(i_tu_1,1)=ones(ni,1)*1;  
                % upper bound             
                B_t(i_tu_2,1)=-p_gmax_i(:);
                sig_t(i_tu_2,1)=ones(ni,1)*1;    
                
                % reference node vpa
                sig_t(i_ref,1)=2;  
                
                % objective function of calculation unit
                C_t(1,i_i)=c_gen_i(i_i);
                C_t(1,i_lsh)=ones(1,nb)*c_lshed; %!!!!!!!!!!!!!!!!!!!!!!!!!!
                C_t(1,i_w)=ones(1,nw)*-c_reCur;                
                C_t=C_t*alpha(is)*dj(ij);%***************************************365********************************************
                %*****************************A_t*********************************************
                B_calu((it-1)*r_t+1:it*r_t,1)=B_t;
                sig_calu((it-1)*r_t+1:it*r_t,1)=sig_t;
                C_calu(1,(it-1)*c_t+1:it*c_t)=C_t;
                
            Bu(:,index_calu)=B_t'; 
            Cu(:,index_calu)=C_t';
            Sigu(:,index_calu)=sig_t;
            %************************CalUnit**************************************************   
            end
        end
    end
end


B0=Bu(1,:);B0([1:nb,i_shed_2,i_wf_2])=0;

index1=find(sig_t~=2); % index of inequality 
index2=find(sig_t==2); % index of equality

index4=find(sig_t==4); % index of constraints of candidate lines
index5=find(sig_t==5); % index of constraints of candidate lines

index_const=find(or(sig_t==1,sig_t==2));
index_coeff=union(index4,index5);

coeff=repmat([T;T;repmat(-mpc.branchcAC(:,6),2,1)]',1,1);


for iy = 1 : ny
    coeffX_mtx([1:ns*nj*nt]+(iy-1)*ns*nj*nt,1:nlc*2*2)=repmat(coeff,ncalu*nt/ny,1);
end

Atv=A_t';
Atv_ineq=sparse(1:size(index1,1),index1,-1*ones(size(index1,1),1),size(index1,1),r_t);

%% BUILD MATRIX FOR MASTER PROBLEM
% Relatoionship of I & X 
mtx_cell_Ameq=cell(ny,1);
for iy= 1 : ny
mtx_cell_Ameq(iy)={sparse([(1:nlc)';(1:nlc)'],[(1:nlc)';nlc+(1:nlc)'],[ones(nlc,1);-ones(nlc,1)],ny*nlc,ny*(nlc*2))};
    if iy >1 
        mtx_cell_Ameq(iy)={sparse([(1:nlc)'+(iy-1)*nlc; (1:nlc)'+(iy-1)*nlc; (1:nlc)'+(iy-1)*nlc],...
            [(1:nlc)'+(iy-1)*nlc*2; nlc+(1:nlc)'+(iy-1)*nlc*2; (1:nlc)'+(iy-2)*nlc*2],...
            [ones(nlc,1); -ones(nlc,1); -ones(nlc,1)],...
            ny*nlc,ny*(nlc*2))};
    end    
end
Am_eq=sparse(ny*nlc,ny*(nlc*2));
for iy = 1 : ny 
    Am_eq=Am_eq+mtx_cell_Ameq{iy};
end
Bm_eq=sparse(ny*nlc,1);

% Lines will be built once at most(BO)
% cons=[cons,(sum(I,2)<=1):'I'];
Am_ineq_BO=[repmat(sparse((1:nlc)',nlc+(1:nlc)',ones(nlc,1),nlc,nlc*2),1,ny)];
Bm_ineq_BO=ones(nlc,1);

% % the minimum demand of line capacity at each load node (BL)
nbl=size(bl,2); % number of pure load nodes 
nbg=size(bl,2); % number of pure generator nodes 
Am_ineq_BL=[];Bm_ineq_BL=[];
for iy =1 : ny
    for  ibl= bl
        lbe=find(or(mpc.branche(:,1)==ibl,mpc.branche(:,2)==ibl));
        lbc=find(or(mpc.branchcAC(:,1)==ibl,mpc.branchcAC(:,2)==ibl));
        nlbc=size(lbc,1);
        if isempty(mpc.branchcAC(lbc,6))==0
            Am_ineq_BL=[Am_ineq_BL;
                sparse(ones(1,nlbc),lbc+(iy-1)*nlc*2,-mpc.branchcAC(lbc,6),1,ny*(nlc*2))];
            Bm_ineq_BL=[Bm_ineq_BL;sum(mpc.branche(lbe,6))-mpc.bus(ibl,3)];
        end
    end
end

Am_ineq_BG=[];Bm_ineq_BG=[];
for iy =1 : ny
    for  ibg= bg
        lgbe=find(or(mpc.branche(:,1)==ibg,mpc.branche(:,2)==ibg));
        lgbc=find(or(mpc.branchcAC(:,1)==ibg,mpc.branchcAC(:,2)==ibg));
        nlgbc=size(lgbc,1);
        if isempty(mpc.branchcAC(lgbc,6))==0         
            Am_ineq_BG=[Am_ineq_BG;
                sparse(ones(1,nlgbc),lgbc+(iy-1)*nlc*2,-mpc.branchcAC(lgbc,6)',1,ny*(nlc*2))];
            Bm_ineq_BG=[Bm_ineq_BG;sum(mpc.branche(lgbe,6))-sum(mpc.gen(find(mpc.gen(:,1)==ibg),10))];        
        end
    end
end

Am_ineq=[Am_ineq_BO;Am_ineq_BL;Am_ineq_BG];
Bm_ineq=[Bm_ineq_BO;Bm_ineq_BL;Bm_ineq_BG];

Cm=[c_line_lc;zeros(nlc,1)];
%% Gurobi III
%--------------Hydro-------------
list_month=unique(idx_ScenMonth);
n_month=length(list_month);
result=[];

for imonth= 1 : n_month
    
    idx_mon=find(idx_ScenMonth==list_month(imonth));
    nj_mon=length(idx_mon);      

    location_in_all_r=(1+(idx_mon(1)-1)*r_calu):(idx_mon(end))*r_calu;
    location_in_all_c=(1+(idx_mon(1)-1)*c_calu):(idx_mon(end))*c_calu;
    
    [A_t_r,A_t_c,A_t_v]=find(A_t);
    temp=repmat(0:r_t:r_t*(nt*nj_mon-1),size(A_t_r,1),1);
    A_calu_r=repmat(A_t_r,nt*nj_mon,1)+temp(:);
    temp=repmat(0:c_t:c_t*(nt*nj_mon-1),size(A_t_c,1),1);
    A_calu_c=repmat(A_t_c,nt*nj_mon,1)+temp(:);
    A_calu_v=repmat(A_t_v,nt*nj_mon,1);    


    A_hydro=sparse(n_hydro*1,c_t*nt*nj_mon); 
    B_hydro=sparse(n_hydro*1,1);

  
    
    for ih= 1 : n_hydro
        n_row_h=ih;
        month_h=list_month(imonth);
        idx_j_h=find(idx_ScenMonth==month_h);
        temp1=((1:length(idx_j_h))-1)*c_calu+i_hi(ih);
        temp2=((1:24)'-1)*c_t;
        temp3=repmat(temp1,24,1)+repmat(temp2,1,length(idx_j_h));
        temp4=repmat(dj(idx_j_h),24,1);
        A_hydro(n_row_h,:)=sparse(1,temp3(:),temp4(:),1,c_t*nt*nj_mon);
        B_hydro(n_row_h,1)=mpc.hydro.dry(ih,imonth)*24*sum(dj(idx_j_h));
    end
    
clear('A_all')
A_all(1:nt*nj_mon*r_t,1:nt*nj_mon*c_t)=sparse(A_calu_r,A_calu_c,A_calu_v);

Ap_all=repmat(Ap_t,nt*nj_mon,1);

A=sparse(size(Am_eq,1)+size(Am_ineq,1)+r_t*nt*nj_mon,2*nlc+c_t*nt*nj_mon);
A(1:size(Am_eq,1)+size(Am_ineq,1),1:2*nlc)=[Am_eq;Am_ineq];
A(size(Am_eq,1)+size(Am_ineq,1)+(1:r_t*nt*nj_mon),1:2*nlc)=Ap_all;
A(size(Am_eq,1)+size(Am_ineq,1)+(1:r_t*nt*nj_mon),2*nlc+(1:c_t*nt*nj_mon))=A_all;

A_ineq1=A(size(Am_eq,1)+(1:size(Am_ineq,1)),:);
A_eq1=A(1:size(Am_eq,1),:);

index_ineq=find(Sigu(location_in_all_r)~=2);
index_eq=find(Sigu(location_in_all_r)==2);

B2=Bu(location_in_all_r);
B_ineq2=B2(index_ineq);
B_eq2=B2(index_eq)';

A2=A(size(Am_eq,1)+size(Am_ineq,1)+(1:r_t*nt*nj_mon),:);
A_ineq2=A2(index_ineq,:);
A_eq2=A2(index_eq,:);

X=binvar(nlc,ny); % State of lines 
I=binvar(nlc,ny); % decision of building lines 

C=sdpvar(c_t*nt*nj_mon,ny);


Xv=VX';
Iv=Xv;


temp=repmat(dj,24,1);
reUpper=diag(temp(:))*-Bu(i_wf_2,:)';

Cu=Cu(:);

%%%%
idx_bin=1:2*nlc;
idx_c=2*nlc+(1:size(C,1));

if size(B_ineq2,1)==1
    B_ineq2=B_ineq2';
end
if size(B_eq2,1)==1
    B_eq2=B_eq2';
end

A_ineq=[A_ineq1(:,idx_c);-A_ineq2(:,idx_c);A_hydro];
B_ineq=[Bm_ineq-A_ineq1(:,idx_bin)*([Xv;Xv]);...
        -B_ineq2+A_ineq2(:,idx_bin)*([Xv;Xv]);...
        B_hydro];
    
A_eq=[A_eq1(:,idx_c);A_eq2(:,idx_c)];
B_eq=[Bm_eq-A_eq1(:,idx_bin)*([Xv;Xv]);...
      B_eq2-A_eq2(:,idx_bin)*([Xv;Xv]);...
        ];
%%%%

%-----separate----and---integrate-----

options = cplexoptimset('Display', 'off', 'Algorithm', 'interior-point');
[cplexS_a,cplexS_b,cplexS_c,~,cplexS_d] = ...
          cplexlp(Cu(location_in_all_c),A_ineq,B_ineq,A_eq,B_eq,[],[],[],options);
fprintf(' Month= %1.0f\n',list_month(imonth));
%-----separate----and---integrate-----
result=[result,reshape(cplexS_a,c_t,nt*nj_mon)];
obj_mon(imonth)=cplexS_b;
end
%--------------Hydro-------------
    
    objsum=sum(obj_mon)+c_line_lc'*Iv+sum(sum(reUpper))*c_reCur;
    
    
    actualOutput=diag(temp(:))*result(i_w,:)';
    reCurtail=reUpper-actualOutput;
    reCur_energy=sum(sum(reCurtail));
    rate=reCur_energy/sum(sum(reUpper));
    reCurPunish=reCur_energy*c_reCur;
    
    loadTotal=diag(temp(:))*Bu(i_pb,:)';
    loadTotal_energy=sum(sum(loadTotal));
    
    loadshed=diag(temp(:))*result(i_lsh,:)';
    loadshed_energy=sum(sum(loadshed));   
    
    fprintf(' total cost = %1.4f\n investment cost = %1.4f\n operating cost = %1.4f Unit:10^9 \n', ...
        objsum/1.0e+05,(c_line_lc'*Iv)/1.0e+05,(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05)
    fprintf(' Rate of RE curtailment= %1.4f\n',rate*100);
    fprintf(' Load shedding MWh= %1.4f\n',loadshed_energy);
    fprintf(' RE curtailment 10^9 CNY= %1.4f\n',reCurPunish/1.0e+05);


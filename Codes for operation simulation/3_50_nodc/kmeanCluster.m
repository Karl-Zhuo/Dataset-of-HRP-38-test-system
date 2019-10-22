%% By five region pu
load wind_region
load solar_region

njA=365;nt=24;
nwind=size(wind_data,3);
nsolar=size(solar_data,3);
rand('seed',9);
LOAD=zeros(nt*5,njA);
WIND=zeros(nt*5,njA);
SOLAR=zeros(nt*5,njA);
datal4draw=[];
dataw4draw=[];
datas4draw=[];
for ij = 1 : njA
    temp=load_data(:,ij,:);
    LOAD(:,ij)=temp(:);
    datal4draw=[datal4draw;reshape(LOAD(:,ij),24,5)];
    
    temp=wind_data(:,ij,:);
    for i = 1 : 5
    WIND((i-1)*24+(1:24),ij)=sum(wind_data(:,ij,wind_region(:,2)==i),3);
    end
    dataw4draw=[dataw4draw;reshape(WIND(:,ij),24,5)];
    
    temp=solar_data(:,ij,:);
    for i = 1 : 5
    SOLAR((i-1)*24+(1:24),ij)=sum(solar_data(:,ij,solar_region(:,2)==i),3);
    end
    datas4draw=[datas4draw;reshape(SOLAR(:,ij),24,5)];   
end
idx_op=1:floor(365/(nj)):365;
temp=[LOAD;WIND;SOLAR]';
% [idx_d,C] = kmeans([LOAD;WIND;SOLAR]',nj,'Start',temp(idx_op(1:nj),:));
% [idx_d,C] = kmedoids([LOAD;WIND;SOLAR]',nj,'Start',temp(idx_op(1:nj),:));
[idx_d,C] = kmedoids([LOAD;WIND;SOLAR]',nj,'Start','plus','Options',statset('MaxIter',10000));


[~,idx_c]=ismember(C,temp,'rows');
[idx_c,idx_cdj]=sort(idx_c);
load_data_new=zeros(nt,nj,5);
wind_data_new=zeros(nt,nj,nwind);
solar_data_new=zeros(nt,nj,nsolar);

    load_data_new=load_data(:,idx_c,:);
    wind_data_new=wind_data(:,idx_c,:);
    solar_data_new=solar_data(:,idx_c,:);


% for i = 1 : nj
% dj(i)=size(find(idx_d==i),1);
% end


for i = 1 : nj
dj(i)=size(find(idx_d==idx_cdj(i)),1);
end

if isExtremeOn==1
    load('index_Extreme.mat')   
    njextreme=size(index_Extreme,1);
    for i = nj+(1 : njextreme)
        load_data_new(:,i,:)=load_data(:,i,:);
        wind_data_new(:,i,:)=wind_data(:,i,:);
        solar_data_new(:,i,:)=solar_data(:,i,:);
    end

    for i = nj + (1 : njextreme)
    dj(i)=1;
    end   
    nj=nj+njextreme;
end


load_data=load_data_new;
wind_data=wind_data_new;
solar_data=solar_data_new;


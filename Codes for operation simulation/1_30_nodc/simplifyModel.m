% merge the candidate lines in the same corridor
% clc,clear
% load NWG 
bmva=mpc.baseMVA;
flag=0;
nlc=size(mpc.branchcAC,1);nle=size(mpc.branche,1);

branche=mpc.branche(:,1:8);
uncheckList=1:nle;
j=0;
for i = 1 : nle
    if ismember(i,uncheckList)
        j=j+1;
        uncheckList=setdiff(uncheckList,i);
        [Lia,Locb]=...
            ismember(branche(uncheckList,[1,2]),branche(i,[1,2]),'rows');
        nSameCor=sum(Lia);
        if nSameCor
            brancheNew(j,[1,2])=branche(i,[1,2]);
            brancheNew(j,[3,4])=ones(1,2)./...
              sum(ones(nSameCor+1,2)./branche([i,uncheckList(Lia)],[3,4]));
            brancheNew(j,[5:8])=...
              sum(branche([i,uncheckList(Lia)],[5:8]));   
            uncheckList=setdiff(uncheckList,uncheckList(Lia));
        else
            brancheNew(j,:)=branche(i,:);
        end
    end
end
mpc.branche=brancheNew;

branchcAC=mpc.branchcAC;
uncheckList=1:nlc;
j=0;
for i = 1 : nlc
    if ismember(i,uncheckList)
        uncheckList=setdiff(uncheckList,i);
        [Lia,Locb]=...
            ismember(branchcAC(uncheckList,[1,2]),branchcAC(i,[1,2]),'rows');
        nSameCor=sum(Lia);
        nLineAdd=ceil((nSameCor+1)/2);
        if mod(nSameCor,2)==1 % even 
            branchcACNew(j+1:j+nLineAdd,[1,2])=repmat(branchcAC(i,[1,2]),nLineAdd,1);
            for il = 1 : nLineAdd
                indexSame=[i,uncheckList(Lia)];
                branchcACNew(j+il,[3,4])=ones(1,2)./...
                  sum(ones(2,2)./branchcAC(indexSame(il*2-1:il*2),[3,4]));
                branchcACNew(j+il,[5:14])=...
                  sum(branchcAC(indexSame(il*2-1:il*2),[5:14]));                
            end
            j=j+nLineAdd;
            uncheckList=setdiff(uncheckList,uncheckList(Lia));
        elseif mod(nSameCor,2)==0 % odd 
            j=j+1;
            branchcACNew(j:j+nLineAdd-1,[1,2])=repmat(branchcAC(i,[1,2]),nLineAdd,1);
            branchcACNew(j,[3,4])=branchcAC(i,[3,4]);
            branchcACNew(j,[5:14])=branchcAC(i,[5:14]);      
            for il = 1 : nLineAdd-1
                indexSame=uncheckList(Lia);
                branchcACNew(j+il,[3,4])=ones(1,2)./...
                  sum(ones(2,2)./branchcAC(indexSame(il*2-1:il*2),[3,4]));
                branchcACNew(j+il,[5:14])=...
                  sum(branchcAC(indexSame(il*2-1:il*2),[5:14]));             
            end
            j=j+nLineAdd-1;
            uncheckList=setdiff(uncheckList,uncheckList(Lia));
        end
    end
end
mpc.branchcAC=branchcACNew;
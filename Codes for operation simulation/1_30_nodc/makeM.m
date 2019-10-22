% make M 
[Lia,Locb]=ismember(mpc.branchcAC(:,[1,2]),...
                    mpc.branche(:,[1,2]),'rows');
for i = 1 : nlc
    if Lia(i)
        index_exist=Locb(i);
        T(i,1)= mpc.branche(index_exist,4)*mpc.branche(index_exist,6)/...
            mpc.branchcAC(i,4)*1.1 ;       
    end
end
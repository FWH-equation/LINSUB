function Cmprd()
global UR UI NR XR XI CR CI NP;
CR=ones(5,NR);
CI=ones(5,NR);
for i=1:NR
    for j=1:NR
        CR(i,j)=0;
        CI(i,j)=0;
        for k=1:NP
            CR(i,j)=CR(i,j)+XR(i,k)*UR(k,j)-XI(i,k)+UI(k,j);
            CI(i,j)=CI(i,j)+XR(i,k)*UI(k,j)+XI(i,k)*UR(k,j);
        end
    end
end
return
end


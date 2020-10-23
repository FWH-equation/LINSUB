function Cmdiv()
% Cmdiv( KR, KI, UR, UI, NP, NR, 15, 15 )
% CMDIV( AR, AI, BR, BI, IRA, ICB, IRAD, IRBD )
global KR KI NP UR UI NR;
%Loop for pivot row
for IP=1:NP
    %Choose largest pivot for pivot column
    Z=0;
    for i=IP:NP
        X=KR(i,IP)^2+KI(i,IP)^2;
        if (X<Z)
            continue;
        end
        Z=X;
        IE=i;
    end
    if (Z==0)
        %Error message
        error('Zero pivot at line %d',IP); 
    end
    %Exchange rows,scaling new pivot row
    ZR=KR(IE,IP)/Z;
    ZI=-KI(IE,IP)/Z;
    if (IP~=NP)
        KR(IE,IP)=KR(IP,IP);
        KI(IE,IP)=KI(IP,IP);
        K=IP+1;
        for j=K:NP
            XR=KR(IE,j)*ZR-KI(IE,j)*ZI;
            XI=KR(IE,j)*ZI+KI(IE,j)*ZR;
            KR(IE,j)=KR(IP,j);
            KI(IE,j)=KI(IP,j);
            KR(IP,j)=XR;
            KI(IP,j)=XI;
        end
    end
    for j=1:NR
            XR=UR(IE,j)*ZR-UI(IE,j)*ZI;
            XI=UR(IE,j)*ZI+UI(IE,j)*ZR;
            UR(IE,j)=UI(IP,j);
            UI(IE,j)=UI(IP,j);
            UR(IP,j)=XR;
            UI(IP,j)=XI;
    end 
    %Multiply rows by multiple of pivot row
    for i=1:NP
        if (i==IP)
            continue;
        end
        if (IP~=NP)
            for j=K:NP
                KR(i,j)=KR(i,j)-KR(IP,j)*KR(i,IP)+KI(IP,j)*KI(i,IP);
                KI(i,j)=KI(i,j)-KR(IP,j)*KI(i,IP)-KI(IP,j)*KR(i,IP);
            end
        end
        for j=1:NR
            UR(i,j)=UR(i,j)-UR(IP,j)*KR(i,IP)+UI(IP,j)*KI(i,IP);
            UI(i,j)=UI(i,j)-UR(IP,j)*KI(i,IP)-UI(IP,j)*KR(i,IP);
        end    
    end
end
return
end


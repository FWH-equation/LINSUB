function Dswx()
global Lam ;
global NP NR XR XI VU VD Apd Apu;
%Lift,Moment,And Wakes
for i=1:NP
    Psi=pi*(i-1)/NP;
    Z=0.5*(1.0-cos(Psi));
    XR(1,i)=-1.0;
    XI(1,i)=0.0;
    XR(2,i)=-Z;
    XI(2,i)=0.0;
    X=Z*Lam;
    XR(3,i)=Lam*sin(X);
    XI(3,i)=-Lam*cos(X);
    if (NR~=5)
        continue;
    end
    %Outgoing Waves
    X=Z*Apu;
    XR(4,i)=-VU*cos(X);
    XI(4,i)=VU*sin(X);
    X=Z*Apd;
    XR(5,i)=-VD*cos(X);
    XI(5,i)=VD*sin(X);
end
return
end


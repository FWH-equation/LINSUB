function Dswu()
global Lam ;
global IR NP UR UI NR ;
global Apu Apd Anu And;

%Wave Properties
if(NR==5)
    IR=0;
    Wave();
    Wavu=Anu/(Lam+Apu);
    Wavd=And/(Lam+Apd);
end
%Bending Torsion And Wakes
for i=1:NP
    Epsil=pi*(2*i-1)/(2*NP);
    Z=0.5*(1.0-cos(Epsil));
    X=Z*Lam;
    UR(i,1)=1.0;
    UI(i,1)=0.0;
    UR(i,2)=1.0;
    UI(i,2)=X;
    UR(i,3)=-cos(X);
    UI(i,3)=sin(X);
    if (NR~=5)
        continue;
    end
    %Upstream and downstream waves
    X=Z*Apu;
    UR(i,4)=Wavu*cos(X);
    UI(i,4)=Wavu*sin(X);
    X=Z*Apd;
    UR(i,5)=Wavd*cos(X);
    UI(i,5)=Wavd*sin(X);
end
return;
end


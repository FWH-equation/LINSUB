function Dswk()
%Calculation of kernel matrix
%constants for vortex sheet calculation
global SC Lam Phase Cosst Sinst Mach2 B B2;
global XU Apu Apd Pur Pdr Pui;
global IR IW KR KI NP;
X=Lam*SC*Cosst;
Y=Lam*SC*Sinst+Phase;
Vort=0.5*Lam*sinh(X)/(cosh(X)-cos(Y));
%constants for log singularity correction
Mach4=Mach2^2;
Mach6=Mach2*Mach4;
B4=B2^2;
B6=B2*B4;
A1=1.0-0.5*Mach2/B2;
A2=1.0-0.5/B2+0.25*Mach2/B4;
A3=0.5*(1.0-1.0/B2+Mach2/(6.0*B4)+1.0/(3.0*B4)-0.375*Mach4/B6+Mach6/(6.0*B6));
%Matching and vortex points
NNP = NP; 
ZE=ones(NNP,1);
ZP=ones(NNP,1);
for i=1:NP
    Epsil=pi*(2*i-1)/(2*NP);
    ZE(i)=0.5*(1.0-cos(Epsil));
    Psi=pi*(i-1)/NP;
    ZP(i)=0.5*(1.0-cos(Psi));
end
%zero counts and arrays
IR=0;
Icount=0;
IW=0;
NP2=NP^2;
Icheck=ones(NNP,NNP);
for i=1:NP
    for j=1:NP
        Icheck(i,j)=0;
        KR(i,j)=0.0;
        KI(i,j)=0.0;
    end
end
while (1)
%assemble matrix
%i(=M+1 in paper) gives vortex position
%j(=L=1 in paper) gives matching point
Wave();
if (IW==1)
    return;
end
for i=1:NP
    for j=1:NP
        if (Icheck(i,j)==1)
            continue;
        end
        Pos=ZE(i)-ZP(j);
        if (Pos>0.0)
            %Downstream point
            XP=exp(-XU*Pos);
            YP=Apd*Pos;
            QR=XP*cos(YP);
            QI=XP*sin(YP);
            Termr=(Pdr*QR-Pui*QI)/SC;
            Termi=(Pdr*QI+Pui*QR)/SC;
        else
            %upstream point
            XP=exp(XU*Pos);
            YP=Apu*Pos;
            QR=XP*cos(YP);
            QI=XP*sin(YP);
            Termr=(Pur*QR-Pui*QI)/SC;
            Termi=(Pur*QI+Pui*QR)/SC;
        end
        %Add to matrix
        KR(i,j)=KR(i,j)+Termr;
        KI(i,j)=KI(i,j)+Termi;
        %Check convergence of series
        X=Termr^2+Termi^2;
        Y=KR(i,j)^2+KI(i,j)^2;
        if ((X/Y)>1.0*10E-10)
            continue;
        end
        Icheck(i,j)=1;
        Icount=Icount+1;
        %Correct for log singularity (last time through)
        sum=0.0;
        Epsil=pi*(2*i-1)/(2*NP);
        Psi=pi*(j-1)/NP;
        Npm1=NP-1;
        for JR=1:Npm1
            Fjr=JR;
            sum=sum+cos(Fjr*Epsil)*cos(Fjr*Psi)/Fjr;
        end
        sum=2.0*sum+log(4.0*abs(Pos));
        sum=sum*Lam/(2.0*pi*B);
        Plam=Lam*Pos;
        Plam2=Plam^2;
        Plam3=Plam2*Plam;
        KR(i,j)=KR(i,j)+sum*(A1*Plam-A3*Plam3);
        KI(i,j)=KI(i,j)+sum*(1.0-A2*Plam2);
        %Add vorticity wave
        if(Pos<=0.0)
            continue;
        end
        KR(i,j)=KR(i,j)+Vort*cos(Plam);
        KI(i,j)=KI(i,j)-Vort*sin(Plam);
        fprintf (' %d, %d, %d, %010.5f, %010.5f, %d\r\n',i,j,IR,KR(i,j),KI(i,j),Icount);
    end
end
%Check for completion
if (Icount==NP2)
    return;
end
if (IR>0.0)
    IR=-IR;
else
    IR=-IR+1;
end
end
end



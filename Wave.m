function Wave()
global SC Lam Phase Cosst Sinst Mach2 BC2;
global XU Apu Apd Anu And Pur Pdr Pui VU VD;
global IR IW;
%Calculation of acoustic wave properties
Betah=(Phase-2.0*pi*IR)/SC;
Betah2=Betah^2;
A=Lam^2+Betah2+2.0*Lam*Betah*Sinst;
D=Mach2*(Lam+Betah*Sinst)*Cosst/BC2;
E=Betah2-Mach2*A;
if (E==0.0)
    fprintf('Resonance at IR=%d\r\n',IR);
    IW=1;
    return;
else
    F=sqrt(abs(E));
    FB=F/BC2;
    H=(Betah+Lam*Sinst)/(2.0*A);
    P=Betah*Lam*Cosst/(F*2.0*A);
    if(E>0.0)
        %Wave numbers,decaying case
        Apu=D*Cosst+Betah*Sinst;
        Apd=Apu;
        Anu=Betah*Cosst-D*Sinst;
        And=Anu;
        Pur=-Anu*H-FB*Sinst*P;
        Pdr=-Pur;
        Pui=Anu*P-FB*Sinst*H;
        XU=FB*Cosst;
        return;
    else
        %Wave numbers, propagating case
        Acui=D+FB;
        Acdi=D-FB;
        Apu=Acui*Cosst+Betah*Sinst;
        Apd=Acdi*Cosst+Betah*Sinst;
        Anu=Betah*Cosst-Acui*Sinst;
        And=Betah*Cosst-Acdi*Sinst;
        Pur=Anu*(P-H);
        Pdr=And*(P+H);
        Pui=0.0;
        XU=0.0;
        VU=(P-H)*(Lam+Apu)/SC;
        VD=(P+H)*(Lam+Apd)/SC;
        return;
    end
end
end


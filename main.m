%linearised subsonic unsteady flow in cascade
%initial set-up
% this code was written by Zhang Ji from Beihang University(BUAA). Cheng Long
% from BUAA debugged the code and verified the solution with standard LINSUB solution. 
% If you have any questions, please contact with Cheng Long. 
% Long.Cheng@buaa.edu.cn
% 2020.10.23
%%  
% input introduction: 
% 1. SC: Space to chord ratio. SC > 0.0; 
% 2. stagger angle(Degree): -90Deg.<Stag<90Deg.;
% 3. Mach: inflow Mach number, 0.0<= Mach <1.0; 
% 4. Lam: frequency parameter, reduced frequency, Lam = omega*c/U, where omega
%         is the angular frequency, radians/sec. ;
%         c is the blade chord; 
%         U is the inflow velocity; 
% 5. PHASE: inter-blade-phase between, Degree.   -180Drg. < Degree <=180Deg.; 
%%----------------------------------------------------------- 
% output: 
% 
% 
clear
clc
global SC Stag Mach Lam Phase Deg Cosst Sinst Mach2 B B2 BC BC2;
global IW NP UR UI NR CR CI ;
Deg=57.29578;
Prnopt=input('Print Presure Jump? Y or N \n','s');
NP=input('Please input the number of points on blade:\n');
while (NP>105) 
    %NP=input('Please re-input the number of points on blade(<15):\n');
end
fp = fopen('./AR1.txt','a');
fprintf(fp,'NP=%d\r\n',NP);
while(1)
%read and print data
SC=input('Space/Chord? \n');
Stag=input('Stagger angle,degrees? \n');
Mach=input('Mach number? \n');
while (Mach>=1.0) 
    Mach=input('Re-input Mach number(<1.0): \n');
end
Lam=input('Frequency Parameter? \n');
fprintf(fp,'Space/Chord=%.5f\r\nStagger angle=%.5f degrees\r\n',SC,Stag);
fprintf(fp,'Mach Number=%.5f\r\nFrequency Parameter=%.5f\r\n',Mach,Lam);
Stag=Stag/Deg;
%constants independent of phase
Cosst=cos(Stag);
Sinst=sin(Stag);
Mach2=Mach^2;
B2=1.0-Mach2;
B=sqrt(B2);
BC2=1.0-Mach2*Cosst^2;
BC=sqrt(BC2);
%calculate phase angles for resonance
Z=Mach*Lam*SC/B2;
Phase1=(Sinst*Mach-BC)*Z;
Phase2=(Sinst*Mach+BC)*Z;
X=Phase1*Deg;
Y=Phase2*Deg;
fprintf('Phase angles for resonance=%.5f,%.5f Degrees\r\n',X,Y);
fprintf(fp,'Phase angles for resonance=%.5f,%.5f Degrees\r\n',X,Y);
while(1)
%Read and print phase
Phase=input('Phase angle= degrees? \n');
fprintf(fp,'Phase angle=%.5f degrees\r\n',Phase);
Phase=Phase/Deg;
NR=3;
if(Phase>Phase1 && Phase<Phase2)
    NR=5;
end
%Matrix generation and algebra
Dswk();
if (IW~=1)
    Dswu();
    Dswx();
    Cmdiv();
    Cmprd();
    %Print results
    fprintf(fp,'                  Bending             Torsion             Wakes\r\n');
    Line={'Force','Moment','Wake','Up Wave','Down Wave'};
    Line = string(Line);
    for i=1:NR
        fprintf(fp,'%s ',Line(i));
        for j=1:3
            fprintf(fp,'%.5f,%.5f,',CR(i,j),CI(i,j));
        end
        fprintf(fp,'\r\n');
    end
    if (NR==5)
        fprintf(fp,'                  Wave up             Wave down\r\n');
        for i=1:NR
            fprintf(fp,'%s ',Line(i));
            for j=4:5
                fprintf(fp,'%.5f,%.5f,',CR(i,j),CI(i,j));
            end
            fprintf(fp,'\r\n');
         end
    end
    %Pring Pressure Jump
    if (Prnopt=='Y')
        fprintf(fp,'Pressure\r\n');
        AN=NP;
        AK=AN*2.0/pi;
        %for i=2:NP
        for i=1:NP
            Y=pi*(i-1)/AN;
            X(i)=(1.0-cos(Y))/2.0;
            AM=AK/sin(Y);
            for j=1:NR
                UR(i,j)=UR(i,j)*AM;
                UI(i,j)=UI(i,j)*AM;              
            end
        end
        for i=2:NP
            if ( i == NP )
                fprintf(fp,' %.5f\r\n',X(i));
            else
                fprintf(fp,' %.5f,',X(i));
            end
        end
        for i=2:NP
            for j=1:3
                if ( j == 3 )
                    fprintf(fp,'%.5f,%.5f;',UR(i,j),UI(i,j));
                else
                    fprintf(fp,'%.5f,%.5f,',UR(i,j),UI(i,j));
                end
            end
            fprintf(fp,'\r\n');
            if (NR~=5)
                continue;
            end
            for j=4:5
                fprintf(fp,'%.5f,%.5f;',UR(i,j),UI(i,j));
            end
            fprintf(fp,'\r\n');
        end
    end
end
%Other Cases
Repeat=input('Another Phase Angle?Y or N\n','s');
if(Repeat~='Y')
    break;
end
end
Repeat=input('Another Whole Set Of Data?Y or N\n','s');
if(Repeat~='Y')
    break;
end
end
%Empty Printer Buffer
fprintf(fp,'\r\n');
fclose(fp); 
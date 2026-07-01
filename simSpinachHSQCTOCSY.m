  function [spectrum,XAxisInd,XAxisFID]=simSpinachHSQCTOCSY(spin_system,w_1,t_90,t_180,npFID,npFIDFT,npInd,npIndFT,swFID_Hz,swInd_Hz,OffsetHzFID,OffsetHzInd,SFO1,SFO2,NUC1)
  %Function to run 1H,31P HSQC-TOCSY experiment within MR_SHIFT using SPINACH         
  
  %Default takes the dwell time of direct and indirect dimensions from experiment
           %but it is possible to enter any arbitrary value, therefore
           %recalculate the axis here:

                dwellFID=1/swFID_Hz; 
                dwell=dwellFID;

                %Calculate indirect dimension from exp SW:
                dwellInd=1/swInd_Hz;
                pointsF1=npInd;%Divide by 2 due to real/imaginary
                CT=dwellInd*pointsF1/2;

                % pointsF1
                % %Alternatively, insert parameters manually:
                % CT=0.074844/2;%Constant time F1 evolution in sec (calculated from d2 in the used pulse sequence)
                % pointsF1=npInd/2;%Points in indirect dimension
                % dwellInd=CT/pointsF1;
                
                %Additional parameters for pulse sequence (taken from 
                np_trans=8;%number of XY16 while magnetization on xy
                np=11;%number of XY16 while magnetization on z
                Delta=201.4e-6*2;

   
                %Define the axis of the direct dimension in ppm and with the offset:
        SW_FID_ppm=swFID_Hz/SFO1; 
        ppmpPT=SW_FID_ppm/npFIDFT;
        XAxisFID=-SW_FID_ppm/2:ppmpPT:SW_FID_ppm/2-ppmpPT;
        XAxisFID=XAxisFID+OffsetHzFID/SFO1; 

                %Define the Axis of the indirect dimension in ppm and with the offset
        SW_Ind_ppm=swInd_Hz/SFO2; 
        ppmpPT=SW_Ind_ppm/npIndFT;
        XAxisInd=-SW_Ind_ppm/2:ppmpPT:SW_Ind_ppm/2-ppmpPT;
        XAxisInd=XAxisInd+OffsetHzInd/SFO2;



R=relaxation(spin_system);
H=hamiltonian(spin_system);
Lx_FH=operator(spin_system,'Lx','1H');
Ly_FH=operator(spin_system,'Ly','1H');
Lx_FP=operator(spin_system,'Lx','31P');
Ly_FP=operator(spin_system,'Ly','31P');

% Transmitter offsets
parameters.spins={'31P','1H'};
parameters.offset=[OffsetHzInd OffsetHzFID];
H=frqoffset(spin_system,H,parameters);

%Pulse Hamiltonian
H_1x=Lx_FH*w_1;
H_1y=Ly_FH*w_1;
H_1xP=Lx_FP*w_1;
H_1yP=Ly_FP*w_1;

H0=hamiltonian(assume(spin_system,'labframe'),'left');
rho0=equilibrium(spin_system,H0);
rho=rho0;
%% %PULSE SEQUENCE START
% Pulse on Proton
rho=step(spin_system,H_1x,rho,t_90);
rho=coherence(spin_system,rho,{{'1H', [+1,-1]}, {'31P',[0]}});

%Initial XY16
for n=1:np_trans
rho=XY16_IS_Evol(spin_system,rho,R,H,H_1yP,H_1y,H_1xP,H_1x,Delta,t_180);%XY16 with pulses on 1H and 31P
end

rho=coherence(spin_system,rho,{{'1H', [+1,-1]}, {'31P',[0]}});
%% Bring proton to z and P to xy
% Pulse on Proton 
rho=step(spin_system,Ly_FH,rho,pi/2);
%Z-spoil
rho=homospoil(spin_system,rho,'keep');
% Pulse on Phosphorous with coherence selection
rho=step(spin_system,Lx_FP,rho,pi/2)-step(spin_system,Lx_FP,rho,-pi/2);
rho=coherence(spin_system,rho,{{'1H', 0}, {'31P',[+1,-1]}});
%% t1 evolution
% Get the time grid for the CT period
if pointsF1 == 1
    t1_grid = 0;
else
    % t1_grid=linspace(0,CT,pointsF1);
    t1_grid = (0:(pointsF1-1)) * dwellInd;

end
% Loop over the value of t1_grid
for n=1:pointsF1
      % Run the delay
    rho_current=evolution(spin_system,H+1i*R,[],rho,(CT-t1_grid(n))/4,1,'final');
    % Pi on X
    rho_current=step(spin_system,H_1xP,rho_current,t_180);
     % Run the delay
    rho_current=evolution(spin_system,H+1i*R,[],rho_current,(CT+t1_grid(n))/4,1,'final');
        % Pi on Proton
    rho_current=step(spin_system,H_1x,rho_current,t_180);
         % Run the delay
    rho_current=evolution(spin_system,H+1i*R,[],rho_current,(CT+t1_grid(n))/4,1,'final');
    % Pi on X
    rho_current=step(spin_system,H_1xP,rho_current,t_180);%-step(spin_system,H_1yP,rho_current,t_180); %Dont see any difference when adding 02 phasecycle here
      % Run the delay
    rho_current=evolution(spin_system,H+1i*R,[],rho_current,(CT-t1_grid(n))/4,1,'final');
    % Assign the stack element
    rho_stack(:,n)=rho_current;
end
% Coherence selection
rho_stack_pos=coherence(spin_system,rho_stack,{{'1H', 0},{'31P',[+1]}});
rho_stack_neg=coherence(spin_system,rho_stack,{{'1H', 0},{'31P',[-1]}});
%% Start Second INEPT
% Pulses on F1 and F2:
rho_stack_pos=step(spin_system,Lx_FP,rho_stack_pos,pi/2);
rho_stack_neg=step(spin_system,Lx_FP,rho_stack_neg,pi/2);
%
rho_stack_pos=step(spin_system,H_1x+1i*R+H,rho_stack_pos,t_90);% 
rho_stack_neg=step(spin_system,H_1x+1i*R+H,rho_stack_neg,t_90);% 
% 
rho_stack_pos=coherence(spin_system,rho_stack_pos,{{'31P', 0}, {'1H',[+1,-1]}});
rho_stack_neg=coherence(spin_system,rho_stack_neg,{{'31P', 0},{'1H',[+1,-1]}});
%% INEPT CPMG IN xy plane
for n=1:np_trans
rho_stack_pos=XY16_IS_Evol(spin_system,rho_stack_pos,R,H,H_1yP,H_1y,H_1xP,H_1x,Delta,t_180);%XY16 with pulses on 1H and 31P
rho_stack_neg=XY16_IS_Evol(spin_system,rho_stack_neg,R,H,H_1yP,H_1y,H_1xP,H_1x,Delta,t_180);%XY16 with pulses on 1H and 31P
end
% 90 pulse to bring magnetization into z
rho_stack_pos=step(spin_system,H_1y+1i*R+H,rho_stack_pos,t_90);% 
rho_stack_neg=step(spin_system,H_1y+1i*R+H,rho_stack_neg,t_90);% 
%% CPMG along z for TOCSY evolution
for n=1:np
rho_stack_pos=XY16Evol(spin_system,rho_stack_pos,R,H,H_1y,H_1x,Delta,t_180);%XY16 with pulses on 1H only
rho_stack_neg=XY16Evol(spin_system,rho_stack_neg,R,H,H_1y,H_1x,Delta,t_180);%XY16 with pulses on 1H only
end
%% Remove multispin terms
%Destroy all undesired ZQ terms (Done experimentally with a composite plus
%z-spoil)
rho_stack_pos=homospoil(spin_system,rho_stack_pos,'destroy');
rho_stack_neg=homospoil(spin_system,rho_stack_neg,'destroy');
%% Bring to xy and acquire
coil=state(spin_system,'L+','1H');

rhoF_stack_pos=step(spin_system,H_1y,rho_stack_pos,t_90);% 
rhoF_stack_neg=step(spin_system,H_1y,rho_stack_neg,t_90);% 

L=H+1i*R;
% Decoupling of P spin
[L,rhoF_stack_pos]=decouple(spin_system,L,rhoF_stack_pos,{'31P'});
[L,rhoF_stack_neg]=decouple(spin_system,L,rhoF_stack_neg,{'31P'});
% Detect
fidpos=evolution(spin_system,L,coil,rhoF_stack_pos,dwell,npFID-1,'observable');
fidneg=evolution(spin_system,L,coil,rhoF_stack_neg,dwell,npFID-1,'observable');
 % Apodisation
    fidpos=apodisation(spin_system,fidpos,{{'sqcos'},{'sqcos'}});
    fidneg=apodisation(spin_system,fidneg,{{'sqcos'},{'sqcos'}});
    % F2 Fourier transform
    f1_pos=fftshift(fft(fidpos,npFIDFT,1),1);
    f1_neg=fftshift(fft(fidneg,npFIDFT,1),1);
    
    % Form States signal
    fid=(f1_neg)+conj(f1_pos);
        [M,Ind] = max((fid(:,1)));
    theta=atan2(imag(fid(Ind,1)),real(fid(Ind,1)))%0th order phase correction
        % F1 Fourier transform
    spectrum=fftshift(fft(fid,npIndFT,2),2);
    parameters.sweep=[1/dwellInd 1/dwell];
parameters.npoints=[pointsF1 npFID];
parameters.zerofill=[1024 npFID];
% parameters.axis_units='ppm';
%     figure(100)
%             plot_2d(spin_system,-real(spectrum),parameters,20,[0.01 1.0 0.01 1.0],2,256,6,'both');

function [rho]=XY16_IS_Evol(spin_system,rho,R,H,H_1yP,H_1y,H_1xP,H_1x,Delta,t_180)
 %xyxy
rho=step(spin_system,1i*R+H,rho,Delta/2);% 
rho=step(spin_system,H_1x+H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 

       %yxyx
rho=step(spin_system,H_1y+H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 


 %m(xyxy)
rho=step(spin_system,-H_1x-H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y-H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x-H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y-H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 

       %m(yxyx)
rho=step(spin_system,-H_1y-H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x-H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y-H_1yP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x-H_1xP+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta/2);% 
end



function [rho]=XY16Evol(spin_system,rho,R,H,H_1y,H_1x,Delta,t_180)

 %xyxy
rho=step(spin_system,1i*R+H,rho,Delta/2);% 
rho=step(spin_system,H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 

       %yxyx
rho=step(spin_system,H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 


 %m(xyxy)
rho=step(spin_system,-H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 

       %m(yxyx)
rho=step(spin_system,-H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1y+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta);% 
rho=step(spin_system,-H_1x+1i*R+H,rho,t_180);% 
rho=step(spin_system,1i*R+H,rho,Delta/2);% 
end
            end
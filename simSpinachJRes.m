        function [spectrum,XAxisInd,XAxisFID]=simSpinachJRes(spin_system,w_1,t_180,npFID,npFIDFT,npInd,npIndFT,swFID_Hz,swInd_Hz,OffsetHzFID,SFO1,NUC1)
    %Function to run JRes experiment within MR_SHIFT using SPINACH         

    %Default takes the dwell time of direct and indirect dimensions from experiment
           %but it is possible to enter any arbitrary value, therefore
           %recalculate the axis here:

                dwellInd=1/swInd_Hz;
                % dwellInd=1e-3;% set manual value for indirect spacing
                dwellFID=1/swFID_Hz; 

   
                %Define the axis of the indirect dimension without offset and in Hz:
                XAxisInd=-1/(2*dwellInd)+0.5/(dwellInd*npIndFT):(1/(dwellInd*npIndFT)):1/(2*dwellInd)-0.5/(dwellInd*npIndFT);         
        
                %Define the Axis of the direct dimension in ppm and with the offset
                SW_FID_ppm=swFID_Hz/SFO1; 
                ppmpPT=SW_FID_ppm/npFIDFT;
                XAxisFID=-SW_FID_ppm/2:ppmpPT:SW_FID_ppm/2-ppmpPT;
                XAxisFID=XAxisFID+OffsetHzFID/SFO1;



                R=relaxation(spin_system);
                H=hamiltonian(spin_system);

                parameters.spins={NUC1};
                parameters.offset=[OffsetHzFID];
                H=frqoffset(spin_system,H,parameters);

                Lx_FH=operator(spin_system,'Lx',NUC1);
                Ly_FH=operator(spin_system,'Ly',NUC1);
                %Pulse Hamiltonian
                H_1x=Lx_FH*w_1;
                H_1y=Ly_FH*w_1;

                H0=hamiltonian(assume(spin_system,'labframe'),'left');
                rho0=equilibrium(spin_system,H0);
                
                rho=rho0;
                % 90 Pulse on Proton
                rho=step(spin_system,Ly_FH,rho,pi/2);
                % Loop over the value of t1_grid
                for n=1:npInd
                      tau=dwellInd*(n-1);
                          % Run the delay
                    rho_current=evolution(spin_system,H,[],rho,tau/2,1,'final');
                
                        % Pi on Proton
                    rho_current=step(spin_system,H_1x,rho_current,t_180);
                
                        % Run the delay
                    rho_current=evolution(spin_system,H,[],rho_current,tau/2,1,'final');
                
                    rho_stack(:,n)=rho_current;
                end
                % Acquire
                coil=state(spin_system,'L+',NUC1);
                
                L=H+1i*R;
                % Detect
                fid=evolution(spin_system,L,coil,rho_stack,dwellFID,npFID-1,'observable');
                % Apodisation
                % fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}});
                fid=apodisation(spin_system,fid,{{'sqcos'},{'sqsin'}});

                t1=(dwellFID:dwellFID:dwellFID*npFID)';
                % F1 Fourier transform
                fid=fftshift(fft(fid,npIndFT,2),2);
                alpha = -1; % shear strength
                fid=fid.*exp(1i*2*pi*alpha.*(t1*XAxisInd));

                % F2 Fourier transform
                spectrum=fftshift(fft(fid,npFIDFT,1),1);
%Symmetrize indirect dimension around zero:
        spectrum=abs(spectrum);
        spectrum = min(spectrum, fliplr(spectrum));
            end
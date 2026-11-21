%electric_guitar.m
clear all
close all
clc

%--------------------------Riff-parameters---------------
%numerating the desired riff/scale to be played 
%Each note is specified by:
%tnote = onset time of note (s)
%anote = relative pluck amplitude
%inote = string index: 1=low E, 2=A, 3=D, 4=G, 5=B, 6=high e
%fnote = fret index: 0=open string, 1=first fret, etc.

%sample data to play a series of power chords
tnote = [ ...
    0.00 0.50 1.00 1.50 2.00 2.50 ...
    3.00 3.00 3.50 3.50 4.00 4.00 ...
    4.50 5.00 5.50 6.00 ...
    6.50 6.50 7.00 7.00];
anote = [ ...
    1.00 1.00 1.00 1.00 1.00 1.00 ...
    0.90 0.90 0.90 0.90 0.90 0.90 ...
    1.00 1.00 1.00 1.00 ...
    1.00 1.00 1.00 1.00];
inote = [ ...
    4 4 4 4 4 4 ...
    4 3 4 3 4 3 ...
    4 4 4 4 ...
    2 3 2 3];
fnote = [ ...
    0 3 5 0 3 6 ...
    5 5 3 3 0 0 ...
    3 5 3 0 ...
    0 2 0 2];

nmax = length(tnote); %number of notes
%------------------------------------------------------


%---------------String-Parameters----------------------
%All guitar strings and their relevant dimensional data and full length fundamental frequency - uniform density is assumed when using gauge (diameter) to calculate mass per unit length for each string
%Everything is normalized to metric units, from low E to high e (strings left to right on the guitar)
%tension values were found online and converted to Newtons as a reference for the calculated values

stringDensity = 8.05*10^3; % density of each string in kg/m^3 - found online
L_0 = 0.62865; %scale length (full length for each string)

%E - string 1
gaugeE = 1.1684*10^(-3); %m
frequencyE = 82.41; %Hz
M_E = stringDensity*pi*((gaugeE/2)^2); %mass per unit length - kg/m
tensionE = M_E*(2*L_0*frequencyE)^2; %~73.218 N

%A - string 2
gaugeA = 9.652*10^(-4); %m
frequencyA = 110; %Hz
M_A = stringDensity*pi*((gaugeA/2)^2); %mass per unit length - kg/m
tensionA = M_A*(2*L_0*frequencyA)^2; %~90.388 N

%D - string 3
gaugeD = 7.112*10^(-4); %m
frequencyD = 146.83; %Hz
M_D = stringDensity*pi*((gaugeD/2)^2); %mass per unit length - kg/m
tensionD = M_D*(2*L_0*frequencyD)^2; %~89.231 N

%G - string 4
gaugeG = 4.318*10^(-4); %m
frequencyG = 196; %Hz
M_G = stringDensity*pi*((gaugeG/2)^2); %mass per unit length - kg/m
tensionG = M_G*(2*L_0*frequencyG)^2; %~69.437 N

%B - string 5
gaugeB = 3.302*10^(-4); %m
frequencyB = 246.94; %Hz
M_B = stringDensity*pi*((gaugeB/2)^2); %mass per unit length - kg/m
tensionB = M_B*(2*L_0*frequencyB)^2; %~64.41 N

%e - string 6
gauge_e = 2.54*10^(-4); %m
frequency_e = 329.63; %Hz
M_e = stringDensity*pi*((gauge_e/2)^2); %mass per unit length - kg/m
tension_e = M_e*(2*L_0*frequency_e)^2; %~67.88 N

%storing all relevant variables in matrices:
M = [M_E M_A M_D M_G M_B M_e];
T = [tensionE tensionA tensionD tensionG tensionB tension_e];
f1 = [frequencyE frequencyA frequencyD frequencyG frequencyB frequency_e]; %matrix of fundamental frequencies for each string 

nstrings = 6; %number of strings
%--------------------------------------------------------


%--------------Fret-placements---------------------------
%Frets are positioned so that shortening the active vibrating segment by one fret raises the pitch by one half-step - these successively shorter distances (going down the neck of the guitar) are calculated by the frequency ratio between notes separated by n half steps

nfretmax = 22; %number of frets included
xfret = zeros(1,nfretmax+1); %xfret(1)=0 for open string
for n=1:nfretmax
    xfret(n+1) = L_0*(1 - 2^(-n/12));
end
%---------------------------------------------------------


%--------Initialize-Parameters---------------------------
%number of considered points on each string and spatial step (121 spatial points chosen in the same ballpark as the notes)
J = 121;
dx = L_0/(J-1);

%convert fret positions to nearest grid indices
jfret = ones(1,nfretmax+1);
for n=1:nfretmax
    jfret(n+1) = round(xfret(n+1)/dx) + 1;
end

%fill in matrix parameters for artificial vibration damping and max time step
R = zeros(1,6);
tau = zeros(1,6);
dtmax = zeros(1,6);

%decay time(seconds)
for n=1:6
    tau(n) = 1.6*(110/f1(n))^(0.30);
    R(n) = (2*M(n)*L_0^2)/(tau(n)*pi^2);
    dtmax(n) = -R(n)/T(n) + sqrt((R(n)/T(n))^2 + dx^2/(T(n)/M(n)));
end

%bounding dt to ensure stability for all string simulations and setting nskip such that dt*nskip=1/8192 - to ensure smooth playback since matlab samples 8192 times per second
dtmaxmin = 0.5*min(dtmax);
nskip = ceil(1/(8192*dtmaxmin));
dt = 1/(8192*nskip);

%total simulation time
ringtime = 3.0; %extra ringing time after final note
tmax = tnote(nmax) + ringtime;
clockmax = ceil(tmax/dt); %set clockmax for main simulation loop

%initialize position and velocity arrays:
H = zeros(6, J);
V = zeros(6, J);

%matrix to track the active fret for each string i: 0 means open string
fretNow = zeros(1,6);

%matrix to track the fixed points for each string given the fret being held (from said fret to the top of the string):
jfixed = ones(1,6);
%-------------------------------------------------------


%--------Pickup/electronics-simulation-------------------
%pickup location along the string:
%bridge pickup option:
pickupFrac = 0.92;

%neck pickup option:
%pickupFrac = 0.78;


jp = round(pickupFrac*(J-1)) + 1; %point the center of the pickup is located at

%finite-width pickup weights:
pickupWeights = [0.20 0.60 0.20];

%preamp gain and effects:
gain = 12;
distDrive = 2.2;
toneAlpha = 0.20;
%-------------------------------------------------------


%-----------Sound-storage--------------------------------
count = 0; %initialize count of samples of recorded soundwave
S_clean = zeros(1,ceil(clockmax/nskip)); %array to record clean pickup signal (pre "electronic modification")
tsave = zeros(1,ceil(clockmax/nskip));   %array for times of samples

%list of interior points
j = 2:(J-1);

%initialize note counter
n = 1;
%-------------------------------------------------------


%--------------------Main-loop---------------------------
%taken almost exactly from the class notes for simplicity

for clock=1:clockmax
    t = clock*dt;

    while((n<=nmax) && tnote(n)<=t)
        i = inote(n);       %string number for this note
        f = fnote(n);       %fret number for this note
        a = anote(n);       %pluck amplitude for this note

        %record which fret is being held on this string
        fretNow(i) = f;

        %determine the point to which the string must be fixed at zero depending on which fret is being held
        if(f==0)
            jfixed(i) = 1; %open string
        else
            jfixed(i) = jfret(f+1); % will set the points to be fixed up until the spatial coordinate of the fret in question
        end

        %clear prior motion on that string before plucking it again
        H(i,:) = zeros(1,J);
        V(i,:) = zeros(1,J);

        %set and hold all points on the fixed fret side to zero
        H(i,1:jfixed(i)) = 0;
        V(i,1:jfixed(i)) = 0;

        %define pluck position within the active vibrating segment - the triangular pluck shape is defined only on the active segment from the fret to the bridge, and its calculated using the gradient from the pluck height to the fixed zero endpoints
        pluckFrac = 0.18; %distance from bridge as fraction of active length - change to do analysis on pluck position

        xleft = (jfixed(i)-1)*dx;     %location of fret boundary
        Lactive = L_0 - xleft;        %active vibrating length
        xp = xleft + (1-pluckFrac)*Lactive; %pluck location in physical coordinates
        Hp = a; %pluck amplitude

        %sets the string position to mirror a string being plucked - has a peak at the point of pluck and straight lines to the two endpoints - calculated using the "slope" from the peak point to each endpoint
        for jj=jfixed(i):J
            x = (jj-1)*dx;
            if(x<xp)
                H(i,jj) = Hp*(x-xleft)/(xp-xleft);
            else
                H(i,jj) = Hp*(L_0-x)/(L_0-xp);
            end
        end

        %re-enforce zero on the fixed fret side
        H(i,1:jfixed(i)) = 0;
        V(i,1:jfixed(i)) = 0;

        n = n+1; %consider next note
    end

    %main updating loop for position and velocity based on wave equation numerical solution. 
    for i=1:6
        %active interior points are only those beyond the fixed fret side
        jactive = max(jfixed(i)+1,2):(J-1);

        %only runs if some part of the string is able to vibrate
        if(~isempty(jactive))
            V(i,jactive) = V(i,jactive) ...
                +(dt/dx^2)*(T(i)/M(i))*(H(i,jactive+1)-2*H(i,jactive)+H(i,jactive-1)) ...
                +(dt/dx^2)*(R(i)/M(i))*(V(i,jactive+1)-2*V(i,jactive)+V(i,jactive-1));

            H(i,jactive) = H(i,jactive) + dt*V(i,jactive);
        end

        %keep the fixed fret side at zero for as long as the fret is held
        H(i,1:jfixed(i)) = 0;
        V(i,1:jfixed(i)) = 0;

        %bridge endpoint is fixed
        H(i,J) = 0;
        V(i,J) = 0;
    end

    %simulates the pickup - converts vibrations into electrical signals with the same frequency
    if(mod(clock,nskip)==0) %only record every nskip timesteps 
        count = count+1;

        pickupNow = 0; 

        %reads the velocity from the designated "pickup point" and the two points to either side. averages the velocity at these points with predesignated weights to simulate the pickup's magnetic field being influenced by the vibration of the ferromagnetic string over a small interval of the string
        for i=1:6
            j1 = max(jp-1,2);
            j2 = jp;
            j3 = min(jp+1,J-1);

            pickupNow = pickupNow ...
                + pickupWeights(1)*V(i,j1) ...
                + pickupWeights(2)*V(i,j2) ...
                + pickupWeights(3)*V(i,j3);

            %alternative simpler single-point pickup:
            %pickupNow = pickupNow + V(i,jp);

            %alternative displacement-based pickup:
            %pickupNow = pickupNow + H(i,jp);
        end

        %save clean soundwave data before simulated "electronic modification" (amplification, distortion, etc)
        S_clean(count) = pickupNow;
        tsave(count) = t;
    end
end
%-------------------------------------------------------


%--------Amplification/tone/distortion------------------
%preamp gain
S_clean_raw = S_clean(1:count);
S_amp = gain*S_clean_raw;

%simple one-pole, low pass tone filter
S_lpf = zeros(1,count);
S_lpf(1) = S_amp(1);
for k=2:count
    S_lpf(k) = toneAlpha*S_amp(k) + (1-toneAlpha)*S_lpf(k-1);
end

%soft clipping distortion applied after LPF - this is the final raw output
S_final = tanh(distDrive*S_lpf);

%normalize final soundwave - scales down wave so amplitude is always between 1 and -1. Changes nothing but overall volume (doesnt change frequency) and makes later data graphing and storage easier
if(max(abs(S_final))>0)
    S_final = S_final/max(abs(S_final));
end

%make separately normalized copies for similar scale fourier transform plotting
S_clean_plot = S_clean_raw;
if(max(abs(S_clean_plot))>0)
    S_clean_plot = S_clean_plot/max(abs(S_clean_plot));
end

S_lpf_plot = S_lpf;
if(max(abs(S_lpf_plot))>0)
    S_lpf_plot = S_lpf_plot/max(abs(S_lpf_plot));
end

S_final_plot = S_final;
if(max(abs(S_final_plot))>0)
    S_final_plot = S_final_plot/max(abs(S_final_plot));
end
%-------------------------------------------------------


%-------Plot and play soundwave--------------------------
%plays the normalized final soundwave
soundsc(S_final,8192)

%alternative function to play soundwave - hard clips the soundwave at 1 and -1
%sound(S_final,8192)

%takes the clean pickup data to graph
figure(1)
plot(tsave(1:count),S_clean_plot)
xlabel('Time (s)')
ylabel('Amplitude')
title('Clean pickup signal (normalized for plotting)')

%graphs post-LPF data
figure(2)
plot(tsave(1:count),S_lpf_plot)
xlabel('Time (s)')
ylabel('Amplitude')
title('Post LPF signal (normalized for plotting)')

%graphs modified, final sound data
figure(3)
plot(tsave(1:count),S_final_plot)
xlabel('Time (s)')
ylabel('Amplitude')
title('Final electric guitar output (normalized for plotting)')

%graphs a fourier transform of the clean signal, post LPF signal, and final signal
figure(4)
Nfft = 2^nextpow2(count); %ensures number of data points is a power of 2 - more efficient for fft
fplot = (0:Nfft-1)*(8192/Nfft);

Yclean = abs(fft(S_clean_plot,Nfft));
Ylpf = abs(fft(S_lpf_plot,Nfft));
Yfinal = abs(fft(S_final_plot,Nfft));

plot(fplot(1:floor(Nfft/2)),Yclean(1:floor(Nfft/2)),'LineWidth',1.0)
hold on
plot(fplot(1:floor(Nfft/2)),Ylpf(1:floor(Nfft/2)),'LineWidth',1.0)
plot(fplot(1:floor(Nfft/2)),Yfinal(1:floor(Nfft/2)),'LineWidth',1.0)
hold off
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum: clean pickup vs post LPF vs final output')
legend('Clean pickup','Post LPF','Final output')
xlim([0 4000])
grid on
%-------------------------------------------------------


%-------Save audio files--------------------------------
%save each stage as a .wav audio file
%audiowrite expects samples typically in the range [-1,1]

audiowrite('clean_pickup.wav', S_clean_plot.', 8192);
audiowrite('post_LPF.wav', S_lpf_plot.', 8192);
audiowrite('final.wav', S_final.', 8192);

%optional: also save all arrays in a MATLAB data file for later analysis
%save('electric_guitar_signals.mat','S_clean_raw','S_lpf','S_final','S_clean_plot','S_lpf_plot','S_final_plot','tsave');
%-------------------------------------------------------


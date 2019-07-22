clc
load chirp.mat;

sound(y, Fs);
inp = pwelch(y);
%inq = periodogram(y);
%inp=psd(spectrum.periodogram,y,'Fs',fs,'NFFT',length(y));
plot(y)
figure
plot(inp)
load dspwlets; % load wavelet coefficients and noisy signal 
Threshold = [3 2 1 0];
hsfw = dsp.SignalSource(y,64);

hdyadanalysis = dsp.DyadicAnalysisFilterBank( ... 
'CustomLowpassFilter', lod, ... 
'CustomHighpassFilter', hid, ... 
'NumLevels', 3);

hdelay1 = dsp.Delay(3*(length(lod)-1));
hdelay2 = dsp.Delay(length(lod)-1);
hdelay3 = dsp.Delay(7*(length(lod)-1));

hdyadsynthesis = dsp.DyadicSynthesisFilterBank( ... 
'CustomLowpassFilter', lor, ... 
'CustomHighpassFilter', hir, ... 
'NumLevels', 3);

hts = dsp.TimeScope('Name', 'Wavelet Denoising', ... 
'SampleRate', fs, ... 
'TimeSpan', 13, ... 
'NumInputPorts', 3, ... 
'LayoutDimensions',[3 1], ... 
'TimeAxisLabels', 'Bottom', ...   
'TimeSpanOverrunAction', 'Scroll'); 

pos = hts.Position;
hts.Position = [pos(1) pos(2)-(0.5*pos(4)) 0.9*pos(3) 2*pos(4)]; 
  
% Set properties for each display
hts.ActiveDisplay = 1;
hts.Title = 'Input Signal'; 
hts.ActiveDisplay = 2;
hts.Title = 'Denoised Signal'; 
hts.ActiveDisplay = 3;
hts.Title = 'Residual Signal';

for ii = 1:length(noisdopp)/64
    sig = step(hsfw);      % Input noisy signal 
    S = step(hdyadanalysis,sig);      % Dyadic analysis   
    
    % separate into four subbands 
    S1 = S(1:32);  S2 = S(33:48);  S3 = S(49:56);  S4 = S(57:64);   
    
    % Delay to compensate for the dyadic analysis filters 
    S1 = step(hdelay1,S1); 
    S2 = step(hdelay2,S2);
    S1 = dspDeadZone(S1, Threshold(1));
    S2 = dspDeadZone(S2, Threshold(2));
    S3 = dspDeadZone(S3, Threshold(3));
    S4 = dspDeadZone(S4, Threshold(4));
  
    % Dyadic synthesis (on concatenated subbands) 
    S = step(hdyadsynthesis,[S1; S2; S3; S4]); 
  
    sig_delay = step(hdelay3,sig);   % Delay to compensate for analysis/synthesis.
    Error = sig_delay - S;
    % Plot the results     
    step(hts,sig_delay, S, Error);
end
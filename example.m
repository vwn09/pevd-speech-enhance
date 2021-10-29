%% This test script demonstrates how the PEVD-based speech enhancement algorithm [1] can be used.

% Author:   Vincent W. Neo
% Date:     26-10-2021

% Copyright (C) 2017-2021 Vincent W. Neo
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

% Please acknowledge our paper if you use our code:
% [1] V. W. Neo, C. Evers, and P. A. Naylor
% Enhancement of noisy reverberant speech using polynomial matrix eigenvalue decomposition
% IEEE/ACM Trans. Audio, Speech and Lang. Process., vol. 28, 2021. doi:
% 10.1109/TASLP.2021.3120630


%% Example 1: Apply PEVD to multi-channel audio signal using default parameters

[x,fs] = audioread('audio/lecture2_babble0dB_noisy.wav');
y = pevd_enhance(x,fs);

audiowrite('audio/pevd_enhanced.wav',y(:,1)./max(abs(y(:))),fs);
% listen to the enhanced signal in the first channel by uncommenting:
% soundsc(x(:,1),fs);
% soundsc(y(:,1),fs);


%% Example 2: Non-default parameters of Example 1

[x,fs] = audioread('audio/lecture2_babble0dB_noisy.wav');

params.L = 1000; % 1000 iterations
params.mu = 10^-2; % trim factor
params.T = 800; % frame size of 800 samples
params.W = 800; % window on the space-time covariance
params.delta = 50; % off-diagonal column norm-2 target

y = pevd_enhance(x,fs,params);

% listen to the enhanced signal in the first channel by uncommenting:
% soundsc(x(:,1),fs);
% soundsc(y(:,1),fs);


%% Example 3: Measured recordings from corpus
% [2] Voicebox: https://github.com/ImperialCollegeLondon/sap-voicebox

[s,fs] = audioread('audio/clean.wav');  % sample speech from TIMIT
[h,~] = audioread('audio/lecture2_rir.wav');    % sample ACE Lecture Room 2 RIR
[v,~] = audioread('audio/lecture2_babble_noise.wav');   % ACE babble noise

SNR = 0;
signal_level_db = 10;
% ideally, should use v_activlev, v_addnoise from [4]
x = zeros(length(s)+size(h,1)-1,3);

for ii = 1:3
    si =  conv(s,h(:,1));

    % scale speech signal level, ideally v_activlev(si) from [2]
    Px = 10*log10((si'*si)/length(si));
    speech_linear_gain = 10^((- Px + signal_level_db)/20);
    si = speech_linear_gain * si;

    % scale noise signal level based on target SNR, ideally v_addnoise from [2]
    vi = v(1:length(si),ii);
    Pv = 10*log10((vi'*vi)/length(vi));
    noise_gain_db = 10^((- Pv + signal_level_db)/20);
    noise_linear_gain = 10.^((signal_level_db - Pv - SNR)/20);
    vi = noise_linear_gain * vi;

    % add signal and noise to generate microphone signal
    x(:,ii) = si+vi;
end

disp('Generated microphone signals.');

y = pevd_enhance(x,fs);

% Uncomment to listen to input and output
% soundsc(x(:,1),fs);
% soundsc(y(:,1),fs);


%% Example 4: Simulation (Pseduo-Code Only)
% Other Useful Simulation Tools - requires external toolboxes and databases.
% Please acknowledge the author(s) if you use their code and toolboxes
% [3] RIR Generator: https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
% [4] Diffuse Noise Generator: https://github.com/ehabets/ANF-Generator

%   Step 1: Read audio files including speech, RIR, noise signals as necessary
%   Step 2: Generate RIR using [3] if measured RIRs are not provided
%   Step 3: Scale the signal levels as necessary using v_activlev from [2] (See Example 3)
%   Step 4: Generate microphone signals by adding signal and noise using v_addnoise from [2]
%   Step 5: Perform speech enhancement using pevd_enhance
%   Step 6: Time align signals and evaluate performance

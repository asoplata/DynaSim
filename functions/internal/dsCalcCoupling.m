function data = dsCalcCoupling(data, varargin)
%CALCPOWER - Compute spectral analysis of DynaSim data
%
% Usage:
%   data = dsCalcPower(data,'option',value)
%
% Inputs:
%   - data: DynaSim data structure (see dsCheckData)
%   - options:
%     'variable'              : name of field containing data on which to
%                               calculate firing rates (default: *_spikes or
%                               first variable in data.labels)
%     'time_limits'           : [beg,end] (units of data.time)
%     'smooth_factor'         : number of samples for smoothing the spectrum (default: 5)
%                               - tip: set to 1 to avoid smoothing
%   - options for peak detection:
%     'min_peak_frequency'    : Hz, min frequency for peak detection (default: 2)
%     'max_peak_frequency'    : Hz, max frequency for peak detection (default: 150)
%     'peak_threshold_prctile': percentile for setting power threshold for peak
%                               detection (default: 95)
%     'peak_area_width'       : Hz, size of frequency bin (centered on peak)
%                               over which to calculate area under spectrum (default: 5)
%     'exclude_data_flag'     : whether to remove simulated data from result
%                               structure (default: 0)
%
% Outputs:
%   - data structure organization:
%     data.VARIABLE_Power_SUA.frequency: TODO
%     data.VARIABLE_Power_SUA.PeakArea: area under spectrum around peak (one value per cell)
%     data.VARIABLE_Power_SUA.PeakFreq: frequency of spectral power (one value per cell)
%     data.VARIABLE_Power_SUA.Pxx: spectral power
%   - if populations present, data also includes:
%     data.VARIABLE_Power_MUA.frequency: TODO
%     data.VARIABLE_Power_MUA.PeakArea: TODO
%     data.VARIABLE_Power_MUA.PeakFreq: TODO
%     data.VARIABLE_Power_MUA.Pxx: spectrum of the mean waveform
%       - population mean spectrum of the individual waveforms can be
%           calculated using "mean(data.VARIABLE_Power_MUA.Pxx,2)".
%   - Note:
%     - "VARIABLE" can be specified as the name of a variable listed in
%         data.labels, a cell array of string listing variable names, or as a
%         regular expression pattern for identifying variables to process. See
%         dsSelectVariables for more info on supported specifications.
%
% Examples:
%   s=[];
%   s.populations(1).name='E';
%   s.populations(1).equations='dv[2]/dt=@current+10; {iNa,iK}; v(0)=-65';
%   s.populations(2).name='I';
%   s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
%   data=dsSimulate(s,'tspan',[0 1000]);
%   data=dsCalcPower(data,'variable','v');
%   % Plot the spectrum of the E-cell average population voltage
%   figure; plot(data.E_v_Power_MUA.frequency,data.E_v_Power_MUA.Pxx);
%   xlabel('frequency (Hz)'); ylabel('power'); xlim([0 200]);
%
% See also: PlotPower, dsAnalyzeStudy, dsSimulate, dsCheckData, dsSelectVariables
%
% Author: Jason Sherfey, PhD <jssherfey@gmail.com>
% Copyright (C) 2016 Jason Sherfey, Boston University, USA

%% 1.0 Check inputs
options=dsCheckOptions(varargin,{...
  'variable',[],[],...
  'time_limits',[-inf inf],[],...
  'smooth_factor',5,[],... % number of samples for smoothing the spectrum
  'min_peak_frequency',1,[],... % Hz, min frequency for peak detection
  'max_peak_frequency',200,[],... % Hz, max frequency for peak detection
  'peak_threshold_prctile',95,[],... % percentile for setting power threshold for peak detection
  'peak_area_width',5,[],... % Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum
  'exclude_data_flag',0,{0,1},...
  'timeBandwidthProduct',[],[],... % time-bandwidth product for multi-taper method
  'output_suffix','',[],...
  'auto_gen_test_data_flag',0,{0,1},...
  'measure', 'mi',{'mi','esc','cfc'},...% Type of coupling measure
  'phase_freqs', [0.01:0.2:2.41],[],...% Hz, Frequencies to analyze for the phase of the "slow" or modulating signal
  'ampl_freqs', [4:1:14],[],...% Hz, Frequencies to analyze for the amplitude of the "fast" or carrier signal
  'plt', 'n',[],...% Don't use internal code to plot data
  'waitbar', 0,[],...% Don't print ongoing progress of significance analysis
  'width', 7,[],...% Width of Morlet wavelets to use for filtering, whatever?
  'nfft', 2500000,[],... % AES TODO Samples to use for each time/freq bin
  'num_shf', 0,[],...% Don't run any statistical significance analysis on coupling
},false);
  %% 'nfft', ceil(Fs/(diff(phase_freqs(1:2)))),[],...

data = dsCheckData(data, varargin{:});
% note: calling dsCheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use dsAnalyzeStudy to recursively call dsCalcPower on each data set
  data=dsAnalyzeStudy(data,@dsCalcCoupling,varargin{:});
  return;
end

% time parameters
time = data.time; % time vector
dt = time(2)-time(1); % time step
%% % ntime=length(time); % number of time points in full data set
t1=nearest(time,options.time_limits(1)); % index to first sample
t2=nearest(time,options.time_limits(2)); % index to last sample
nsamp=t2-t1+1; % number of samples for spectral estimate

% frequency parameters
Fs = fix(1/(dt/1000)); % effective sampling frequency
%% Fmin=options.min_peak_frequency; % min frequency for peak detection
%% Fmax=options.max_peak_frequency; % max frequency for peak detection
%% thresh_prctile=options.peak_threshold_prctile; % percentile for setting power threshold for peak detection
%% smooth_factor=options.smooth_factor; % number of samples to smooth spectrum
%% Fwin=options.peak_area_width; % size of frequency bin (centered on peak) for calculating area under spectrum
%% FreqRange=[max(Fmin,2/time(end)) Fmax]; % range to look for spectral peaks
%% NFFT=2^(nextpow2(nsamp-1)-1);%2); % <-- use higher resolution to capture STO freq variation
%% % WINDOW=2^(nextpow2(NFFT-1)-3);
%% % NOVERLAP=[]; % spectral parameters
%% NW = options.timeBandwidthProduct;


%% 2.0 set list of variables to process as cell array of strings
options.variable=dsSelectVariables(data(1).labels,options.variable, varargin{:});

%% 3.0 calculate power spectrum for each variable
if ~isfield(data,'results')
  data.results={};
end

warning off
for v=1:length(options.variable)
  % extract this data set
  var=options.variable{v};
  dat=data.(var);

  % determine how many cells are in this data set
  ncells=size(dat,2);

  %% % preallocation
  %% PeakFreq=nan(1,ncells);
  %% PeakArea=nan(1,ncells);

  % SUA spectra: loop over cells
  %% for i=1:ncells
  %%   % select data
  %%   X=detrend(dat(t1:t2,i)); % detrend the data
  %%   % calculate spectral estimate
  %%   if strcmp(reportUI,'matlab')
  %%     [tmpPxx,f] = pmtm(X, NW, NFFT, Fs); % calculate power
  %%   elseif exist('pwelch') == 2 % 'pwelch is in Octave's path
  %%     [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power in octave (pmtm is not implemented yet)
  %%   elseif exist('pwelch') ~= 2 % 'pwelch is not in Octave's path
  %%     try
  %%       pkg load signal; % trying to load octave forge 'signal' package before using pwelch function
  %%       fprintf('''pmtm'' function for spectral analysis not available in Octave, using pwelch.\n')
  %%       [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power in octave (pmtm is not implemented yet)
  %%     catch
  %%       error('pwelch function is needed for spectral analysis in Octave, please install the signal package from Octave Forge');
  %%     end
  %%   end

  %%   if i==1
  %%     % get size of spectrum and preallocate result matrix
  %%     nfreq=length(f);
  %%     Pxx=nan(nfreq,ncells);
  %%   end

  %%   if all(isnan(tmpPxx(:)))
  %%     tmpPxx=zeros(size(tmpPxx));
  %%   end

  %%   if ~isa(tmpPxx,'double')
  %%     % convert to double precision
  %%     tmpPxx=double(tmpPxx);
  %%   end

  %%   % smooth the spectrum
  %%   if smooth_factor>1 && strcmp(reportUI,'matlab')
  %%     tmpPxx=smooth(tmpPxx,smooth_factor);
  %%   else
  %%     tmpPxx=lsmooth(tmpPxx,smooth_factor);
  %%   end

  %%   % Peak Detection:
  %%   % select range of frequencies over which to look for peaks
  %%   sel = find(FreqRange(1)<=f & f<=FreqRange(end));

  %%   % set threshold for peak detection
  %%   ht=prctile(tmpPxx(sel),thresh_prctile); % ht=prctile(log10(tmpPxx(sel)),thresh_prctile);

  %%   if ~isnan(ht)
  %%     % get index of peaks in range over threshold
  %%     if strcmp(reportUI,'matlab')
  %%       [linPeakPower,PPind]=findpeaks(tmpPxx(sel),'MinPeakHeight',ht,'NPeaks',3); % [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
  %%     else
  %%       [linPeakPower,PPind]=findpeaks(tmpPxx(sel),'MinPeakHeight',ht,'MinPeakDistance',0,'MinPeakWidth',0);
  %%     end
  %%     PeakPower = log10(linPeakPower);
  %%   else
  %%     PPind=[];
  %%   end

  %%   if ~isempty(PPind)
  %%     % if multiple peaks, only consider the largest
  %%     if numel(PPind)>1
  %%       PPind=PPind(max(PeakPower)==PeakPower); %PPind=PPind(1);
  %%     end

  %%     % get frequency at that index
  %%     PeakFreq(i) = f(sel(PPind));

  %%     % set limits for calculating area under spectrum around peak
  %%     flo=PeakFreq(i)-Fwin/2;
  %%     fhi=PeakFreq(i)+Fwin/2;
  %%     sel2=(flo<=f & f<=fhi);
  %%     % calculate area under spectrum around peak
  %%     PeakArea(i) = sum(tmpPxx(sel2))*(f(2)-f(1));
  %%   else
  %%     PeakFreq(i)=nan;
  %%     PeakArea(i)=nan;
  %%   end
  %%   % Store results
  %%   Pxx(:,i)=tmpPxx;
  %% end
  % -----------------------------------------------------
  % Repeat spectral estimate for MUA:
  %% if ncells==1
  %%   % same as SUA
  %%   Pxx_mean=Pxx;
  %%   Pxx_mean_PeakFreq=PeakFreq;
  %%   Pxx_mean_PeakArea=PeakArea;
  %% else
    % calculate MUA
    % TODO AES try with and without detrend
    X=detrend(nanmean(dat(t1:t2,:),2)); % detrend the data

    %% [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power

    fprintf('About to start running coupling analysis')
    [pacmat, freqvec_ph, freqvec_amp, pmat, pac_angles] = ...
    find_pac_shf(X, Fs, options.measure, X, ...
    options.phase_freqs, options.ampl_freqs, options.plt,...
    options.waitbar, options.width, options.nfft, ...
    options.num_shf);

    if all(isnan(pacmat(:)))
      pacmat=zeros(size(pacmat));
    end
    if all(isnan(pac_angles(:)))
      pac_angles=zeros(size(pac_angles));
    end

    % TODO AES add back peak finding if necessary

    %% Add resulting power spectra to data structure
    % organization scheme:
    % data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
    % data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
    data.([var '_Coupling_MUA' options.output_suffix]).amplitudes=pacmat;
    data.([var '_Coupling_MUA' options.output_suffix]).angles=pac_angles;
    data.([var '_Coupling_MUA' options.output_suffix]).ampl_freq_axis=freqvec_amp;
    data.([var '_Coupling_MUA' options.output_suffix]).ph_freq_axis=freqvec_ph;

    if ~ismember([var '_Coupling_MUA' options.output_suffix],data.results)
      data.results{end+1}=[var '_Coupling_MUA' options.output_suffix];
    end

    if options.exclude_data_flag
      for l=1:length(data.labels)
        data=rmfield(data,data.labels{l});
      end
    end
end % end of options.variable loop

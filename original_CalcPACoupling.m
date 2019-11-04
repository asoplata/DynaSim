function data = original_CalcPACoupling(data,varargin)
%% data = CalcPACoupling(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   options:
%     'variable' - name of field containing data on which to calculate firing
%                rates (default: *_spikes or first variable in data.labels)
%     'time_limits' - [beg,end] (units of data.time)
%     'smooth_factor' - number of samples for smoothing the spectrum (default: 5)
%                       tip: set to 1 to avoid smoothing.
%   options for peak detection:
%     'min_peak_frequency' - Hz, min frequency for peak detection (default: 2)
%     'max_peak_frequency' - Hz, max frequency for peak detection (default: 150)
%     'peak_threshold_prctile' percentile for setting power threshold for peak detection (default: 95)
%     'peak_area_width' - Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum (default: 5)
%     'exclude_data_flag' - whether to remove simulated data from result structure (default: 0)
% 
% Outputs:
%   data: data structure with spectral power in data.VARIABLE_Power_SUA.Pxx
%         data.VARIABLE_Power_SUA.PeakFreq: frequency of spectral power (one value per cell)
%         data.VARIABLE_Power_SUA.PeakArea: area under spectrum around peak (one value per cell)
%   NOTE: for populations: spectrum of the mean waveform is stored in 
%         data.VARIABLE_Power_MUA.Pxx. population mean spectrum of the individual
%         waveforms can be calculated as mean(data.VARIABLE_Power_MUA.Pxx,2).
% 
% organization scheme for spectral results:
% data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
% data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
% 
% note:
% "variable" can be specified as the name of a variable listed in
% data.labels, a cell array of string listing variable names, or as a 
% regular expression pattern for identifying variables to process.
% See SelectVariables for more info on supported specifications.
% 
% Examples:
% s=[];
% s.populations(1).name='E';
% s.populations(1).equations='dv[2]/dt=@current+10; {iNa,iK}; v(0)=-65';
% s.populations(2).name='I';
% s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
% data=SimulateModel(s,'tspan',[0 1000]);
% data=CalcPACoupling(data,'variable','v');
% % Plot the spectrum of the E-cell average population voltage
% figure; plot(data.E_v_Power_MUA.frequency,data.E_v_Power_MUA.Pxx); 
% xlabel('frequency (Hz)'); ylabel('power'); xlim([0 200]);
% 
% See also: PlotPower, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0 Check inputs
options=dsCheckOptions(varargin,{...
  'variable',[],[],...        
  'time_limits',[-inf inf],[],...
  'smooth_factor',5,[],... % number of samples for smoothing the spectrum
  'min_peak_frequency',2,[],... % Hz, min frequency for peak detection
  'max_peak_frequency',200,[],... % Hz, max frequency for peak detection
  'peak_threshold_prctile',95,[],... % percentile for setting power threshold for peak detection
  'peak_area_width',5,[],... % Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum
  'exclude_data_flag',0,{0,1},...
  'plot_flag',0,[],...
  'slow_freq_range',[0 1.5],[],...
  'fast_freq_range',[8 13],[],...
  'window_length',2.0,[],...
  'window_overlap',1.0,[],...
  'number_bins',18,[],...
  'varlabel','V',[],...
  },false);

data = dsCheckData(data);
% note: calling CheckData() at beginning enables analysis function to
% accept data matrix [time x cells] in addition to DynaSim data structure.

if numel(data)>1
  % use AnalyzeStudy to recursively call CalcPACoupling on each data set
  data=dsAnalyzeStudy(data,@CalcPACoupling,varargin{:});
  return;
end

% time parameters
time = data.time; % time vector
dt = time(2)-time(1); % time step
ntime=length(time); % number of time points in full data set
t1=nearest(time,options.time_limits(1)); % index to first sample
t2=nearest(time,options.time_limits(2)); % index to last sample
nsamp=t2-t1+1; % number of samples for spectral estimate

% frequency parameters
Fs = fix(1/(dt/1000)); % effective sampling frequency
Fmin=options.min_peak_frequency; % min frequency for peak detection
Fmax=options.max_peak_frequency; % max frequency for peak detection
Fwin=options.peak_area_width; % size of frequency bin around peak for calculating area under spectrum
thresh_prctile=options.peak_threshold_prctile; % percentile for setting power threshold for peak detection
smooth_factor=options.smooth_factor; % number of samples to smooth spectrum
Fwin=options.peak_area_width; % size of frequency bin (centered on peak) over which to calculate area under spectrum
FreqRange=[max(Fmin,2/time(end)) Fmax]; % range to look for spectral peaks
NFFT=2^(nextpow2(nsamp-1)-1);%2); % <-- use higher resolution to capture STO freq variation
WINDOW=2^(nextpow2(NFFT-1)-3);
NOVERLAP=[]; % spectral parameters

%% 2.0 set list of variables to process as cell array of strings
options.variable=dsSelectVariables(data(1).labels,options.variable);

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
  % preallocation
  PeakFreq=nan(1,ncells);
  PeakArea=nan(1,ncells);   
  % SUA spectra: loop over cells
  for i=1:ncells
    % select data
    X=detrend(dat(t1:t2,i)); % detrend the data
    % calculate spectral estimate
    [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power
    if i==1
      % get size of spectrum and preallocate result matrix
      nfreq=length(f);
      Pxx=nan(nfreq,ncells);
    end
    if all(isnan(tmpPxx(:)))
      tmpPxx=zeros(size(tmpPxx));
    end
    if ~isa(tmpPxx,'double')
      % convert to double precision
      tmpPxx=double(tmpPxx);
    end
    if smooth_factor>1
      % smooth the spectrum
      tmpPxx=smooth(tmpPxx,smooth_factor);
    end
    % Peak Detection:
    % select range of frequencies over which to look for peaks
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    % set threshold for peak detection
    ht=prctile(log10(tmpPxx(sel)),thresh_prctile);
    if ~isnan(ht)
      % get index of peaks in range over threshold
      [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    else
      PPind=[];
    end
    if ~isempty(PPind)
      % if multiple peaks, only consider the largest
      if numel(PPind)>1
        PPind=PPind(max(PeakPower)==PeakPower); %PPind=PPind(1);
      end
      % get frequency at that index
      PeakFreq(i) = f(sel(PPind));
      % set limits for calculating area under spectrum around peak
      flo=PeakFreq(i)-Fwin/2;
      fhi=PeakFreq(i)+Fwin/2;
      sel2=(flo<=f & f<=fhi);
      % calculate area under spectrum around peak
      PeakArea(i) = sum(tmpPxx(sel2))*(f(2)-f(1));
    else
      PeakFreq(i)=nan;
      PeakArea(i)=nan;
    end
    % Store results
    Pxx(:,i)=tmpPxx;
  end
  % -----------------------------------------------------
  % Repeat spectral estimate for MUA: 
  if ncells==1
    % same as SUA
    Pxx_mean=Pxx;
    Pxx_mean_PeakFreq=PeakFreq;
    Pxx_mean_PeakArea=PeakArea;
  else
    % calculate MUA
    X=detrend(nanmean(dat(t1:t2,:),2)); % detrend the data
    % calculate spectral estimate
    [tmpPxx,f] = pwelch(X,NFFT,[],NFFT,Fs); % calculate power
    if all(isnan(tmpPxx(:)))
      tmpPxx=zeros(size(tmpPxx));
    end
    if ~isa(tmpPxx,'double')
      % convert to double precision
      tmpPxx=double(tmpPxx);
    end
    if smooth_factor>1
      % smooth the spectrum
      tmpPxx=smooth(tmpPxx,smooth_factor);
    end
    % Peak Detection:
    % select range of frequencies over which to look for peaks
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    % set threshold for peak detection
    ht=prctile(log10(tmpPxx(sel)),thresh_prctile);
    if ~isnan(ht)
      % get index of peaks in range over threshold
      [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    else
      PPind=[];
    end    
    if ~isempty(PPind)
      % if multiple peaks, only consider the largest
      if numel(PPind)>1
        PPind=PPind(max(PeakPower)==PeakPower); %PPind=PPind(1);
      end
      % get frequency at that index
      Pxx_mean_PeakFreq = f(sel(PPind));
      % set limits for calculating area under spectrum around peak
      flo=Pxx_mean_PeakFreq-Fwin/2;
      fhi=Pxx_mean_PeakFreq+Fwin/2;
      sel2=(flo<=f & f<=fhi);
      % calculate area under spectrum around peak
      Pxx_mean_PeakArea = sum(tmpPxx(sel2))*(f(2)-f(1));
    else
      Pxx_mean_PeakFreq=nan;
      Pxx_mean_PeakArea=nan;
    end    
    Pxx_mean=tmpPxx;
  end
  
  % Add resulting power spectra to data structure
  % organization scheme:
  % data.VARIABLE_Power_SUA.(Pxx,PeakFreq,PeakArea,frequency)
  % data.VARIABLE_Power_MUA.(Pxx,PeakFreq,PeakArea,frequency)
  data.([var '_Power_SUA']).Pxx=Pxx;
  data.([var '_Power_SUA']).PeakFreq=PeakFreq;
  data.([var '_Power_SUA']).PeakArea=PeakArea;
  data.([var '_Power_SUA']).frequency=f;
  data.([var '_Power_MUA']).Pxx=Pxx_mean;
  data.([var '_Power_MUA']).PeakFreq=Pxx_mean_PeakFreq;
  data.([var '_Power_MUA']).PeakArea=Pxx_mean_PeakArea;
  data.([var '_Power_MUA']).frequency=f;
  if ~ismember([var '_Power_SUA'],data.results)
    data.results{end+1}=[var '_Power_SUA'];
  end
  if ~ismember([var '_Power_MUA'],data.results)
    data.results{end+1}=[var '_Power_MUA'];
  end
  if options.exclude_data_flag
    for l=1:length(data.labels)
      data=rmfield(data,data.labels{l});
    end
  end
end


if ~isfield(spec,'entities') && isfield(spec,'cells')
  spec.entities=spec.cells;
elseif ~isfield(spec,'entities') && isfield(spec,'nodes')
  spec.entities=spec.nodes;
end

npop = length(spec.entities);
Fs = fix(data(1).sfreq);

for pop=1:npop
  labels = {data(pop).sensor_info.label};
  var = find(~cellfun(@isempty,regexp(labels,['_' options.varlabel '$'])),1,'first');
  if spec.entities(pop).multiplicity <= data(pop).epochs.num_trials
    n = spec.entities(pop).multiplicity;
  else
    n = data(pop).epochs.num_trials;
  end
  t = data(pop).epochs.time;
  dat = double(squeeze(data(pop).epochs.data(var,:,1:n))');
  if size(dat,2)>1
    lfp = mean(dat,1)';
  else
    lfp = dat;
  end
  try
    %% Filter the data
    % Note that `slow/fast_data` are thus TIMESERIES objects, not regular arrays.
    slow_data = idealfilter(timeseries(dat), [options.slow_freq_range(1)/(Fs/2.0), options.slow_freq_range(2)/(Fs/2.0)], 'pass');
    fast_data = idealfilter(timeseries(dat), [options.fast_freq_range(1)/(Fs/2.0), options.fast_freq_range(2)/(Fs/2.0)], 'pass');

    % HERE is probably where you'd want to use any of the available PAC
    % libraries, like the Eden/Kramer GLMCFC, since the signal has been
    % filtered, but not transformed into the analytic signal, and not windowed
    % yet (to get a phase-amplitude coupling "time-series" so to speak).

    %% Get angle and amplitude from the respective signals via Hilbert Transform
    % Kramer's GLMCFC code is espcially concise/informative. Note the '.Data'
    phi = angle(hilbert(slow_data.Data));
    amp = abs(hilbert(fast_data.Data));

    %% Construct windows across the time series
    %    Adapted & taken from Angela Onslow's 'pac_code_best/window_data.m' of
    %    her MATLAB Toolbox for Estimating Phase-Amplitude Coupling from
    %
    %    http://www.cs.bris.ac.uk/Research/MachineLearning/pac/
    %
    % Compute indices
    number_amp = size(amp,1);
    number_trials = size(amp,2);
    number_windows = ceil(options.window_length*Fs);
    number_overlaps = ceil(options.window_overlap*Fs);
    idx = bsxfun(@plus, (1:number_windows)', 1+(0:(fix((number_amp-number_overlaps)/(number_windows-number_overlaps))-1))*(number_windows-number_overlaps))-1;

    % Initialize the main data objects
    modulation_index_timeseries = [];
    modulogram_matrix = [];
    %% Loop over sliding windows
    for k=1:size(idx,2)
        amp_window = [];
        phi_window = [];
        % Loop over trials
        for j = 1:number_trials
            amp_window = [amp_window, amp(idx(:,k),j)];
            phi_window = [phi_window, phi(idx(:,k),j)];
        end

        %% Bin the faster frequency's amplitude in the slower's phase bins
        %    Adapted & taken from Adriano Tort's
        %    'Neurodynamics-master/16ch/Comodulation/ModIndex_v1.m' of the
        %    'Neurodynamics-Toolbox' repo on Github, at
        %
        %    https://github.com/cineguerrilha/Neurodynamics
        %
        phi_bin_beginnings = zeros(1,options.number_bins); % this variable will get the beginning (not the center) of each bin (in rads)
        bin_size = 2*pi/options.number_bins;
        for j=1:options.number_bins
            phi_bin_beginnings(j) = -pi+(j-1)*bin_size;
        end

        % Now we compute the mean amplitude in each phase:
        amp_means = zeros(1,options.number_bins);
        for j=1:options.number_bins
            phi_indices = find((phi_window >= phi_bin_beginnings(j)) & (phi_window < phi_bin_beginnings(j)+bin_size));
            amp_means(j) = mean(amp_window(phi_indices));
        end
        modulogram_matrix = [modulogram_matrix, (amp_means/sum(amp_means))'];

        % Quantify the amount of amp modulation by means of a normalized entropy index (Tort et al PNAS 2008):
        modulation_index=(log(options.number_bins)-(-sum((amp_means/sum(amp_means)).*log((amp_means/sum(amp_means))))))/log(options.number_bins);
        modulation_index_timeseries = [modulation_index_timeseries, modulation_index];

        % % Debug, for mid-function plotting:
        % % So note that the center of each bin (for plotting purposes) is phi_bin_beginnings+bin_size/2
        % % at this point you might want to plot the result to see if there's any amplitude modulation
        % figure(10)
        % bar(10:10:360,(amp_means/sum(amp_means)),'k')
        % xlim([0 360])
        % set(gca,'xtick',0:180:360)
        % xlabel('Phase (Deg)')
        % ylabel('Amplitude')
    end
  catch
    fprintf('\n!\n!\n!\n whoa somethings wrong, debug\n!\n!\n!')
  end
end

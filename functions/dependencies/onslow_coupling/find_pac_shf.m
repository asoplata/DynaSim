function [pacmat, freqvec_ph, freqvec_amp, pmat, pac_angles] = find_pac_shf (sig_pac, Fs, measure, ...
sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha,...
dataname, sig_pac_name, sig_mod_name)
% This function calculates a matrix of PAC values using either the ESC, MI 
% or CFC measure. 
% It uses shuffled datasets to conduct a significance analysis of the PAC
% values found
%
% AES notes:
%   - Function signature has been changed from its original! The
%     documentation has NOT been updated to reflect these changes.
%   - Significance detection has been commented out
%
% Basic function call:
% pacmat = find_pac_shf(sig_pac, Fs, measure)
%
% REQUIRED INPUTS:
%  sig_pac - signal suspected of containing PAC. Can be either a vector
%  (single trial) or a matrix (multiple trials). If input as a matrix,
%  trials are stored along the columns.
%  Fs - sampling frequency
%  measure - measure to be used - it should be: 'esc', 'mi' or 'cfc'
%
% The function can be executed with many optional inputs
%  function pacmat = find_pac_shf(sig_pac, Fs, measure, ...
%    sig_mod, ph_freq_vec, amp_freq_vec, plt, width, nfft, num_shf, alpha,...
%       dataname, sig_pac_name, sig_mod_name)
%
% OPTIONAL INPUTS:
%  sig_mod - signal containing modulating frequency; default = sig_pac
%  ph_freq_vec - range of frequencies for modulating signal; default = 1:5:101
%  amp_freq_vec - range of frequencies for PAC signal; default = 1:5:101
%  plt - flag indicating whether the output should be plotted - it should
%        be 'y' or 'n'; default = 'y'
%  waitbar - display progress in the command window; default = 1 (yes)
%  width - width of the wavelet filter; default = 7
%  nfft - the number of points in fft; default = 200
%  num_shf - the number of shuffled data sets to use during significance
%  testing; default = 50
%  alpha - significance value to use; default = 0.05 
%  dataname - the name of the dataset to be included as a graph title;
%             default = ''
%  sig_pac_name - the name of sig_pac to be printed on the y-axis; default = ''
%  sig_mod_name - the name of sig_mod to be printed on the x-axis; default = ''
% SHUFFLING ALGORITHM:
% If input is a single vector, shuffling will be done via random insertion
% type shuffling. Else, it will be done by randomly shuffling the trials.
%
% Author: Angela Onslow, May 2010
% 
% Modified by David Stanley, March 2015 to work with data in the form of
% multiple trials.

% Checks of input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    sig_mod = sig_pac;
end
if nargin < 5
    ph_freq_vec = 1:5:101;
end
if nargin < 6
    amp_freq_vec = 1:5:101;
end
if nargin < 7
    plt = 'y';
end
if nargin < 8
    waitbar = 1;
end
if nargin < 9
    width = 7;
end
if nargin < 10
    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
end
if nargin < 11
    num_shf = 50;
end
if nargin < 12
    alpha = 0.05;
end
if nargin < 13
    dataname = '';
end
if nargin < 14
    sig_pac_name = '';
end
if nargin < 15
    sig_mod_name = '';
end

% If data is a vector, make sure it's columnwise
if isvector(sig_pac)
    sig_pac = sig_pac(:);
end

if isvector(sig_mod)
    sig_mod = sig_mod(:);
end

if (size(sig_pac,2) ~= size(sig_mod,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some parameters for clarity
xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));
alpha = alpha/(xbins*ybins); % Uncomment to use Bonferonni Correction

% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each cell array element of ph_filt_signals and amp_filt_signals has the
% same dimensions as the original signals i.e. number of columns = number of
% trials
if (strcmp(measure, 'esc')) ||(strcmp(measure, 'mi')) 
[filt_sig_mod, filt_sig_pac] = filt_signalsWAV(sig_pac, sig_mod, Fs, ...
    ph_freq_vec, amp_freq_vec, measure, width);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Create shuffled datasets and distribution of PAC values %%%%%%%%%%%%%%%%%
% if num_shf ~= 0
% if waitbar == 1
%     fprintf('\nCreating shuffled data sets\n');
% end
% for s = 1:num_shf
%     
%     if strcmp(measure, 'esc')
%            
%         if isvector(sig_pac)
%             shuffled_sig_amp = shuffle_esc(filt_sig_pac, Fs);
%         else
%             shuffled_sig_amp = shuffle_trials_no_id(filt_sig_pac);
%         end
%         shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig_amp, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec,'n');
%         
%     end
%      
% 
%     if strcmp(measure, 'mi')
%         
%         if isvector(sig_pac)
%             shuffled_sig_amp = shuffle_esc(filt_sig_pac, Fs);
%         else
%             shuffled_sig_amp = shuffle_trials_no_id(filt_sig_pac);
%         end
%         shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig_amp, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec,'n');
%         
%     end
%     
%     if strcmp(measure, 'cfc')
%          
%         if isvector(sig_pac)
%             shuffled_sig1 = shuffle_esc(sig_pac, Fs);
%         else
%             shuffled_sig1 = shuffle_trials_no_id(sig_pac);
%         end
%         shf_pacmat_final(s,:,:) = find_pac_nofilt(shuffled_sig1, Fs,measure, sig_mod, ph_freq_vec, amp_freq_vec,'n', 0, width, nfft);
%         
%     end
%     
%     % Display current computational step to user
%     if waitbar == 1
%         if s == 1
%             fprintf('%03i%% ', floor((s/num_shf)*100));
%         else
%             fprintf('\b\b\b\b\b%03i%% ', floor((s/num_shf)*100));
%         end
%         if s == num_shf
%             fprintf('\n');
%         end
%     end
% end
% 
% %Find mean and standard deviation of shuffled data sets
% if strcmp(measure, 'mi')
%     for i =1:ybins
%         for j=1:xbins
%             [shf_data_mean(i,j), shf_data_std(i,j)] = normfit(shf_pacmat_final(:,i,j));
%         end
%     end
%     
% else
%     shf_data_mean = squeeze (mean (shf_pacmat_final, 1));
%     shf_data_std = squeeze (std (shf_pacmat_final, 1));
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate PAC measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(measure, 'esc')
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt(filt_sig_pac, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec, 'n', waitbar);
end

if strcmp(measure, 'mi')
    [pacmat, freqvec_ph, freqvec_amp, pac_angles] = find_pac_nofilt(filt_sig_pac, Fs, measure, filt_sig_mod, ph_freq_vec, amp_freq_vec, 'n', waitbar);
end

if strcmp(measure, 'cfc')
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_nofilt(sig_pac, Fs, measure, sig_mod, ph_freq_vec, amp_freq_vec, 'n', waitbar, width, nfft);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Compute significance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pacmat_plt = pacmat;

pmat = zeros(size(pacmat));

% if num_shf ~= 0 
% for i = 1:size(pacmat,1)
%     for j = 1:size(pacmat,2)
%         [h, p] = my_sig_test(pacmat(i,j), squeeze(shf_pacmat_final(:,i,j)), alpha);
%         pmat(i,j) = p;
%         if h == 0
%             pacmat_plt(i,j) = 0;
%         else
%             j;      % Dummy command to drop breakpoint
%         end
%     end
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% % Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if strcmp(plt, 'y')
%%     pac_plot_fun(pacmat, freqvec_ph, freqvec_amp, measure, sig_pac_name,...
%%                  sig_mod_name, dataname, pac_angles);
%% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

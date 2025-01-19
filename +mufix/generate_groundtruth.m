function [gndtruth, tvec] = generate_groundtruth(nsig, length__sec, samplingrate__hz, variationrate__hz, sig_range)
arguments
  nsig        (1,1) {mustBePositive, mustBeInteger}
  length__sec (1,1) {mustBePositive}
  samplingrate__hz  (1,1) {mustBePositive}
  variationrate__hz (1,1) {mustBePositive, mustBeLessThanOrEqual(variationrate__hz,samplingrate__hz)}
  sig_range   (:,2)
end 
%% 
%   [gndtruth, tvec] = generate_groundtruth(nsig, length__sec, fs__hz, nodes__hz, sig_range)
%
% Generate ground truth signals for photometry crosstalk simulations.
% To simulate the slower dynamics of the photometry response, grounnd truth
% signals are first generated in with a sampling rate of variationrate__hz, then
% up-sampled to samplingrate__hz.
% 
% Input:
%   nsig        - number of ground truth signals to generate
%   length__sec - length of the ground truth in seconds
%   samplingrate__hz  - sampling rate of the generated signals in Hz
%   variationrate__hz - number of signal values per second to simulate the
%                 slower dynamics of the photometry response
%   sig_range   - signal value ranges [min max]. Use a vertical vector for
%                 specifying different ranges for each ground truth signals.
%                 Otherwise, all signals share the same range.
%                 Default: [0 1]
% Output:
%   gndtruth    - nsig x time vector of generated ground truth
%   tvec        - corresponding timestamp vector

    %% Input processing
    if ~exist('sig_range','var')
      sig_range = [0 1];
    end
  
    if size(sig_range,1) > nsig
      sig_range = sig_range(1:nsig,:);
    elseif size(sig_range,1) < nsig
      if size(sig_range,1) == 1
        sig_range = repmat(sig_range, nsig, 1);
      else
        error('Number of signal ranges does not match number of ground truth signals specified.')
      end
    end

    %% Generate signals
    nNodes   = ceil(length__sec * variationrate__hz) + 1;
    gndtruth = rand(nsig,nNodes);
    
    % Scale to range
    sig_range(:,3) = sig_range(:,2) - sig_range(:,1);
    gndtruth = gndtruth .* sig_range(:,3) + sig_range(:,1);

    % Upsample
    Tx = (0:size(gndtruth,2)-1) / variationrate__hz;
    [gndtruth,tvec] = resample(gndtruth', Tx, samplingrate__hz);
    gndtruth = gndtruth';

    % Time vector
    tmask = tvec <= length__sec;
    tvec  = tvec(tmask);
    gndtruth = gndtruth(:,tmask);
   
end
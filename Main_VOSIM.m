function [results, SDANN, SDNNI, mse] = Main_VOSIM(InputSig,t,InputFormat,HRVparams,subjectID,annotations)
%  ====================== VOSIM Toolbox Main Script ======================
%
%   Main_VOSIM(InputSig,t,annotations,InputFormat,ProjectName,subjectID)
%	OVERVIEW:
%       Main "Validated Open-Source Integrated Matlab" VOSIM Toolbox script
%       Configured to accept RR intervals as well as raw data as input file
%
%   INPUT:
%       InputSig    - Vector containing RR interval data or ECG waveform  
%       t           - Time indices of the rr interval data (seconds) or
%                     ECG time
%       InputFormat - String that specifiy if the input vector is: 
%                     'RRinetrvals' for RR interval data 
%                     'ECGWaveform' for ECG waveform 
%       HRVparams   - struct of settings for hrv_toolbox analysis that can
%                     be obtained using InitializeHRVparams.m function 
%                     HRVparams = InitializeHRVparams();
%       subjectID   - (optional) string to identify current subject
%       annotations - (optional) annotations of the RR data at each point
%                     indicating the quality of the beat 
%
%   OUTPUT
%       results - HRV time and frequency domain metrics as well as AC and
%                 DC
%       SDANN   - standard deviation of the averages of values
%       SDNNI   - mean of the standard deviations of all values
%       mse     - Multiscale Entropy computed using SampleEntropy
%
%       NOTE: before running this script review and modifiy the parameters
%             in "initialize_HRVparams.m" file accordingly with the specific
%             of the new project (see the readme.txt file for further details)   
%   EXAMPLE
%       - rr interval input
%       Main_VOSIM(RR,t,'RRintervals',HRVparams)
%       - MIT Arrhythmia ECG wavefrom input
%       Main_VOSIM(ECGsig,t,'ECGWavefrom',HRVparams,'101')
%    
%   OUTPUT:
%       HRV Metrics 
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin < 4
    error('Wrong number of input arguments')
end
if nargin < 5
    subjectID = '0000';
    annotations = [];
end
if nargin < 6
    annotations = [];
end



%try    
    if strcmp(InputFormat, 'ECGWaveform')
        % Convert ECG waveform in rr intervals
        [t, rr, sqi] = ConvertRawDataToRRIntervals(InputSig, HRVparams, subjectID);
        avgLeadSQI = mean(sqi);
        GenerateHRVresultsOutput(subjectID,[],avgLeadSQI,'SQI','SQI',HRVparams,[],[]);    
    else
        rr = InputSig; 
        sqi = [];
    end

    % Exlude undesiderable data from RR series (i.e., arrhytmia, low SQI, ectopy, artefact, noise)

    [NN, tNN, fbeats] = RRIntervalPreprocess(rr,t,annotations, sqi, HRVparams);
    RRwindowStartIndices = CreateWindowRRintervals(tNN, NN, HRVparams);
    
    % 1. Atrial Fibrillation Detection
    if HRVparams.af.on == 1
        [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectID,tNN,NN,HRVparams);
        % Create RRAnalysisWindows contating AF segments
        RRwindowStartIndices = RomoveAFsegments(RRwindowStartIndices,AfAnalysisWindows, AFtest,HRVparams);
        fprintf('AF analysis completed for patient %s \n', subjectID);
    end
    
    % 2. Calculate time domain HRV metrics - Using VOSIM Toolbox Functions        

    [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
        RMSSD,pnn50,btsdet,avgsqi,fbeatw] = ...
        EvalTimeDomainHRVstats(NN,tNN,[],HRVparams,RRwindowStartIndices,fbeats);

    % 3. Frequency domain  metrics (LF HF TotPow) - Using VOSIM Toolbox Functions

     [ulf, vlf, lf, hf, lfhf, ttlpwr, methods, fdflag] = ...
         EvalFrequencyDomainHRVstats(NN,tNN, [],HRVparams,RRwindowStartIndices);
     
    % 4. PRSA
    try
        [ac,dc,~] = prsa(NN, tNN, [], RRwindowStartIndices, HRVparams);
    catch
        ac = NaN; 
        dc = NaN;
    end

    % 5.Export HRV Metrics as CSV File
    results = [RRwindowStartIndices(:), ac(:),dc(:),ulf(:),vlf(:),lf(:),hf(:), ...
               lfhf(:),ttlpwr(:),fdflag(:), NNmean(:),NNmedian(:), ...
               NNmode(:),NNvariance(:),NNskew(:),NNkurt(:),SDNN(:),...
               NNiqr(:),RMSSD(:),pnn50(:),btsdet(:),fbeatw(:)];

    col_titles = {'t_win','ac','dc','ulf','vlf','lf','hf',...
                  'lfhf','ttlpwr','fdflag','NNmean','NNmedian',...
                  'NNmode','NNvar','NNskew','NNkurt','SDNN',...
                  'NNiqr','RMSSD','pnn50','beatsdetected','corrected_beats'};

    % Save results
    ResultsFileName = GenerateHRVresultsOutput(subjectID,RRwindowStartIndices,results,col_titles, [],HRVparams, tNN, NN);
    
    fprintf('HRV metrics for file ID %s saved in the output folder in %s \n', subjectID, ResultsFileName);
    
    % 5. SDANN and SDNNi
    [SDANN, SDNNI] = ClalcSDANN(RRwindowStartIndices, tNN, NN(:),HRVparams); 

    
    
    % 6. Multiscale Entropy
    try
        mse = ComputeMultiscaleEntropy(NN,HRVparams.MSE.MSEpatternLength, HRVparams.MSE.RadiusOfSimilarity, HRVparams.MSE.maxCoarseGrainings);  
        % Save Results for MSE
        results = mse;
        col_titles = {'MSE'};
        % Generates Output - Never comment out
        GenerateHRVresultsOutput(subjectID,[],results,col_titles, 'MSE', HRVparams, tNN, NN);
    catch
        fprintf('MSE failed for file ID %s \n', subjectID);
    end

    
    fprintf('HRV Analysis completed for file ID %s \n',subjectID )
%catch
    
    results = NaN;
    col_titles = {'NaN'};
    GenerateHRVresultsOutput(subjectID,RRwindowStartIndices,results,col_titles, [],HRVparams, tNN, NN);    
    fprintf('Analysis not performed for file ID %s \n', subjectID);
%end

end %== function ================================================================
%




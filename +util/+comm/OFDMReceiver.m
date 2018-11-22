classdef OFDMReceiver < matlab.System
    % Untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    % Public, tunable properties
    properties
        
    end
    
    properties(DiscreteState)
        frequencyError
        symbolInterval
        lockIndicator
    end
    
    % Pre-computed constants
    properties(Access = private)
        sampleNumberFirstSample
        oscillatorPhaseFirstSample
        symbolPosition
        sampleBuffer
        symbolDurationExtended
        samplePointIndex
        trackingWindowPreLength
    end
    
    properties(Nontunable)
        guardIntervalDuration = 2048;
        symbolDuration = 8192;
        trackingWindow = 3;
        samplePoint = 0.8;
        filterLength = 1600;
    end
    
    methods
        function obj = OFDMReceiver(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.lockIndicator = false;
            obj.symbolDurationExtended = obj.guardIntervalDuration + obj.symbolDuration;
            obj.sampleBuffer = [];
            obj.sampleNumberFirstSample = 1;
            obj.samplePointIndex = round(obj.samplePoint*obj.guardIntervalDuration);
            obj.trackingWindowPreLength = ceil((obj.trackingWindow-1)/2);
            obj.symbolInterval = obj.symbolDurationExtended;
            obj.symbolPosition = [];
        end
        
        function [y, symbolPosition, symbolTimeshift,symbolIntervalCurrent2 ] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            trackingWindowWidth = obj.trackingWindowPreLength;
            symbolDurExt = obj.symbolDurationExtended;
            symbolDur = obj.symbolDuration;
            samplePointInd = obj.samplePointIndex;
            sampleGIRem = obj.guardIntervalDuration-samplePointInd;
            symbolPosition = [];
            symbolPositionOld = obj.symbolPosition;
            frequencyErr = obj.frequencyError;
            frequencyErrOld = frequencyErr;
            oscillatorPhaseStart = obj.oscillatorPhaseFirstSample;
            y = [];
            buffer = cat(1,obj.sampleBuffer,u);
            symbolIntervalOld = obj.symbolInterval;
            filterCoefficient = 1/obj.filterLength;
            symbolTimeshift = -sampleGIRem;
            symbolIntervalCurrent2=0;
            if obj.lockIndicator
                vmax = 0;
                for ind = 1:obj.trackingWindow
                    gi = buffer(ind:ind+2239);
                    sym = buffer(ind+8960:ind+11199);
                    gipower = sum(abs(gi).^2);
                    sympower = sum(abs(sym).^2);
                    vval = sum(conj(gi).*sym)/(gipower+sympower)*2;
                    if abs(vmax) < abs(vval)
                        vmax = vval;
                        indmax = ind;
                    end
                end
                if abs(vmax)>0.6
                        indmax = indmax + samplePointInd;
                        symbolPosition = obj.sampleNumberFirstSample + indmax-1;
                        symbolIntervalCurrent2 = symbolPosition-symbolPositionOld;
                        symbolIntervalCurrent = symbolIntervalOld*(1-filterCoefficient)+symbolIntervalCurrent2*filterCoefficient;
                        obj.symbolPosition = symbolPositionOld + symbolDurExt;
                        obj.symbolInterval = symbolIntervalCurrent;
                        symbolTimeshift = symbolPosition - symbolPositionOld - symbolIntervalCurrent + symbolTimeshift;
                        buffer = buffer(indmax:end);
                        %frequencyErr = angle(vmax)/symbolDur;
                        %obj.lockIndicator = true;
                else
                    error('test')
                end
            else
                if length(buffer)>=symbolDurExt*2
                    [vmax,indmax] = obj.syncOfdmFrame(buffer);
                    if abs(vmax)>0.6
                        indmax = indmax + samplePointInd;
                        symbolPosition = obj.sampleNumberFirstSample + indmax-1;
                        obj.symbolPosition = symbolPosition;
                        buffer = buffer(indmax:end);
                        frequencyErr = angle(vmax)/symbolDur;
                        obj.lockIndicator = true;
                    else
                        error('test')
                    end
                else
                    obj.sampleBuffer = buffer;
                    return
                end
            end
            oscillatorPhase = oscillatorPhaseStart+(symbolPosition-obj.sampleNumberFirstSample)*2*pi*frequencyErrOld;
            modulator = exp(-1i*frequencyErr*(symbolPosition:symbolPosition+symbolDur-1).');
            
            freqa=[0:ceil(symbolDur/2)-1,-floor(symbolDur/2):-1]';
            obj.frequencyError = frequencyErr;
            y = fft(modulator.*buffer(1:symbolDur));
            y = y.*exp(-2i*pi/symbolDur*symbolTimeshift*freqa);
            obj.sampleBuffer = buffer(symbolDur+1 + sampleGIRem - trackingWindowWidth:end);
            sampleNumberFirstSampleNew = symbolPosition + symbolDur + sampleGIRem - trackingWindowWidth;
            oscillatorPhaseStart = oscillatorPhase + (sampleNumberFirstSampleNew-symbolPosition)*2*pi*frequencyErr;
            obj.sampleNumberFirstSample =sampleNumberFirstSampleNew;
            obj.oscillatorPhaseFirstSample = oscillatorPhaseStart;
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        function [vmax,indmax] = syncOfdmFrame(obj,kbuff)
            mval = 0;
            giend = obj.guardIntervalDuration-1;
            symstart = obj.symbolDuration;
            symend = obj.symbolDurationExtended-1;
            for ind = 1:length(kbuff)-symend
                gi = kbuff(ind:ind+giend);
                sym = kbuff(ind+symstart:ind+symend);
                gipower = sum(abs(gi).^2);
                sympower = sum(abs(sym).^2);
                vval = sum(conj(gi).*sym)/(gipower+sympower)*2;
                if abs(mval) < abs(vval)
                    mval = vval;
                    framepos = ind;
                end
            end
            vmax = mval;
            indmax = framepos;
        end
    end
end

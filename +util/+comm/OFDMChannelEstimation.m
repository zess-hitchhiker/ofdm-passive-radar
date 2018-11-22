classdef OFDMChannelEstimation < matlab.System
    % Untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        
    end

    properties(DiscreteState)
        pilotState
    end

    % Pre-computed constants
    properties(Access = private)
        buffer
        nbrPilotState
        window
    end
    
    properties(Nontunable)
        pilotIndices
        pilotValues
        guardInterval = 0.25;
    end
    
    methods
        function obj = OFDMChannelEstimation(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj,u)
            % Perform one-time calculations, such as computing constants
            
            obj.buffer = zeros(size(u));
            
            obj.nbrPilotState = length(obj.pilotIndices);
            localBuffer = zeros(length(u),obj.nbrPilotState);
            
            for ind = 1:obj.nbrPilotState
                pindArr = obj.pilotIndices{ind}+1;
                pvalArr = obj.pilotValues{ind};
                localBuffer(pindArr,ind) = u(pindArr)./pvalArr;
            end
            [~,ind] = max(max(ifft(localBuffer)));
            obj.pilotState = ind;
            nwin = floor(length(u)*obj.guardInterval)/4;
            obj.window = circshift(cat(2,tukeywin(nwin,0.5).',zeros(1,length(u)-nwin)),-round(nwin/4));
        end

        function [y, pstate, localBuffer2] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            localBuffer = obj.buffer;
            localBuffer2 = zeros(size(obj.buffer));
            pstate = obj.pilotState;
            pindArr = obj.pilotIndices{pstate}+1;
            pvalArr = obj.pilotValues{pstate};
            
            localBuffer(pindArr) = u(pindArr)./pvalArr;
            localBuffer2(pindArr) = u(pindArr)./pvalArr;
            obj.buffer = localBuffer;
            windowfun = obj.window;
            localBuffer = fft(ifft(localBuffer).*windowfun)*3;
            %localBuffer = interp1(1:3:length(localBuffer),localBuffer(1:3:end),1:length(localBuffer));
            obj.pilotState= mod(pstate,obj.nbrPilotState)+1;
            
            y = localBuffer;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end

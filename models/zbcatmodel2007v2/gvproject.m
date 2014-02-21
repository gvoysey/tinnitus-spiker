function [result] =gvproject(freq)
    lowSpont = modelPopulation(10,1,freq);
    highSpont = modelPopulation(80,1,freq);
    lowSpont50 = modelPopulation(10,.5,freq);
    highSpont50 = modelPopulation(80,.5,freq);
    lowSpont25 = modelPopulation(10,.25,freq);
    highSpont25= modelPopulation(80,.25,freq);
    
    result = struct;
    result.lowSpont = lowSpont;
    result.highSpont = highSpont;
    result.lowSpont25 = lowSpont25;
    result.highSpont25 = highSpont25;
    result.lowSpont50 = lowSpont50;
    result.highSpont50 = highSpont50;
    %fibers at CF, lowspont undamaged vs highspont undamaged
    figure;hold on;
    subplot(2,2,1);plot(result.lowSpont.times(3,:),result.lowSpont.psths(3,:));xlabel('time(s)');title('high spont')
    subplot(2,2,3);spectrogram(result.lowSpont.psths(3,:),hamming(64),32,64,200e3);
    
    subplot(2,2,2);plot(result.highSpont.times(3,:),result.highSpont.psths(3,:));xlabel('time(s)');title('low spont')
    subplot(2,2,4);spectrogram(result.highSpont.psths(3,:),hamming(64),32,64,200e3);
    [ax,h3] = suplabel('Undamaged fibers; CF = F0 = 150 Hz','t');
    
    %fibers at CF, lowspont undamaged vs lowspont damaged
    figure;hold on;
    subplot(2,2,1);plot(result.lowSpont.times(3,:),result.lowSpont.psths(3,:));xlabel('time(s)');title('low spont undamaged')
    subplot(2,2,3);spectrogram(result.lowSpont.psths(3,:),hamming(64),32,64,200e3);
    
    subplot(2,2,2);plot(result.lowSpont25.times(3,:),result.lowSpont25.psths(3,:));xlabel('time(s)');title('low spont 75% damaged')
    subplot(2,2,4);spectrogram(result.lowSpont25.psths(3,:),hamming(64),32,64,200e3);
    [ax,h3] = suplabel('CF = F0 = 150 Hz','t');
    
    % mean spont rates for fibers across -2,-1,0,1,2 CF. 25% and
    % 100% functional
    figure; hold on;
    lowbars = zeros(2,5);
    for i = 1:5,
        lowbars(1,i)=mean(result.lowSpont.psths(i,:));
        lowbars(2,i)=mean(result.lowSpont25.psths(i,:));
    end        
    highbars = zeros(2,5);
    for i = 1:5,
        highbars(1,i)=mean(result.highSpont.psths(i,:));
        highbars(2,i)=mean(result.highSpont25.psths(i,:));
    end 
    subplot(1,2,1);bar(lowbars); xlabel('fibers around CF');title('low-spont');
        set(gca,'XTickLabel',{'undamaged','75% damaged'});
    subplot(1,2,2);bar(highbars); xlabel('fibers around CF');title('high-spont');
        set(gca,'XTickLabel',{'undamaged','75% damaged'});
    
end

function [psths] = modelPopulation(spont,ihc, freq)
%define a vector of 24 CFs as third-octave CFs:
    %cfs = (10.0^3) * ((2.0) .^ ([-18:13]./3));
    %cfs = cfs(9:end); %model restictions <80 Hz
    cfs = [99.2126  125.0000  157.4901  198.4251  250.0000];    
    cohcs = ones(1,length(cfs));
    cihcs = ihc.*ones(1,length(cfs));
    %sponts = 50*ones(1,24);
    %f0s = cfs;
    psths = struct;   
    for i = 1:length(cfs)
        cf = cfs(i);
        cohc = cohcs(i);
        cihc = cihcs(i);
        %spont = sponts(i);
        f0 = freq;%f0s(i);
        graphs = false;
        [psth,pstime] = modelFiber(cf, cohc,cihc,spont,f0,graphs);
        psths.psths(i,:)= psth;
        psths.times(i,:)= pstime;
        figure;
        hold on;
        subplot(2,1,1);plot(pstime,psth);xlabel('time, s');
        subplot(2,1,2);spectrogram(psth);
        title(['CF= ' num2str(cf) ' Hz. Spont rate' num2str(spont) 'spikes/sec, IHC percentage' num2str(ihc)]);        
    end  
end

function [psth,psthtime] = modelFiber(Cf, Cohc, Cihc, Spont, f0, graphs)
    % model fiber parameters
    CF    = Cf; % 1e3;  % CF in Hz;   
    cohc  = Cohc; %1.0;  % normal ohc function
    cihc  = Cihc; %1.0;  % normal ihc function
    spont = Spont; % 50; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects
    % stimulus parameters
    F0 = f0;     % stimulus frequency in Hz
    Fs = 200e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
    T  = 100e-3;  % stimulus duration in seconds
    rt = 5e-3;   % rise/fall time in seconds
    stimdb = 50; % stimulus intensity in dB SPL
    % PSTH parameters
    nrep = 50;            % number of stimulus repetitions (e.g., 50);
    psthbinwidth = 0.2e-3; % binwidth in seconds;

    t = 0:1/Fs:T-1/Fs; % time vector
    mxpts = length(t);
    irpts = rt*Fs;

    pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

    [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth500k] ...
        = zbcatmodel(pin,CF,nrep,1/Fs,T*1.5,cohc,cihc,spont);

    psthbins = round(psthbinwidth*Fs);  % number of psth500k bins per psth bin
    psthtime = timeout(1:psthbins:end); % time vector for psth
    pr = sum(reshape(psth500k,psthbins,length(psth500k)/psthbins))/nrep; % pr of spike in each bin
    psth = pr/psthbinwidth; % psth in units of spikes/s
    
    if(graphs)
        figure
        subplot(2,1,1)
        plot(timeout,[pin zeros(1,length(timeout)-length(pin))])
        title('pin')
        yl1 = ylim;
        subplot(2,1,2)
        plot(timeout,meout)
        title('meout')
        xlabel('Time (s)')
        yl2 = ylim;
        yl = [min(yl1(1),yl2(1)) max(yl1(2),yl2(2))];
        subplot(2,1,1)
        ylim(yl)
        subplot(2,1,2)
        ylim(yl)

        figure
        subplot(2,1,1)
        plot(timeout,c1filterout)
        title('c1filterout')
        yl1 = ylim;
        subplot(2,1,2)
        plot(timeout,c2filterout)
        title('c2filterout')
        xlabel('Time (s)')
        yl2 = ylim;
        yl = [min(yl1(1),yl2(1)) max(yl1(2),yl2(2))];
        subplot(2,1,1)
        ylim(yl)
        subplot(2,1,2)
        ylim(yl)

        figure
        subplot(3,1,1)
        plot(timeout,c1vihc)
        title('c1vihc')
        yl1 = ylim;
        subplot(3,1,2)
        plot(timeout,c2vihc)
        title('c2vihc')
        xlabel('Time (s)')
        yl2 = ylim;
        subplot(3,1,3)
        plot(timeout,vihc)
        title('vihc')
        yl3 = ylim;
        yl = [min([yl1(1) yl2(1) yl3(1)]) max([yl1(2) yl2(2) yl3(2)])];
        subplot(3,1,1)
        ylim(yl)
        subplot(3,1,2)
        ylim(yl)
        subplot(3,1,3)
        ylim(yl)

        figure
        subplot(2,1,1)
        plot(timeout,synout)
        title('synout')
        xlabel('Time (s)')
        subplot(2,1,2)
        plot(psthtime,psth)
        title('psth')
        xlabel('Time (s)')
    end
end

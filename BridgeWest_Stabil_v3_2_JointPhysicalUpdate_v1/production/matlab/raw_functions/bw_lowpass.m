function y=bw_lowpass(y,fs,fc),[b,a]=butter(4,fc/(fs/2),'low');y=filtfilt(b,a,y);end


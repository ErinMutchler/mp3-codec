fs = 48000;
t = 0:1/fs:10;
                   
insig = sin(2*pi*200*t) + sin(2*pi*562*t) + sin(2*pi*1956*t);


frame = zeros(1, 256);


for i = 1:256
    frame(i) = insig(i + 256);
end

winFrame = frame'.*hann(256);

nfft = 8192;
FR = fft(winFrame, nfft);
FR = FR(1:nfft/2);
mx = abs(FR);
p = (10 * log10(mx)) + 90.3;
f = (0:nfft/2-1)*fs/nfft;

semilogx(p);

% find tonal maskers
tonalFlags = zeros(length(p));
pTM = zeros(length(p));
for k = 1:(length(p) - 1)
    pdk = 0; mdk = 0;
    
    if (k > 2) && (k < 63)
        pdk = 2;
        mdk = 0;
    elseif (k > 62) && (k < 127)
        pdk = 3;
        mdk = 2;
    elseif (k > 126) && (k < 256)
        pdk = 6;
        mdk = 2;
    end
    if k + pdk > length(p)
        pdk = 0;
    end
    
    if (p(k) > p(k + 1)) && (p(k) > p(k - 1))
        tonalFlags(k) = 1;
    elseif (p(k) > p(k + pdk) + 7) && (p(k) > p(k - mdk) + 7)
        tonalFlags(k) = 1;
    end
    
    if tonalFlags(k) == 1
       sumPeaks = 0;
       for j = -1:1
           sumPeaks = sumPeaks + 10^(0.1 * p(k + j));
       end
       pTM(k) = 10 * log10(sumPeaks);
    end
end


% find noise maskers
pNM = zeros(length(p));
for k = 1:length(p)
    pdk = 0;
    mdk = 0;
    if (k > 2) && (k < 63)
        pdk = 2;
        mdk = 0;
    elseif (k > 62) && (k < 127)
        pdk = 3;
        mdk = 2;
    elseif (k > 126) && (k < 256)
        pdk = 6;
        mdk = 2;
    end
    
    for j = 1:length(p)
        if (p(j) == pTM(k)) || (p(j) == pTM(k + 1)) || (p(j) == pTM(k - 1)) || (p(j) == pTM(k + pdk)) || (p(j) == pTM(k - pdk))
           pNM(k) = 0;
           break;
        else
            pNM(k) = pNM(k) + 10^(0.1 * p(j));
        end
    end
    pNM(k) = 10 * log10(pNM(k));
end


% decimate maskers under threshold of hearing
for k = 1:length(p)
    Tq = 3.64*(k/1000)^-0.8 - (6.5 * exp(1))^(-0.6 * (k/1000 - 3.3)^2) + (10^-3) * (k/1000)^4;
    if pTM(k) <= Tq
        pTM(k) = 0;
    end
    if pNM(k) <= Tq
        pNM(k) = 0;
    end
end

% decimate maskers in the same bark band

barkOfLastTonalMasker = 0;
barkOfLastNoiseMasker = 0;
lastTonalMasker = 0;
lastNoiseMasker = 0;
for k = 1:length(p)
    bark = (13 * atan(0.00076 * k)) + (3.5 * atan((k/7500)^2));
    if pTM(k) ~= 0
        if barkOfLastTonalMasker == bark
            if pTM(k) >= pTM(lastTonalMasker)
                pTM(lastTonalMasker) = 0;
            elseif pTM(k) < pTM(lastTonalMasker)
                pTM(k) = 0;
            end
        end
        lastTonalMasker = k;
        barkOfLastTonalMasker = bark;
    end
end


% decimate 2:1 and 4:1 maskers for mid and then high freqs
midFlag = 0;
highFlag = 0;
for k = 1:232
    if (k >= 49) && (k <= 96)
        if (pTM(k) ~= 0) && (midFlag == 1)
            pTM(k) = 0;
            midFlag = 0;
        elseif (pTM(k) ~= 0) && (midFlag == 0)
            midFlag = 1;
        end
    end
    if (k >= 97) && (k <= 232)
        if (pTM(k) ~= 0) && (highFlag == 3)
            pTM(k) = 0;
            highFlag = 0;
        elseif (pTM(k) ~= 0) && (highFlag ~= 4)
            highFlag = highFlag + 1;
        end
    end
end

% calculate individual masking curves for tonal maskers
ttm = zeros(length(p), 256);
for i = 1:256
    for j = 1:256
        z = (13 * atan(0.00076 * j)) + (3.5 * atan((j/7500)^2));
        SF = 0;
        deltaZ = (13 * atan(0.00076 * j)) + (3.5 * atan((j/7500)^2)) - (13 * atan(0.00076 * i)) + (3.5 * atan((i/7500)^2));
        if (deltaZ >= -3 && deltaZ < -1)
            SF = (17 * deltaZ) - (0.4 * pTM(j)) + 11;
        elseif (deltaZ >= -1 && deltaZ < 0)
            SF = ((0.4 * pTM(j)) + 6) * deltaZ;
        elseif (deltaZ >= 0 && deltaZ < 1)
            SF = -17 * deltaZ;
        elseif (deltaZ >= 1 && deltaZ < 8)
            SF = (((0.15 * pTM(j)) - 17) * deltaZ) - (0.15 * pTM(j));
        end
        
        ttm(i, j) = pTM(j) - (0.275 * z) + SF - 6.025;
    end
end

% calculate individual masking curves for noise maskers
tnm = zeros(length(p), 256);
for i = 1:256
    for j = 1:256
        z = (13 * atan(0.00076 * j)) + (3.5 * atan((j/7500)^2));
        SF = 0;
        deltaZ = (13 * atan(0.00076 * j)) + (3.5 * atan((j/7500)^2)) - (13 * atan(0.00076 * i)) + (3.5 * atan((i/7500)^2));
        if (deltaZ >= -3 && deltaZ < -1)
            SF = (17 * deltaZ) - (0.4 * pNM(j)) + 11;
        elseif (deltaZ >= -1 && deltaZ < 0)
            SF = ((0.4 * pNM(j)) + 6) * deltaZ;
        elseif (deltaZ >= 0 && deltaZ < 1)
            SF = -17 * deltaZ;
        elseif (deltaZ >= 1 && deltaZ < 8)
            SF = (((0.15 * pTM(j)) - 17) * deltaZ) - (0.15 * pNM(j));
        end
        
        ttm(i, j) = pNM(j) - (0.275 * z) + SF - 6.025;
    end
end

sumOfTTM = 0;
sumOfTNM = 0;
Tg = zeros(length(p));

for i = 1:256
    Tq = 3.64*(i/1000)^-0.8 - (6.5 * exp(1))^(-0.6 * (i/1000 - 3.3)^2) + (10^-3) * (i/1000)^4;
    for l = 1:256
        sumOfTTM = sumOfTTM + 10^(0.1 * ttm(i, l));
        sumOfTNM = sumOfTNM + 10^(0.1 * tnm(i, l));
    end
    Tg(i) = 10 * log10(10^(0.1 * Tq) + sumOfTTM + sumOfTNM);
end

%semilogx(Tg);

# DSP Lab Report — LAB(1) Assignment

**University:** Zagazig University — Faculty of Engineering  
**Student:** Abdelrahman Husseini Abdallah Khalil  
**Year:** 3rd Year — Lab Section  
**Topic & Tool:** IIR Filter Design | GNU Octave

---

## Table of Contents

| Section | Topic |
|---------|-------|
| [Q1](#q1--butterworth-bandpass-filter) | Butterworth Bandpass Filter (n=4 and n=21) |
| [S1](#s1--stability-analysis) | Stability Analysis — Pole-Zero Check |
| [Q2](#q2--whistlewav-noise-rejection) | whistle.wav — Band-Stop Noise Rejection |
| [Q3](#q3--fir-vs-iir-moving-average-filters) | FIR vs IIR Moving Average Filters |
| [C3](#c3--comments-on-the-plots) | Comments on the Plots |

---

## Q1 — Butterworth Bandpass Filter

**Fs = 10 kHz | fc1 = 200 Hz | fc2 = 2000 Hz**

Design a Butterworth Bandpass filter for n=4 and n=21. Plot Frequency Response and Impulse Response with labelled axes. Check stability.

> **Key Rule:** Normalise frequency to Nyquist: `W = f / (Fs/2)`

---

### Order n = 4 — Complete Code

```octave
clear; clc; close all
Fs=10000; fc1=200; fc2=2000; n=4;
[b,a] = butter(n, [fc1 fc2]/(Fs/2), 'bandpass');

%% Impulse Response
figure; impz(b,a);
title('Impulse Response - Butterworth Bandpass (n=4)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

%% Frequency Response
f = (0:0.001:1)*(Fs/2);
H = freqz(b,a,f,Fs);
figure; plot(f,abs(H));
title('Frequency Response - Butterworth Bandpass (n=4)');
xlabel('Frequency (Hz)'); ylabel('|H(f)|'); grid on;
```

### Plots — n = 4

| Impulse Response | Frequency Response |
|:---:|:---:|
| ![Impulse Response n=4](Q1_1at_n4.jpg) | ![Frequency Response n=4](Q1_2at_n4.jpg) |
| *Figure 1a: Impulse Response (n=4)* | *Figure 1b: Frequency Response (n=4)* |

---

### Order n = 21 — Complete Code

```octave
clear; clc; close all
Fs=10000; fc1=200; fc2=2000; n=21;
[b,a] = butter(n, [fc1 fc2]/(Fs/2), 'bandpass');

%% Impulse Response
figure; impz(b,a);
title('Impulse Response - Butterworth Bandpass (n=21)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

%% Frequency Response
f = (0:0.001:1)*(Fs/2);
H = freqz(b,a,f,Fs);
figure; plot(f,abs(H));
title('Frequency Response - Butterworth Bandpass (n=21)');
xlabel('Frequency (Hz)'); ylabel('|H(f)|'); grid on;
```

### Plots — n = 21

| Impulse Response | Frequency Response |
|:---:|:---:|
| ![Impulse Response n=21](Q1_3at_n21.jpg) | ![Frequency Response n=21](Q1_4at_n21.jpg) |
| *Figure 2a: Impulse Response (n=21)* | *Figure 2b: Frequency Response (n=21)* |

---

### Observations

| Property | n = 4 | n = 21 |
|----------|-------|--------|
| Transition band | Wide — gradual rolloff | Very narrow — sharp rolloff |
| Stopband rejection | Moderate | Very strong |
| Impulse response | Short | Long |
| Computation cost | Low | High |

---

## S1 — Stability Analysis

Butterworth is an IIR filter. **Stability condition:** all poles inside the unit circle.  
MATLAB/Octave `butter()` always places poles inside the unit circle by design.  
Verify practically with `zplane(b,a)` — all **×** marks must be inside the circle.

```octave
zplane(b,a) % all x (poles) inside unit circle => STABLE
```

| Order | Stable? | Reason |
|-------|---------|--------|
| n = 4 | ✅ YES — all poles inside unit circle | Butterworth design guarantees poles inside unit circle |
| n = 21 | ✅ YES — all poles inside unit circle | Butterworth design guarantees poles inside unit circle |

---

## Q2 — whistle.wav Noise Rejection

Read `whistle.wav`, identify noise frequencies via FFT, design two Butterworth Band-Stop filters to reject them, verify stability, compare energy before/after, play the result.

### Steps 1 & 2 — Read File, Info & Plot Spectrum

```octave
pkg load signal
[x, fs] = audioread('D:\whistle.wav');
n = length(x); time = n/fs; t = linspace(0,time,n);
fprintf('Samples: %d | Duration: %.4f s | Fs: %d Hz', n, time, fs);
% Results: Fs=24000 Hz, N=41472, Duration=1.728 s

sound(x, fs, 16); pause(time+1);
X = fft(x); f = (0:n-1)*(fs/n);
figure; plot(f(1:n/2), abs(X(1:n/2)));
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
title('Spectrum of Original Signal x'); grid on;
% --> Two sharp peaks visible at 500 Hz and 1500 Hz
```

| Parameter | Value |
|-----------|-------|
| Fs | 24,000 Hz |
| N | 41,472 samples |
| Duration | 1.7280 s |
| Noise peaks | 500 Hz and 1500 Hz |

### Steps 3 & 4 — Design Filters, Apply, Energy & Playback

```octave
% -- Filter 1: reject 500 Hz
fc1=480; fc2=520;
[b1,a1] = butter(4, [fc1 fc2]/(fs/2), 'stop');
[H1,w1] = freqz(b1,a1,1024,fs);
figure; plot(w1,abs(H1)); title('Freq Response - Band-Stop @ 500 Hz');
figure; impz(b1,a1,100); title('Impulse Response - Band-Stop @ 500 Hz');
poles1=roots(a1); fprintf('Filter1 stable: %d', all(abs(poles1)<1));

% -- Filter 2: reject 1500 Hz
fc3=1450; fc4=1550;
[b2,a2] = butter(4, [fc3 fc4]/(fs/2), 'stop');
[H2,w2] = freqz(b2,a2,1024,fs);
figure; plot(w2,abs(H2)); title('Freq Response - Band-Stop @ 1500 Hz');
figure; impz(b2,a2,100); title('Impulse Response - Band-Stop @ 1500 Hz');
poles2=roots(a2); fprintf('Filter2 stable: %d', all(abs(poles2)<1));

% -- Apply cascade & compute energy
y = filter(b1,a1,x); y = filter(b2,a2,y);
Ex = sum(x.^2); Ey = sum(y.^2);
fprintf('Energy x: %.4f | Energy y: %.4f', Ex, Ey);
fprintf('Removed: %.4f (%.2f%%)', Ex-Ey, (Ex-Ey)/Ex*100);

% -- Spectrum of filtered signal
Y = fft(y);
figure; plot(f(1:n/2), abs(Y(1:n/2))); title('Spectrum of Filtered Signal y');

% -- Play filtered signal
y_play = y/max(abs(y)); sound(y_play, fs, 16); pause(time+1);
```

### Plots — Q2 Results

![Q2 All Figures](Q2_1.png)
*Figure 3: Original Spectrum, Band-Stop Filter Responses (500 Hz & 1500 Hz), Impulse Responses, and Filtered Spectrum*

### Results Summary

| Parameter | Value |
|-----------|-------|
| Filter 1 stopband | 480–520 Hz (rejects 500 Hz) — ✅ STABLE |
| Filter 2 stopband | 1450–1550 Hz (rejects 1500 Hz) — ✅ STABLE |
| Pole magnitudes | All < 1 (inside unit circle) |
| Energy original x | 585.8219 |
| Energy filtered y | 324.8403 |
| Energy removed | 260.9817 (44.55%) |
| Whistle after filter | Gone — not audible ✅ |

---

## Q3 — FIR vs IIR Moving Average Filters

For each filter plot the Frequency Response and Impulse Response with labelled axes and titles.

| Part | Difference Equation |
|------|---------------------|
| (a) | `y[n] = 1/8 ( x[n] + x[n-1] + x[n-2] + x[n-3] + x[n-4] + x[n-5] + x[n-6] + x[n-7] )` |
| (b) | `y[n] = 1/8 x[n] - 1/8 x[n-8] + y[n-1]` |

---

### Part (a) — FIR Moving Average LPF

Non-recursive filter. Equal weight 1/8 for 8 samples. All `a_k = 0` except `a_0`.

```octave
a = 1; Fs = 8000;
b = 0.125 * [1 1 1 1 1 1 1 1]; % FIR coefficients

f = linspace(0, Fs/2, 1024);
H = freqz(b, a, f, Fs);
figure; plot(f, abs(H));
title('Frequency Response - FIR Moving Average LPF (n=7)');
xlabel('Frequency (Hz)'); ylabel('Magnitude |H(f)|'); grid on;

figure; impz(b, a);
title('Impulse Response - FIR Moving Average LPF (n=7)');
xlabel('Samples (n)'); ylabel('Amplitude'); grid on;
```

### Plots — Part (a) FIR

![FIR Moving Average](Q3_1.png)
*Figure 4: FIR Moving Average — Frequency Response (left) and Impulse Response (right)*

| Property | Value |
|----------|-------|
| Type | FIR — Non-recursive |
| b vector | `[0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]` |
| Stable | ✅ Always stable — FIR filters are unconditionally stable |
| Impulse | Exactly 8 non-zero samples then zero **(FINITE)** |

---

### Part (b) — IIR Recursive Moving Average

Recursive equivalent of (a). Same averaging result with only 2 multiplications instead of 8.

```octave
b = [1/8, 0, 0, 0, 0, 0, 0, 0, -1/8];
a = [1, -1];
Fs = 8000;

[H, f] = freqz(b, a, 1024, Fs);
figure; plot(f, abs(H));
title('Frequency Response - IIR Recursive Moving Average');
xlabel('Frequency (Hz)'); ylabel('Magnitude |H(f)|'); grid on;

figure; impz(b, a, 50);
title('Impulse Response - IIR Recursive Moving Average');
xlabel('Samples (n)'); ylabel('Amplitude'); grid on;
```

### Plots — Part (b) IIR

![IIR Recursive Moving Average](Q3_2.png)
*Figure 5: IIR Recursive Moving Average — Frequency Response (left) and Impulse Response (right)*

| Property | Value |
|----------|-------|
| Type | IIR — Recursive |
| b vector | `[1/8, 0, 0, 0, 0, 0, 0, 0, -1/8]` |
| a vector | `[1, -1]` |
| Efficient | Only 2 multiplications vs 8 for FIR |
| Impulse | Decays gradually — theoretically **INFINITE** |

---

## C3 — Comments on the Plots

| Property | Part (a) FIR | Part (b) IIR |
|----------|-------------|-------------|
| Filter type | FIR — non-recursive | IIR — recursive |
| Impulse response | Finite: 8 non-zero samples | Infinite: decays, never truly zero |
| Freq response | LPF with ripple in stopband | Nearly identical LPF shape |
| Computation | 8 multiplications per sample | 2 multiplications per sample |
| Stability | Always stable (unconditional) | Marginally stable (pole at z=1) |
| Phase | Linear phase | Non-linear phase |
| Conclusion | Simple but costly at high order | Efficient recursive shortcut |

> **Key conclusion:** Both parts produce the same LPF moving average and nearly identical frequency responses. The critical difference is the impulse response — FIR has exactly 8 non-zero samples (finite), while IIR decays but never truly reaches zero (infinite). Part (b) is also **4× more computationally efficient.**

---

*Zagazig University · Faculty of Engineering · CSE Department*

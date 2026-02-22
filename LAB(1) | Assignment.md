# 📡 DSP Lab Report — LAB(1) | Assignment

> **Digital Signal Processing** · Zagazig University · Faculty of Engineering · CSE Department

| Field | Details |
|---|---|
| **Student** | Abdelrahman Husseini Abdallah Khalil |
| **Year** | 3rd Year — Lab Section |
| **Tool** | GNU Octave |
| **Topic** | IIR Filter Design & Audio Processing |

---

## 📋 Table of Contents

| # | Section | Description |
|---|---|---|
| [Q1](#q1--butterworth-bandpass-filter) | Butterworth Bandpass Filter | Design for n=4 and n=21 |
| [S1](#s1--stability-analysis) | Stability Analysis | Pole-Zero Check |
| [Q2](#q2--whistlewav-noise-rejection) | whistle.wav Noise Rejection | Band-Stop Filter Design |
| [Q3](#q3--fir-vs-iir-moving-average-filters) | FIR vs IIR Moving Average | Comparison & Comments |
| [C3](#c3--comments-on-the-plots) | Comments on the Plots | Analysis & Conclusions |

---

## Q1 — Butterworth Bandpass Filter

> Design a Butterworth Bandpass filter for **n=4** and **n=21**. Plot Frequency Response and Impulse Response with labelled axes. Check stability.

**Key Rule:** Normalise frequency to Nyquist → `W = f / (Fs/2)`

### Order n = 4 — Complete Code

```matlab
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

### Order n = 21 — Complete Code

```matlab
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

### Observations

| Property | n = 4 | n = 21 |
|---|---|---|
| Transition band | Wide — gradual rolloff | Very narrow — sharp rolloff |
| Stopband rejection | Moderate | Very strong |
| Impulse response | Short | Long |
| Computation cost | Low | High |

---

## S1 — Stability Analysis

> **Stability Condition:** All poles must lie **inside the unit circle** in the z-plane.

Butterworth is an **IIR** filter. MATLAB/Octave `butter()` always places poles inside the unit circle by design. Verify practically with `zplane(b,a)` — all `×` marks must be inside the circle.

```matlab
zplane(b,a)   % all x (poles) inside unit circle => STABLE
```

| Check | Result |
|---|---|
| **n = 4 Stable?** | ✅ YES — all poles inside unit circle |
| **n = 21 Stable?** | ✅ YES — all poles inside unit circle |
| **Reason** | Butterworth design guarantees poles inside unit circle |

---

## Q2 — whistle.wav Noise Rejection

> Read `whistle.wav`, identify noise frequencies via FFT, design two Butterworth Band-Stop filters to reject them, verify stability, compare energy before/after, play the result.

### Signal Properties

| Parameter | Value |
|---|---|
| **Fs** | 24,000 Hz |
| **N (samples)** | 41,472 |
| **Duration** | 1.7280 s |
| **Noise peaks** | 500 Hz and 1500 Hz |

### Complete Solution — Full Code

```matlab
pkg load signal

%% -- Read the audio file
[x, fs] = audioread('D:\whistle.wav');
n    = length(x);
time = n / fs;
t    = linspace(0, time, n);
fprintf('Number of samples : %d\n', n);
fprintf('Duration          : %.4f sec\n', time);
fprintf('Sampling Freq     : %d Hz\n\n', fs);

%% -- Play original signal
sound(x, fs, 16);

%% -- Frequency spectrum of x
X = fft(x);
f = (0:n-1) * (fs / n);
figure;
plot(f(1:n/2), abs(X(1:n/2)));
xlabel('Frequency (Hz)');  ylabel('|X(f)|');
title('Spectrum of Original Signal x');  grid on;
% Notice the two sharp peaks at 500 Hz and 1500 Hz

%% -- Filter 1: Butterworth Band-Stop to reject 500 Hz
fc1 = 480;  fc2 = 520;  order = 4;
[b1, a1] = butter(order, [fc1 fc2] / (fs/2), 'stop');
[H1, w1] = freqz(b1, a1, 1024, fs);
figure; plot(w1, abs(H1));
xlabel('Frequency (Hz)');  ylabel('|H(f)|');
title('Frequency Response - Butterworth Band-Stop @ 500 Hz');  grid on;
figure; impz(b1, a1, 100);
xlabel('Samples (n)');  ylabel('Amplitude');
title('Impulse Response - Band-Stop @ 500 Hz');  grid on;
poles1  = roots(a1);  stable1 = all(abs(poles1) < 1);
fprintf('Filter 1 stable: %d\n\n', stable1);

%% -- Filter 2: Butterworth Band-Stop to reject 1500 Hz
fc3 = 1450;  fc4 = 1550;  order2 = 4;
[b2, a2] = butter(order2, [fc3 fc4] / (fs/2), 'stop');
[H2, w2] = freqz(b2, a2, 1024, fs);
figure; plot(w2, abs(H2));
xlabel('Frequency (Hz)');  ylabel('|H(f)|');
title('Frequency Response - Butterworth Band-Stop @ 1500 Hz');  grid on;
figure; impz(b2, a2, 100);
xlabel('Samples (n)');  ylabel('Amplitude');
title('Impulse Response - Band-Stop @ 1500 Hz');  grid on;
poles2  = roots(a2);  stable2 = all(abs(poles2) < 1);
fprintf('Filter 2 stable: %d\n\n', stable2);

%% -- Apply both filters (cascade)
y = filter(b1, a1, x);
y = filter(b2, a2, y);

%% -- Energy calculation
Ex = sum(x .^ 2);
Ey = sum(y .^ 2);
fprintf('Energy of original signal x : %.4f\n', Ex);
fprintf('Energy of filtered signal y : %.4f\n', Ey);
fprintf('Energy removed              : %.4f  (%.2f%%)\n\n', ...
        Ex - Ey, (Ex - Ey) / Ex * 100);

%% -- Frequency spectrum of filtered signal y
Y = fft(y);
figure; plot(f(1:n/2), abs(Y(1:n/2)));
xlabel('Frequency (Hz)');  ylabel('|Y(f)|');
title('Spectrum of Filtered Signal y');  grid on;

%% -- Play filtered signal
y_play = y / max(abs(y));
sound(y_play, fs, 16);
```

### Results

| Parameter | Value |
|---|---|
| **Filter 1 stopband** | 480–520 Hz (rejects 500 Hz) — STABLE ✅ |
| **Filter 2 stopband** | 1450–1550 Hz (rejects 1500 Hz) — STABLE ✅ |
| **Pole magnitudes** | All < 1 (inside unit circle) |
| **Energy original x** | 585.8219 |
| **Energy filtered y** | 324.8403 |
| **Energy removed** | 260.9817 (44.55%) |
| **Whistle after filter** | Gone — not audible ✅ |

---

## Q3 — FIR vs IIR Moving Average Filters

> For each filter plot the Frequency Response and Impulse Response with labelled axes and titles. Comment on the plots.

### Difference Equations

| Part | Equation |
|---|---|
| **(a)** | `y[n] = 1/8 ( x[n] + x[n-1] + x[n-2] + x[n-3] + x[n-4] + x[n-5] + x[n-6] + x[n-7] )` |
| **(b)** | `y[n] = 1/8 x[n]  -  1/8 x[n-8]  +  y[n-1]` |

---

### Part (a) — FIR Moving Average LPF

Non-recursive filter. Equal weight `1/8` for 8 samples. All `a_k = 0` except `a_0`.

```matlab
a = 1;   Fs = 8000;
b = 0.125 * [1 1 1 1 1 1 1 1];   % FIR coefficients

f = linspace(0, Fs/2, 1024);
H = freqz(b, a, f, Fs);
figure; plot(f, abs(H));
title('Frequency Response - FIR Moving Average LPF (n=7)');
xlabel('Frequency (Hz)');  ylabel('Magnitude |H(f)|');  grid on;

figure; impz(b, a);
title('Impulse Response - FIR Moving Average LPF (n=7)');
xlabel('Samples (n)');  ylabel('Amplitude');  grid on;
```

| Property | Value |
|---|---|
| **Type** | FIR — Non-recursive |
| **b vector** | `[0.125  0.125  0.125  0.125  0.125  0.125  0.125  0.125]` |
| **Stable** | Always stable — FIR filters are unconditionally stable ✅ |
| **Impulse** | Exactly 8 non-zero samples then zero (FINITE) |

---

### Part (b) — IIR Recursive Moving Average

Recursive equivalent of (a). Same averaging result with only **2 multiplications** instead of 8.

```matlab
b = [1/8, 0, 0, 0, 0, 0, 0, 0, -1/8];
a = [1, -1];
Fs = 8000;

[H, f] = freqz(b, a, 1024, Fs);
figure; plot(f, abs(H));
title('Frequency Response - IIR Recursive Moving Average');
xlabel('Frequency (Hz)');  ylabel('Magnitude |H(f)|');  grid on;

figure; impz(b, a, 50);
title('Impulse Response - IIR Recursive Moving Average');
xlabel('Samples (n)');  ylabel('Amplitude');  grid on;
```

| Property | Value |
|---|---|
| **Type** | IIR — Recursive |
| **b vector** | `[1/8, 0, 0, 0, 0, 0, 0, 0, -1/8]` |
| **a vector** | `[1, -1]` |
| **Efficient** | Only 2 multiplications vs 8 for FIR |
| **Impulse** | Decays gradually — theoretically INFINITE |

---

## C3 — Comments on the Plots

| Property | Part (a) FIR | Part (b) IIR |
|---|---|---|
| Filter type | FIR — non-recursive | IIR — recursive |
| Impulse response | Finite: 8 non-zero samples | Infinite: decays, never truly zero |
| Freq response | LPF with ripple in stopband | Nearly identical LPF shape |
| Computation | 8 multiplications per sample | 2 multiplications per sample |
| Stability | Always stable (unconditional) | Marginally stable (pole at z=1) |
| Phase | Linear phase | Non-linear phase |
| **Conclusion** | **Simple but costly at high order** | **Efficient recursive shortcut** |

> **Key conclusion:** Both parts produce the same LPF moving average and nearly identical frequency responses. The critical difference is the impulse response — FIR has exactly 8 non-zero samples (finite), while IIR decays but never truly reaches zero (infinite). Part (b) is also **4× more computationally efficient**.

---

*Zagazig University · Faculty of Engineering · CSE Department*

# 📡 Active IRS-Aided MU-MISO System Optimization
**Performance Analysis and Joint Beamforming Framework in MATLAB**

---

### 📄 Overview  
This repository provides the **MATLAB implementation** of the work presented in the paper:  

> **“Performance Analysis of Active IRS Assisted MU-MISO System”**  
> *Gitartha Gogoi, Tanushree Bose, Pradeep Vishwakarma, Samarendra Nath Sur*  
> *International Journal of Information Technology (Springer, 2025)*  
> DOI: [10.1007/s41870-025-02663-1](https://doi.org/10.1007/s41870-025-02663-1)

This project focuses on the **performance modeling and optimization** of an **Active Intelligent Reflecting Surface (IRS)-aided Multi-User Multiple-Input Single-Output (MU-MISO)** wireless communication system.  
It introduces a **joint beamforming and precoding framework** to maximize the **sum-rate** while minimizing **power consumption**, under both *weak* and *strong* line-of-sight (LoS) channel conditions.

---

### ⚙️ Key Objectives  
- Model the impact of **active amplification** and **phase shift quantization** in IRS-assisted systems.  
- Formulate a **joint optimization problem** to enhance **sum-rate** under **BS and IRS power constraints**.  
- Apply **Alternating Optimization (AO)**, **Weighted Minimum Mean Square Error (WMMSE)**, and **Fractional Programming (FP)** techniques for efficient convergence.  
- Evaluate system performance in comparison with **Passive IRS**, **Ideal Active IRS**, and **Practical Active IRS** setups.

---

### 🧠 System Description  
The system considers:
- A **multi-user MISO downlink** network with *M* BS antennas and *K* single-antenna users.  
- An **N-element Active IRS** placed between the BS and users to overcome multiplicative fading.  
- Each active IRS element integrates **reflection-type amplifiers** providing signal gain and phase control.

#### Key signal models:
- **Ideal Active IRS:**
  \[
  \Phi = \text{diag}(q_1 e^{j\theta_1}, \dots, q_N e^{j\theta_N})
  \]
- **Practical Active IRS:**
  \[
  \theta_n \in \{0, \Delta\theta, 2\Delta\theta, \dots, (L-1)\Delta\theta\}
  \]
  where \(L\) represents quantization levels and \(q_n > 0\) is the amplification factor.

---

### 🧩 Optimization Framework  

#### 🎯 **Objective Function**
Maximize the total system **sum-rate**:
\[
\max_{v, \Phi} \; R_{\text{sum}} = \sum_{k=1}^{K} \log_2(1 + \eta_k)
\]
subject to:
- BS transmit power constraint  
- Active IRS power constraint  
- Hardware power consumption inclusion

#### 🧮 **Techniques Used**
- **Fractional Programming (FP):** Reformulates non-convex logarithmic and fractional terms.  
- **Alternating Optimization (AO):** Iteratively updates BS beamforming and IRS reflection matrices.  
- **WMMSE (Weighted Minimum Mean Square Error):** Optimizes user-level SINR iteratively for convergence.  
- **QCQP (Quadratically Constrained Quadratic Programming):** Used in optimizing the IRS precoding matrix.

#### 🔁 **Joint Beamforming–Precoding Algorithm**
1. Initialize \(v\), \(\Phi\), and auxiliary variables.  
2. Optimize \(\eta_k\) and auxiliary parameters using FP.  
3. Update BS beamforming vector \(v\) via **Lagrange Multiplier method**.  
4. Update IRS precoding \(\Phi\) using **QCQP** under power constraints.  
5. Repeat until \(R_{\text{sum}}\) convergence.

---

### 🧪 Simulation Setup  

| Parameter | Value / Description |
|------------|---------------------|
| Bandwidth | 5 GHz |
| BS antennas (M) | 4 |
| IRS elements (N) | 512 |
| Channel model | Ricean fading and Rayleigh Fading |
| Power budget | 10 dBm total |
| Noise power | −100 dBm |
| Quantization levels | 2–4 bits |
| Scenarios | Weak and Strong LoS |
| Algorithms | AO, WMMSE, FP, QCQP |

The simulation compares:
1. **Without IRS**
2. **Passive IRS**
3. **Ideal Active IRS**
4. **Practical Active IRS**

---

### 📈 Results Summary  

#### 📡 **Sum-Rate Performance**
| Case | Model | Sum-Rate (bps/Hz) | Gain |
|------|--------|------------------:|------|
| Weak LoS | Without IRS | 3.02 | – |
| Weak LoS | Passive IRS | 14.10 | 4.66× |
| Weak LoS | Practical Active IRS | **16.70** | **5.52×** |
| Weak LoS | Ideal Active IRS | 32.97 | – |
| Strong LoS | Passive IRS | 20.60 | – |
| Strong LoS | Practical Active IRS | **22.27** | 1.22× |

#### 🔋 **Energy Efficiency Improvement**
- Achieves same throughput with **3–5 dBm lower total power** compared to passive IRS.  
- Indicates a **10–15% reduction** in total energy consumption.  
- Demonstrates better EE–SR tradeoff due to active amplification and optimized precoding.

#### 🧮 **Quantization Analysis**
- Increasing quantization levels from **2-bit to 4-bit** yields:  
  - **+5.02 bps/Hz** gain under weak LoS  
  - **+1.38 bps/Hz** gain under strong LoS  
- Hardware imperfections limit ideal performance, but practical gains remain significant.

---

### 🗂️ Repository Structure

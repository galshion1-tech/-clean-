# OAM-OFDM 通导一体化论文 — 公式全拆解

> 论文: A Novel Underwater Integrated Communication and Positioning Algorithm Based on OAM-OFDM
> Liu, Zhao, Zou — HIT Weihai + Dalian UT — ICASSP 2026

---

## 1. 子载波分配开关 u[m] (Eq.1)

$$u[m] = \begin{cases} 1 & \text{第 } m \text{ 个子载波用于通信} \\ 0 & \text{第 } m \text{ 个子载波用于定位} \end{cases}$$

**直觉**：OFDM 有 53 个子载波 (Nc=53)，DC 子载波 (m=0) 分配给定位，其余 52 个分配给通信。

**为什么这么做**：
- DC 子载波频率最低，传播衰减最小，适合定位参考信号
- 定位只需要一个子载波 → 不会显著降低通信吞吐量
- 通信和定位在频域正交 → 互不干扰

**论文没有说明的**：u[m] 是在发射端就固定分配好了的，接收端也知道这个分配表（"定位参考信号已知"的前提）。

---

## 2. OFDM 基带信号 (Eq.2)

$$s_{\text{OFDM}}(t) = \underbrace{\sum_{n_s=1}^{N_s}\sum_{m=0}^{N_c-1} \text{Re}\left(a_{n_s,m} \cdot e^{j2\pi(f_c + m\Delta f)t}\right) \cdot u[m] \cdot \text{rect}\left(\frac{t - T_0 - (k-1)T}{T}\right)}_{\text{通信部分}} + \underbrace{\sum_{n_s=1}^{N_s}\sum_{m=0}^{N_c-1} \text{Re}\left(b_{n_s,m} \cdot e^{j2\pi(f_c + m\Delta f)t}\right) \cdot \bar{u}[m] \cdot \text{rect}\left(\frac{t - T_0 - (k-1)T}{T}\right)}_{\text{定位部分}}$$

**符号表**：

| 符号 | 含义 | 典型值 |
|------|------|--------|
| $N_s$ | 每帧 OFDM 符号数 | 未给出 |
| $N_c$ | 子载波总数 | 53 |
| $f_c$ | 中心频率 (DC 子载波) | 未给出 |
| $\Delta f$ | 子载波间隔 | $B/N_c$ |
| $a_{n_s,m}$ | 通信星座点 (BPSK/QPSK) | $\in \{\pm1\}$ 或 $\in \{\pm1\pm j\}$ |
| $b_{n_s,m}$ | 定位星座点 (固定 BPSK) | $\pm1$ |
| $T$ | OFDM 符号长度 | $1/\Delta f + T_{CP}$ |
| $T_0$ | 起始时刻 | — |
| $\bar{u}[m]$ | $1 - u[m]$ | 定位子载波指示 |

**关键设计选择**：$b_{n_s,m}$ 使用固定 BPSK 调制。意味着：
- 定位参考信号是**已知的确定性信号** — 这是 TOA 互相关估计的前提
- 若 $b_{n_s,m}$ 也携带未知数据 → TOA 估计变成盲估计 → 精度大幅下降

**rect() 窗**：$$\text{rect}\left(\frac{t - T_0 - (k-1)T}{T}\right) = \begin{cases} 1 & T_0+(k-1)T \leq t \leq T_0+kT \\ 0 & \text{otherwise} \end{cases}$$

每个 OFDM 符号持续 T 秒，第 k 个符号的起始时刻是 $T_0 + (k-1)T$。

---

## 3. OAM 螺旋相位因子 (Eq.3)

$$\psi_n = l \cdot \phi_n$$

其中：
- $\phi_n = \frac{2\pi(n-1)}{N}$ — 第 n 个阵元的方位角（UCA 均匀分布）
- $l \in \mathbb{Z}$ — OAM 模式数（拓扑荷）
- $N$ — UCA 阵元数

**物理意义**：
- $l = 0$：所有阵元同相 → 普通波束，无涡旋
- $l = 1$：相位绕一圈变化 $2\pi$ → 单螺旋波前
- $l = -2$：反向双螺旋
- $l = 3$：三螺旋

**正交性**：$\int_0^{2\pi} e^{jl_1\phi} \cdot e^{-jl_2\phi} \, d\phi = 2\pi \cdot \delta_{l_1,l_2}$

不同模式的 OAM 波在理想自由空间中正交 → 可以在同一频率同时传输多个独立数据流。

**论文用了三个模式**：{1, -2, 3} → 传输速率 ≈ 3× 单模式 OFDM。

---

## 4. UCA 单阵元远场 (Eq.4)

$$E_n = \frac{e^{-jk d_n}}{d_n}$$

- $k = 2\pi/\lambda$ — 波数
- $d_n$ — 第 n 个阵元到远场观测点的距离

这是**球面波传播的格林函数**（三维 Helmholtz 方程的基本解）。幅度衰减 $1/d_n$，相位旋转 $e^{-jk d_n}$。

负号 $e^{-jk d_n}$ 表示波向外传播（outgoing wave）。

---

## 5. OAM 远场总声场 (Eq.5 → Eq.7)

### 5.1 叠加求和 (Eq.5)

$$E_{\text{sum}} = \sum_{n=1}^{N} E_n \cdot e^{j l \phi_n} = \sum_{n=1}^{N} \frac{e^{-jk d_n}}{d_n} \cdot e^{j l \phi_n}$$

每个阵元的信号先乘 OAM 相位因子 $e^{j l \phi_n}$，再空间叠加。

### 5.2 几何近似 (Eq.6)

$$d \approx d_n - R_t \sin\theta \cos(\varphi - \phi_n)$$

**几何推导**（论文省略了这一步）：

```
设 UCA 中心在原点，UCA 在 xy 平面：
  第 n 个阵元位置：r_n = (R_t cos φ_n, R_t sin φ_n, 0)
  远场点方向单位向量：û = (sinθ cosφ, sinθ sinφ, cosθ)

则 d_n = |r_far - r_n| ≈ r_far - û · r_n  （远场近似，泰勒展开一阶）

û · r_n = (sinθ cosφ)(R_t cosφ_n) + (sinθ sinφ)(R_t sinφ_n)
        = R_t sinθ (cosφ cosφ_n + sinφ sinφ_n)
        = R_t sinθ cos(φ - φ_n)

所以 d_n ≈ d - R_t sinθ cos(φ - φ_n)
```

注意：$d$ 是 UCA 中心到远场点的距离。论文把等号写成 $d = d_n - R_t \sin\theta \cos(\varphi - \phi_n)$，实际上应该是 $d_n = d - R_t \sin\theta \cos(\varphi - \phi_n)$，这里有个符号小问题但不影响后续推导。

### 5.3 近似代入 (Eq.7 推导)

将 Eq.6 代入 Eq.5：

$$E(d,\theta,\varphi) \approx \sum_{n=1}^{N} \frac{e^{-jk[d - R_t \sin\theta \cos(\varphi - \phi_n)]}}{d} \cdot e^{j l \phi_n}$$

幅度项 $d_n \approx d$（远场，$d \gg R_t$）。相位项不能近似，因为 $k R_t$ 不是小量。

提取公共因子：

$$E(d,\theta,\varphi) \approx \frac{e^{-jk d}}{d} \sum_{n=1}^{N} e^{j l \phi_n} \cdot e^{j k R_t \sin\theta \cos(\varphi - \phi_n)}$$

### 5.4 Bessel 函数的出现

令 $\xi = k R_t \sin\theta$，则被积式的和式为：

$$S = \sum_{n=1}^{N} e^{j l \phi_n} \cdot e^{j \xi \cos(\varphi - \phi_n)}$$

当 $N \to \infty$（连续圆环近似），用积分替换求和：

$$S \approx \frac{N}{2\pi} \int_0^{2\pi} e^{j l \alpha} \cdot e^{j \xi \cos(\varphi - \alpha)} \, d\alpha$$

令 $\beta = \alpha - \varphi + \pi/2$（旋转坐标系使余弦变成正弦），经过 Jacobi-Anger 展开：

$$e^{j \xi \cos\theta} = \sum_{m=-\infty}^{\infty} j^m J_m(\xi) \cdot e^{j m \theta}$$

代入积分后，正交性 $\int_0^{2\pi} e^{j(l+m)\beta} d\beta = 2\pi \delta_{m,-l}$ 只保留 $m = -l$ 项：

$$\int_0^{2\pi} e^{j l \alpha} \cdot e^{j \xi \cos(\varphi-\alpha)} d\alpha = 2\pi \cdot j^{-l} \cdot J_{-l}(\xi) \cdot e^{j l \varphi}$$

利用 $J_{-l}(\xi) = (-1)^l J_l(\xi)$ 和 $j^{-l} = (-j)^l$：

$$S = N \cdot j^l \cdot J_l(k R_t \sin\theta) \cdot e^{j l \varphi}$$

所以：

$$\boxed{E(d,\theta,\varphi) = \frac{N \cdot e^{-jk d}}{d} \cdot j^l \cdot J_l(k R_t \sin\theta) \cdot e^{j l \varphi}}$$

**关键物理含义**：
- $J_l(k R_t \sin\theta)$ — Bessel 包络，决定了涡旋波束的**环形强度分布**（中心为零，环状峰值）
- $e^{j l \varphi}$ — 螺旋相位波前，携带 OAM
- $j^l$ — 常数相位偏置，对应论文 Eq.7 中的 $j^l$

论文的 Eq.7 写的是 $N \cdot \frac{e^{-jkd}}{d} \cdot j^l \cdot J_l(k R_t \sin\theta) \cdot e^{jl\varphi}$，一致。

---

## 6. OAM-OFDM 发射信号 (Eq.8)

$$s_{\text{OAM-OFDM}}(t) = \sum_{n_s=1}^{N_s}\sum_{m=0}^{N_c-1} \text{Re}\left(a_{n_s,m} \cdot e^{j 2\pi(f_c + m\Delta f)t + j l \phi_n}\right) \cdot u[m] \cdot \text{rect}(\cdots) + \sum_{n_s=1}^{N_s}\sum_{m=0}^{N_c-1} \text{Re}\left(b_{n_s,m} \cdot e^{j 2\pi(f_c + m\Delta f)t + j l \phi_n}\right) \cdot \bar{u}[m] \cdot \text{rect}(\cdots)$$

**Eq.8 相对于 Eq.2 的唯一变化**：每个阵元的信号多乘了一个 $e^{j l \phi_n}$。

对发射 UCA 的第 n 个阵元：
$$x_n(t) = s_{\text{OFDM}}(t) \cdot e^{j l \phi_n}$$

全部 N 个阵元同时发射不同相位的同一 OFDM 信号 → 空间叠加形成涡旋波束。

**物理实现**：每个阵元用独立的 DAC + 功放，数字域乘以 $e^{j l \phi_n}$ 再上变频。

---

## 7. 接收信号模型 (Eq.9)

$$S_r(t) = U \cdot s_{\text{OAM-OFDM}}(t - \tau) * h(t) + n(t)$$

- $U$ — 传播损耗（标量，常数衰减）
- $\tau$ — 传播时延（待估计的 TOA）
- $h(t)$ — 多径信道冲激响应
- $n(t)$ — 复 AWGN
- $*$ — 卷积

**这是标准的水声信道模型**：大尺度衰减 $U$ + 时延 $\tau$ + 多径 $h(t)$ + 噪声。

**论文假设**：$h(t)$ 是已知或可估计的。实际上水声多径 $h(t)$ 的估计本身就是一个难题，论文没讨论。

---

## 8. AOA 估计 — UCA 阵列流形 (Eq.10)

$$a(\theta_k, \varphi_k) = \begin{bmatrix} e^{j \frac{2\pi R_t}{\lambda_k} \sin\theta_k \cos(\varphi_k - \phi_1)} \\ e^{j \frac{2\pi R_t}{\lambda_k} \sin\theta_k \cos(\varphi_k - \phi_2)} \\ \vdots \\ e^{j \frac{2\pi R_t}{\lambda_k} \sin\theta_k \cos(\varphi_k - \phi_N)} \end{bmatrix}$$

**物理含义**：第 k 个信号从方向 $(\theta_k, \varphi_k)$ 入射到 UCA，第 n 个阵元相对于 UCA 中心的相位差是：

$$\Delta \psi_n = -\frac{2\pi}{\lambda_k} \cdot (\text{阵元 n 到参考点的程差}) = -\frac{2\pi}{\lambda_k} \cdot R_t \sin\theta_k \cos(\varphi_k - \phi_n)$$

注意这里的负号被约掉了（因为论文定义的 $a(\theta_k,\varphi_k)$ 实际上用了 $+j$ 而非 $-j$），这是因为 UCA 接收时的时间反演约定不同。

**与你的 ULA 导向矢量对比**：

| | ULA (你的项目) | UCA (这篇论文) |
|---|---|---|
| 阵元排列 | 沿 y 轴等距 d | 圆周等角距 $2\pi/N$ |
| 导向矢量 | $a_m = e^{-j k y_m \sin\theta}$ | $a_n = e^{j k R_t \sin\theta \cos(\varphi-\phi_n)}$ |
| 自由度 | 1D DOA (水平角 θ) | 2D DOA (仰角 θ + 方位角 φ) |
| 阵列流形结构 | Vandermonde | 非 Vandermonde → 需要 beamspace 变换 |

**UCA 的核心困难**：$a(\theta_k, \varphi_k)$ 不是 Vandermonde 矩阵 → 不能直接用 root-MUSIC（root-MUSIC 要求阵列流形是 Vandermonde）。

---

## 9. Beamspace 变换 (Eq.12)

$$y(t) = F_r \cdot x(t) = B \cdot s(t) + F_r \cdot n(t)$$

**核心思想**：通过西变换 $F_r$ 将 UCA 阵元空间映射到"相位模式空间"（phase mode domain）。

### 9.1 $F_r$ 的构造

$F_r$ 是一个 $N \times (2K+1)$ 的波束形成矩阵（实际通常用 DFT 矩阵的变体），使得变换后的导向矢量变为：

$$b(\theta_k, \varphi_k) = F_r^H \cdot a(\theta_k, \varphi_k)$$

变换后 $b(\theta_k, \varphi_k)$ 具有近似 Vandermonde 结构。

### 9.2 为什么可以这样做

UCA 的导向矢量可以展开为 Fourier 级数（利用 Bessel 函数的 Jacobi-Anger 展开）：

$$e^{j \xi \cos(\varphi - \phi_n)} = \sum_{m=-\infty}^{\infty} j^m J_m(\xi) \cdot e^{j m(\varphi - \phi_n)}$$

其中 $\xi = k R_t \sin\theta$。

$F_r$ 的每一行对应一个相位模式 m → 阵元空间 → 模式空间的变换：

$$[F_r]_{m,n} = \frac{1}{\sqrt{N}} e^{j m \phi_n}$$

实际实现时取 $m \in \{-K, \dots, K\}$，其中 $K = \lfloor k R_t \rfloor$（Bessel 函数的有效模式数）。

### 9.3 论文省略的关键细节

- $K$ 的选取：$J_m(\xi)$ 当 $|m| > \xi$ 时迅速衰减 → 取 $K \approx k R_t$ 足够
- Beamspace 变换后维度从 $N \times 1$ 变成 $(2K+1) \times 1$ — 降维了
- 这个变换是频率依赖的 → 不同频率需要不同的 $F_r$ → 对 OFDM 多载波是潜在问题

---

## 10. root-MUSIC AOA 估计 (Eq.13-17)

### 10.1 复变量表示 (Eq.13)

定义：
$$\xi_k = \frac{2\pi R_t \sin\theta_k \cos\varphi_k}{\lambda_k}, \quad \psi_k = \frac{2\pi R_t \sin\theta_k \sin\varphi_k}{\lambda_k}$$

则导向矢量变为：

$$a_k(z_{1k}, z_{2k}) = \begin{bmatrix} z_{1k}^{\cos\phi_1} \cdot z_{2k}^{\sin\phi_1} \\ z_{1k}^{\cos\phi_2} \cdot z_{2k}^{\sin\phi_2} \\ \vdots \\ z_{1k}^{\cos\phi_N} \cdot z_{2k}^{\sin\phi_N} \end{bmatrix}$$

其中 $z_{1k} = e^{j\xi_k}$，$z_{2k} = e^{j\psi_k}$。

经过 beamspace 变换后，每一行对应一个相位模式 m，导向矢量变成：

$$b_k(z_{1k}, z_{2k}) \propto \begin{bmatrix} z_{1k}^{-K} z_{2k}^{-K} \\ \vdots \\ z_{1k}^0 z_{2k}^0 \\ \vdots \\ z_{1k}^{K} z_{2k}^{K} \end{bmatrix}$$

→ **现在具有 Vandermonde 结构了！**

### 10.2 协方差矩阵与特征分解 (Eq.14)

$$R_y = E[y(t) y^H(t)] = U_s \Sigma_s U_s^H + U_n \Sigma_n U_n^H$$

- $U_s$：信号子空间（前 p 个特征向量，p = 信号数）
- $U_n$：噪声子空间（后 $(2K+1)-p$ 个特征向量）
- 关键性质：$U_s \perp U_n$，且 $b_k \in \text{span}(U_s)$ → $b_k^H U_n = 0$

### 10.3 二元多项式求根 (Eq.15)

$$P(z_1, z_2) = a^H(z_1, z_2) \cdot U_n U_n^H \cdot a(z_1, z_2) = 0$$

在单位圆上 $(z_1, z_2) = (e^{j\xi}, e^{j\psi})$ 请求根 → 得到 $(\xi_k, \psi_k)$ 对。

**为什么是二元多项式**：因为 UCA 的 2D DOA 需要两个参数 → 两个根。

**实现细节**（论文省略）：
- 实际上 $P(z_1, z_2)$ 展开后是 $\sum_{p,q} c_{pq} z_1^p z_2^q$ 形式
- 对 $z_1$ 而言它是多项式（$z_2$ 固定时），可以用 root-finding
- 通常做法：交替固定一个求另一个 → 迭代收敛

### 10.4 角度反解 (Eq.16-17)

$$\varphi_k = \arctan\left(\frac{\psi_k}{\xi_k}\right)$$

$$\theta_k = \arcsin\left(\frac{\lambda_k \sqrt{\xi_k^2 + \psi_k^2}}{2\pi R_t}\right)$$

**推导**：

从定义：
$$\xi_k = \frac{2\pi R_t}{\lambda_k} \sin\theta_k \cos\varphi_k$$
$$\psi_k = \frac{2\pi R_t}{\lambda_k} \sin\theta_k \sin\varphi_k$$

比值得 $\tan\varphi_k = \psi_k/\xi_k$ → Eq.16。

平方和得 $\xi_k^2 + \psi_k^2 = (2\pi R_t/\lambda_k)^2 \sin^2\theta_k$ → Eq.17。

## 11. 解螺旋恢复 OFDM (Eq.19)

$$S_r'(t) = S_r(t) \cdot e^{-j l \phi_1}$$

**为什么 $\phi_1$**：接收 UCA 的第 1 个阵元处的信号作为参考，乘 $e^{-j l \phi_1}$ 去掉螺旋相位。

**物理解释**：
- 发射时每个阵元乘 $e^{j l \phi_n^{\text{tx}}}$ 形成涡旋
- 接收时乘 $e^{-j l \phi_1^{\text{rx}}}$ 解涡旋
- 接收 UCA 的 $\phi_1^{\text{rx}}$ 不一定等于发射 UCA 的 $\phi_n^{\text{tx}}$ → 解螺旋不完美 → TOA 估计有残留误差

**论文省略的问题**：接收端解螺旋只用了第一个阵元的相位 $\phi_1$，但信号到达每个接收阵元时 OAM 相位分布已经因为传播距离和方向发生了变化。完整的解螺旋应该按接收 UCA 的每个阵元独立做 $S_r^{(n)}(t) \cdot e^{-j l \phi_n^{\text{rx}}}$。

---

## 12. 互相关 TOA 估计 (Eq.20-21)

### 12.1 定位参考信号 (Eq.18)

$$S_p(t) = e^{j 2\pi f_c t}$$

只在 DC 子载波上调制 BPSK +1 的恒包络信号。

### 12.2 离散互相关 (Eq.20)

$$R[m] = \sum_{j=0}^{T \cdot f_s} S_r'[j] \cdot S_p^*[j - m]$$

- $S_r'[j]$：解螺旋后的接收信号（已下变频到基带）
- $S_p^*[j-m]$：本地参考信号的共轭延迟版本
- 互相关在 $m = \tau \cdot f_s$ 处取得峰值

### 12.3 峰值检测 (Eq.21)

$$\hat{\tau} = \frac{1}{f_s} \cdot \arg\max_m |R[m]|$$

TOA 分辨率 = $1/f_s$。如果只有 DC 子载波 → 带宽极小 → TOA 分辨率很差。

**这可能是论文中定位误差 12m 的根本原因**：
- DC 子载波带宽 ≈ 0 → TOA 分辨率由符号长度限制 → 距离分辨率 ≈ $c \cdot T \approx 1500 \times (\text{几ms})$ → 米级甚至十米级
- 如果用宽带 LFM 做 TOA → 距离分辨率 ≈ $c/B \approx 1500/5000 = 0.3\text{m}$

---

## 13. 3D 位置重建 (Eq.22)

$$\begin{cases} x_1 = \hat{d} \cdot \sin\hat{\theta}_1 \cdot \cos\hat{\varphi}_1 \\ y_1 = \hat{d} \cdot \sin\hat{\theta}_1 \cdot \sin\hat{\varphi}_1 \\ z_1 = \hat{d} \cdot \cos\hat{\theta}_1 \end{cases}$$

其中 $\hat{d} = c \cdot \hat{\tau}$。

**球坐标 → 直角坐标**，以接收 UCA 中心为原点。这是标准的远场定位几何。

---

## 14. 公式链总图

```
                  发射端                                接收端
                   ====                                ====

  ┌──────────────────────────┐           ┌──────────────────────────┐
  │ 通信数据 a_{ns,m}        │           │                          │
  │ 定位数据 b_{ns,m}        │           │  UCA 接收 x(t) = A s(t)  │
  │      ↓                   │           │       ↓                  │
  │ 子载波分配 u[m]          │           │  ┌─────────┴─────────┐   │
  │      ↓                   │           │  ↓                   ↓   │
  │ OFDM 调制 (IFFT)         │           │ beamspace F_r       e^{-jlφ₁}
  │ s_OFDM(t)               │           │  ↓                   ↓   │
  │      ↓                   │           │ b(θ,φ)            S_r'(t) │
  │ 乘 OAM 相位 e^{jlφ_n}    │           │  ↓                   ↓   │
  │      ↓                   │           │ R_y = E[yy^H]      互相关│
  │ UCA N阵元同时发射         │           │  ↓                   ↓   │
  │ s_OAM-OFDM(t)           │  ═══════► │ U_n 噪声子空间    τ̂ = argmax│
  └──────────────────────────┘  水下信道  │  ↓                   ↓   │
                                         │ P(z₁,z₂)=0        d̂ = cτ̂ │
                                         │  ↓                   ↓   │
                                         │ (ξ_k, ψ_k) → (θ̂, φ̂)      │
                                         │  └─────────┬─────────┘   │
                                         │            ↓             │
                                         │  (x̂, ŷ, ẑ) = d̂(sinθ̂cosφ̂, │
                                         │           sinθ̂sinφ̂, cosθ̂) │
                                         └──────────────────────────┘
```

---

## 15. 各步骤的 MATLAB 可实现性

| 步骤 | 难度 | 所需工具箱 | 备注 |
|------|------|-----------|------|
| OFDM 调制 | 低 | 无 (手写 IFFT) | 与你的 LFM 生成类似 |
| OAM 相位叠加 | 低 | 无 | `exp(1j*l*phi_n)` |
| 多径信道 | 中 | 无 | 用你已有的虚源法 |
| UCA beamspace | 中 | 无 | 需要实现 $F_r$ 矩阵 |
| root-MUSIC | 高 | 无 | 二元多项式求根复杂，可用 2D MUSIC 谱搜索替代 |
| 解螺旋 | 低 | 无 | 乘 $e^{-jl\phi_1}$ |
| 互相关 TOA | 低 | xcorr | MATLAB 内置 |
| 全流程仿真 | 高 | — | 整合 OFDM + OAM + UCA + 多径 + MUSIC |

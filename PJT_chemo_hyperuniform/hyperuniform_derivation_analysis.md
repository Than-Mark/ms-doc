# Hyperuniformity理论推导的完整分析
## 从Boltz & Ihle的反平行活性模型到Chemotactic模型的系统推导与对比

---

## 目录

1. [引言：问题与方法](#1-引言)
2. [Boltz & Ihle模型：从Master方程到Structure Factor](#2-boltz--ihle模型)
3. [为什么Boltz & Ihle能得到涨落抑制？核心机制分析](#3-核心机制分析)
4. [Chemotactic模型的标准推导](#4-chemotactic模型的标准推导)
5. [对Chemotactic模型的改进尝试](#5-对chemotactic模型的改进尝试)
6. [综合对比与结论](#6-综合对比与结论)
7. [附录](#7-附录)

---

## 1. 引言

### 1.1 什么是Hyperuniformity？

超均匀（hyperuniform）态是一类特殊的点构型，其大尺度密度涨落被异常抑制。对于通常的流体，在一个半径为 $R$ 的窗口内计数的粒子数 $N$ 的方差满足 $\langle \delta N^2 \rangle \sim R^d$（$d$ 为空间维度），即与体积成正比。而在超均匀态中：

$$\langle \delta N^2 \rangle \sim R^{\beta}, \quad \text{其中 } \beta < d$$

等价地，用静态结构因子 $S(\mathbf{q}) = \langle \delta\rho_{\mathbf{q}} \delta\rho_{-\mathbf{q}} \rangle$ 来描述：

$$S(\mathbf{q} \to 0) \to 0 \quad \text{（超均匀）}$$

对比之下，泊松（随机）点过程的 $S(0) = \rho_0$（平均密度），普通流体的 $S(0)$ 为有限正值。

### 1.2 Boltz & Ihle (PRE 2025) 的发现

Boltz & Ihle研究了一维晶格上的反平行活性粒子模型，发现长波密度涨落被抑制。他们的核心方法和技术路线是：

- **模型**：双组分粒子（向右/左运动），通过反平行相互作用改变方向
- **方法**：Poisson Representation——将离散Master方程转化为连续Fokker-Planck方程
- **关键发现**：出现了虚噪声（imaginary noise），导致非泊松统计
- **结果**：$S(k \to 0) = \rho_0 - c^2/(2D) < \rho_0$

### 1.3 本文的目标

本文要完成三项任务：

**任务一**：完整复现Boltz & Ihle的推导——不是仅仅陈述结果，而是展示每一步的数学细节，特别是解释Poisson Representation的替换规则为什么成立。

**任务二**：提炼出Boltz & Ihle成功的"核心机制"——不是描述性的，而是结构性的：什么样的数学结构导致了涨落抑制。

**任务三**：将这个核心机制应用到Chemotactic模型上——系统地尝试多种可能的推导路径，分析为什么某些路径行不通，以及如果要实现hyperuniformity需要什么样的结构修改。

---

## 2. Boltz & Ihle模型：从Master方程到Structure Factor

### 2.1 模型定义

#### 2.1.1 晶格与粒子

考虑一维晶格，位点标记为 $i = 1, 2, ..., L$（周期边界条件）。每个位点上有两种粒子：
- $n_i^+ \in \{0, 1, 2, ...\}$：向右（$+$）运动的粒子数
- $n_i^- \in \{0, 1, 2, ...\}$：向左（$-$）运动的粒子数

#### 2.1.2 动力学规则

**Streaming（流动）**：$+$ 粒子以速率 $\lambda$ 从 $i$ 跳到 $i+1$；$-$ 粒子以速率 $\lambda$ 从 $i$ 跳到 $i-1$。

**Conversion（转化）**：$+$ 粒子以速率 $r_+(n_i^+, n_i^-)$ 变为 $-$ 粒子；$-$ 粒子以速率 $r_-(n_i^+, n_i^-)$ 变为 $+$ 粒子。

**反平行相互作用**（核心设计）：

$$r_\pm(n^+, n^-) = r_0 + (n^\pm - 1)r_1$$

其中 $r_0$ 是自发转化速率，$r_1$ 是相互作用驱动的转化速率。物理含义：局部有更多同向粒子时，它们更频繁地碰撞并翻转方向，抑制净流动。

#### 2.1.3 Master方程

概率分布 $P(\{n_i^+, n_i^-\}, t)$ 满足：

$$\frac{\partial P}{\partial t} = S_{\text{out}}^\lambda + C_{\text{out}}^r + S_{\text{in}}^\lambda + C_{\text{in}}^r$$

**流出项**：

$$S_{\text{out}}^\lambda = -\lambda\sum_i(n_i^+ + n_i^-)P$$

$$C_{\text{out}}^r = -\sum_i[r_+(n_i^+, n_i^-)n_i^+ + r_-(n_i^+, n_i^-)n_i^-]P$$

**流入项（Streaming）**：

$$S_{\text{in}}^\lambda = \lambda\sum_i\left[(n_{i-1}^+ + 1)P(..., n_{i-1}^++1, n_i^+-1, ...) + (n_{i+1}^- + 1)P(..., n_{i+1}^-+1, n_i^--1, ...)\right]$$

**流入项（Conversion）**：

$$C_{\text{in}}^r = \sum_i\left[(n_i^+ + 1)r_+(n_i^++1, n_i^--1)P(..., n_i^++1, n_i^--1, ...) + (n_i^- + 1)r_-(n_i^+-1, n_i^-+1)P(..., n_i^+-1, n_i^-+1, ...)\right]$$

### 2.2 Poisson Representation

现在面临的核心问题是：Master方程定义在离散的整数格点上，包含非线性项（如 $n_i^+ n_i^-$），直接求解极其困难。我们需要一种系统的方法将离散方程转化为连续方程。

#### 2.2.1 为什么选择Poisson核？

Poisson分布 $\Pi(n|\lambda) = e^{-\lambda}\lambda^n/n!$ 有一个独特的性质——**位移恒等式**：

$$n \cdot \frac{e^{-\alpha}\alpha^n}{n!} = \frac{e^{-\alpha}\alpha^n}{(n-1)!} = \alpha \cdot \frac{e^{-\alpha}\alpha^{n-1}}{(n-1)!}$$

这意味着算子 $n$ 作用在Poisson核上等价于乘以 $\alpha$ 后指标下移一位。更关键的是，通过分部积分，离散算子可以被**精确替换**为连续微分算子。

正式定义Poisson Representation：将离散概率分布表示为连续伪分布 $f$ 的积分变换：

$$\boxed{P(\{n_i^+, n_i^-\}) = \int \prod_i \left[\frac{e^{-\alpha_i}\alpha_i^{n_i^+}}{n_i^+!} \frac{e^{-\beta_i}\beta_i^{n_i^-}}{n_i^-!}\right] f(\{\alpha_i, \beta_i\}) \, d\alpha_i d\beta_i}$$

其中 $\alpha_i, \beta_i$ 是连续变量（通常是复值），$f$ 是伪分布（可以取负值或复值）。

#### 2.2.2 算子替换规则的严格证明

**规则一：$n_i^+ P \to (\alpha_i - \partial_{\alpha_i}\alpha_i)f$**

我们需要证明：
$$\int d\alpha \, \frac{e^{-\alpha}\alpha^n}{n!} \cdot n \cdot g(\alpha) = \int d\alpha \, \frac{e^{-\alpha}\alpha^n}{n!} \cdot (1-\partial_\alpha)\alpha \cdot g(\alpha)$$

展开右边：
$$\text{RHS} = \underbrace{\int d\alpha \, \frac{e^{-\alpha}\alpha^{n+1}}{n!} g}_{\text{term A}} - \underbrace{\int d\alpha \, \frac{e^{-\alpha}\alpha^n}{n!} \partial_\alpha(\alpha g)}_{\text{term B}}$$

对term B分部积分（假设边界项为零）：
$$\text{term B} = -\int d\alpha \, \partial_\alpha\left(\frac{e^{-\alpha}\alpha^n}{n!}\right) \alpha g = \int d\alpha \left[\frac{e^{-\alpha}\alpha^n}{n!} - \frac{e^{-\alpha}\alpha^{n-1}}{(n-1)!}\right] \alpha g$$

$$= \text{term A} - n \int d\alpha \frac{e^{-\alpha}\alpha^n}{n!} g$$

因此：
$$\text{RHS} = \text{term A} - [\text{term A} - n \cdot \text{(积分)}] = n \int d\alpha \frac{e^{-\alpha}\alpha^n}{n!} g = \text{LHS}$$

证毕。$$\boxed{n_i P \to (1 - \partial_{\alpha_i})\alpha_i f = (\alpha_i - \partial_{\alpha_i}\alpha_i)f}$$

**规则二：$(n_{i-1}^+ + 1)P(..., n_{i-1}^++1, n_i^+-1, ...) \to (\alpha_{i-1} - \partial_{\alpha_i}\alpha_{i-1})f$**

流入项的Poisson核为 $\pi(n_{i-1}^++1|\alpha_{i-1})\pi(n_i^+-1|\alpha_i)$。利用 $(n+1)\pi(n+1|\alpha) = \alpha\pi(n|\alpha)$：

$$(n_{i-1}^+ + 1) \cdot \pi(n_{i-1}^++1|\alpha_{i-1}) \cdot \pi(n_i^+-1|\alpha_i) = \alpha_{i-1} \cdot \pi(n_{i-1}^+|\alpha_{i-1}) \cdot \frac{n_i^+}{\alpha_i} \pi(n_i^+|\alpha_i)$$

用规则一替换 $n_i^+$，对 $\alpha_i$ 分部积分，最终得到 $(\alpha_{i-1} - \partial_{\alpha_i}\alpha_{i-1})f$。

**规则三：$n_i^+(n_i^+-1)P \to (1-\partial_{\alpha_i})^2\alpha_i^2 f$**

两次应用规则一即得。

#### 2.2.3 为什么 $f$ 可以是复值？

Poisson表示是一个**积分变换**，不是概率分布的重新参数化。只要最终物理量 $P(\{n_i\})$ 及其矩是实的和非负的，$f$ 可以取任何值。

虚部编码了**非泊松统计**。例如，参数为 $\lambda + i\xi$ 的Poisson变量，方差为 $\lambda - \xi^2 < \lambda$（sub-Poissonian）。

### 2.3 Fokker-Planck方程

将替换规则应用到Master方程的各项中。

**Streaming贡献** $F_\lambda$：

$$F_\lambda = \lambda\sum_i\left[\partial_{\alpha_i}(\alpha_i - \alpha_{i-1}) - \partial_{\beta_i}(\beta_i - \beta_{i+1})\right]f$$

**Conversion贡献** $F_r$ 和 **噪声项** $J$：

展开 $r_\pm = r_0 + (n^\pm - 1)r_1$，只保留二阶导数项（一阶导数贡献漂移，二阶导数贡献噪声）：

$$\boxed{J = -r_1\sum_i\left[\partial_{\alpha_i}^2\alpha_i^2 + \partial_{\beta_i}^2\beta_i^2\right]f}$$

$$F_r = r_0\sum_i(\partial_{\alpha_i} - \partial_{\beta_i})(\alpha_i - \beta_i)f + r_1\sum_i(\partial_{\alpha_i} - \partial_{\beta_i})(\alpha_i + \beta_i)(\alpha_i - \beta_i)f$$

完整Fokker-Planck方程：
$$\partial_t f = F_\lambda + F_r + J$$

### 2.4 Langevin方程

#### 2.4.1 确定性部分

$$\partial_t \alpha_i|_{\text{det}} = \lambda(\alpha_{i-1} - \alpha_i) - r_0(\alpha_i - \beta_i) - r_1(\alpha_i + \beta_i)(\alpha_i - \beta_i)$$
$$\partial_t \beta_i|_{\text{det}} = \lambda(\beta_{i+1} - \beta_i) + r_0(\alpha_i - \beta_i) + r_1(\alpha_i + \beta_i)(\alpha_i - \beta_i)$$

#### 2.4.2 扩散矩阵与虚噪声

从噪声项 $J = -r_1\sum_i(\partial_{\alpha_i}^2\alpha_i^2 + \partial_{\beta_i}^2\beta_i^2)f$，对比标准Fokker-Planck形式：

$$D_{\alpha_i\alpha_i} = -2r_1\alpha_i^2, \quad D_{\beta_i\beta_i} = -2r_1\beta_i^2, \quad D_{\alpha_i\beta_i} = 0$$

#### 2.4.3 变量变换：$(\alpha, \beta) \to (\tilde\rho, \tilde m)$

引入组合变量：$\tilde\rho = \alpha + \beta$（总密度），$\tilde m = \alpha - \beta$（磁化强度）。

扩散矩阵在新变量下（Jacobian $J = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$）：

$$D' = JDJ^T = -r_1\begin{pmatrix} \tilde\rho^2 + \tilde m^2 & 2\tilde\rho\tilde m \\ 2\tilde\rho\tilde m & \tilde\rho^2 + \tilde m^2 \end{pmatrix}$$

在 $|\tilde m| \ll \tilde\rho \approx \rho_0$ 极限下：$D' \approx -r_1\rho_0^2 I$。

**虚数矩阵平方根出现了**：
$$B = \sqrt{D'} = i\rho_0\sqrt{r_1} \cdot I$$

#### 2.4.4 最终Langevin方程

确定性方程变为：
$$\boxed{\partial_t \tilde\rho = -v\nabla\tilde m}$$
$$\boxed{\partial_t \tilde m = -v\nabla\tilde\rho - r_0\tilde m - r_1\tilde\rho\tilde m + i\sqrt{r_1\tilde\rho}\,\xi}$$

核心发现：
- **密度方程无直接噪声**！噪声只进入磁化方程。
- 出现了**虚噪声** $i\sqrt{r_1\tilde\rho}\,\xi$——非泊松统计的标志。

### 2.5 Structure Factor的计算

#### 2.5.1 线性化与Fourier空间

$\tilde\rho = \rho_0 + \delta\tilde\rho$，$|\tilde m| \ll \rho_0$：

$$\partial_t \delta\tilde\rho_k = -ikv\tilde m_k$$
$$\partial_t \tilde m_k = -ikv\delta\tilde\rho_k - (r_0 + r_1\rho_0)\tilde m_k + i\sqrt{r_1\rho_0}\,\xi_k$$

#### 2.5.2 过阻尼极限

合并为二阶方程，特征值 $\lambda_+ \approx -Dk^2$，$\lambda_- \approx -(r_0+r_1\rho_0)$，其中 $D = v^2/(r_0+r_1\rho_0)$。在过阻尼极限下：

$$\partial_t\delta\tilde\rho_k = -Dk^2\delta\tilde\rho_k + ic\eta_k$$

其中 $c = v\sqrt{r_1\rho_0/(r_0+r_1\rho_0)}$，$\eta_k$ 是空间关联噪声，$\langle\eta_k\eta_{-k}\rangle = 2[1-\cos(k)] \approx k^2$。

#### 2.5.3 稳态Structure Factor

由于噪声是纯虚数（$ic\eta_k$）：

$$\langle\delta\tilde\rho_k\delta\tilde\rho_{-k}\rangle = -\frac{c^2[1-\cos(k)]}{Dk^2}$$

物理Structure Factor（$N_S = \rho_0 + $ Poisson场涨落）：

$$\boxed{N_S(k) = \rho_0 - \frac{c^2[1-\cos(k)]}{Dk^2} \approx \rho_0 - \frac{c^2}{2D} + O(k^2)}$$

关键结果：$S_0 = \rho_0(1 - \frac{r_1\rho_0}{2(r_0+r_1\rho_0)})$。在 $r_1\rho_0 \gg r_0$ 极限下 $S_0/\rho_0 \approx 1/2$。

---

## 3. 为什么Boltz & Ihle能得到涨落抑制？核心机制分析

### 3.1 密度模式的"保护"结构

观察Langevin方程的结构：

| 方程 | 噪声入口 |
|------|---------|
| $\partial_t \tilde\rho = -v\nabla\tilde m$ | **无直接噪声** |
| $\partial_t \tilde m = ... + i\sqrt{r_1\tilde\rho}\,\xi$ | 有噪声 |

噪声路径：$\xi \rightarrow \tilde m \xrightarrow{-v\nabla} \tilde\rho$

密度涨落通过磁化模式的梯度**间接**接收噪声，形成**双重过滤**：
1. 噪声进入 $\tilde m$（步骤1：局域耦合）
2. $\tilde m$ 通过 $-v\nabla\tilde m$ 影响 $\tilde\rho$（步骤2：梯度过滤）

在Fourier空间：$\text{noise}_\rho \sim k \times \text{noise}_m \sim k^2$。

### 3.2 虚噪声的非泊松统计效应

在Poisson representation中：
- 虚噪声 $i\sqrt{...}\,\xi$ 驱动 $\tilde\rho$ 的虚部涨落
- 物理密度涨落 = 泊松基线 $\rho_0$ + **虚数涨落的负贡献**
- 类比：参数 $\lambda = \rho_0 + i\xi$ 的泊松变量，方差 $= \rho_0 - \langle\xi^2\rangle < \rho_0$
- 结果：**sub-Poissonian统计**（方差小于均值）

### 3.3 与守恒律的联系

在该模型中：
- 总粒子数守恒（streaming不改变局域总数）
- 磁化（极化）不守恒（conversion改变 $m$）
- 这种"部分守恒"结构是产生超均匀涨落抑制的关键

---

## 4. Chemotactic模型的标准推导

现在，我们试图将上述分析框架应用到Chemotactic模型上。本节先建立标准模型，展示为什么它不是hyperuniform。

### 4.1 模型定义

粒子Langevin方程：
$$\frac{d\mathbf{r}_i}{dt} = \alpha\nabla c(\mathbf{r}_i, t) + \sqrt{2D}\,\boldsymbol{\xi}_i(t)$$

化学场方程：
$$\partial_t c = D_c\nabla^2c + \mu c + k\sum_{j=1}^N\delta(\mathbf{r} - \mathbf{r}_j(t))$$

参数：$\alpha = -|\alpha| < 0$（负趋化），$k = -|k| < 0$（消耗），$\mu = -|\mu| < 0$（衰变）。

### 4.2 Dean-Kawasaki随机密度方程

从粒子描述到场描述（使用Ito公式）：

$$\boxed{\partial_t\rho = -\nabla\cdot(\alpha\rho\nabla c) + D\nabla^2\rho + \nabla\cdot\left[\sqrt{2D\rho}\,\boldsymbol{\eta}\right]}$$

化学场准静态响应：$\nabla c_q = -i\mathbf{q}\frac{|k|}{D_cq^2+|\mu|}\rho_q$。

### 4.3 线性化与有效弛豫率

$$\partial_t \delta\rho_q = \frac{|\alpha||k|\rho_0 q^2}{D_c q^2 + |\mu|}\delta\rho_q - Dq^2\delta\rho_q + \sqrt{2D\rho_0}\,i\mathbf{q}\cdot\boldsymbol{\eta}_q$$

$$\boxed{\gamma(q) = Dq^2 - \frac{|\alpha||k|\rho_0 q^2}{D_c q^2 + |\mu|}}$$

**稳定性条件**：$D|\mu| > |\alpha||k|\rho_0$。

### 4.4 Structure Factor——发散！

$$S(q) = \frac{D\rho_0 q^2}{\gamma(q)^2}$$

对 $q \to 0$：$\gamma(q) \approx \frac{D|\mu| - |\alpha||k|\rho_0}{|\mu|}q^2 = O(q^2)$

$$\boxed{S(q) \sim \frac{q^2}{(q^2)^2} = O\left(\frac{1}{q^2}\right) \to \infty}$$

**不是hyperuniform。** 原因：Dean-Kawasaki噪声 $\sim q^2$，弛豫率 $\gamma(q) \sim q^4$，所以 $S \sim q^2/q^4 = 1/q^2$。

### 4.5 与Boltz & Ihle的初步对比

| 特征 | Boltz & Ihle | 标准Chemotaxis |
|------|-------------|---------------|
| 场变量 | 2个（$\tilde\rho$, $\tilde m$） | 1个（$\rho$） |
| 噪声入口 | 只进入辅助场 $\tilde m$ | 直接进入密度场 |
| 噪声过滤 | 双重（$\sim k^2$） | 单次（$\sim q$） |
| S(k→0) | 有限（$< \rho_0$） | 发散（$\sim 1/k^2$） |

关键差异：Boltz & Ihle有双场结构，密度方程无直接噪声。Chemotaxis中噪声直接进入密度方程。

---

## 5. 对Chemotactic模型的改进尝试

标准模型不是hyperuniform，但其根源是什么？本节系统地尝试四种改进方案，每种方案都针对上一方案暴露出的问题。

### 5.1 方案一：化学场独立涨落

**动机**：标准模型中噪声只来自粒子扩散（Dean-Kawasaki）。化学场由离散粒子的$\delta$函数源驱动，本身应该有shot noise。这种额外的噪声源可能改变Structure Factor的行为。

**实现**：化学场方程加入shot noise项：
$$\partial_t c = D_c\nabla^2c + \mu c + k\rho + \sqrt{|k|\rho_0}\,\eta_c$$

**分析**：化学噪声通过 $\nabla c$ 进入密度方程。Fourier空间：
- 化学噪声强度：$\Gamma_c(q) = \frac{|\alpha|^2\rho_0^3|k|\,q^4}{(D_cq^2+|\mu|)^2} = O(q^4)$（$q \to 0$）
- 扩散噪声强度：$\Gamma_\rho(q) = 2D\rho_0 q^2 = O(q^2)$（$q \to 0$）

**结果**：小q极限下扩散噪声主导，$S(q) \sim 1/q^2 \to \infty$。

**教训**：化学场的shot noise在 $q \to 0$ 时为 $O(q^4)$，被 $O(q^2)$ 的扩散噪声掩盖。需要一种噪声源在 $q \to 0$ 时比 $q^2$ 更强，或者需要一个"保护机制"让扩散噪声不直接作用于密度。

### 5.2 方案二：晶格化 + Poisson Representation

**动机**：Boltz & Ihle的成功依赖于Poisson Representation产生的虚噪声和非泊松统计。也许将Chemotactic模型放到晶格上并用Poisson表示，可以恢复类似的结构。

**实现**：一维晶格，化学场快速弛豫 $c_i \approx -|k|n_i/|\mu|$。Hopping速率：
$$\lambda_{i\to i+1} = \lambda_0[1 + \tilde{\varepsilon}(n_{i+1} - n_i)]$$

Poisson表示 $P(\{n_i\}) = \int \prod_i \frac{e^{-\alpha_i}\alpha_i^{n_i}}{n_i!} f(\{\alpha_i\}) d\alpha_i$。

**关键计算**：线性hopping不产生噪声。非线性hopping（$O(\tilde{\varepsilon})$）产生二阶导数项：

$$\boxed{J = -\lambda_0\tilde{\varepsilon}\sum_i\left[\partial_{\alpha_i}^2(\alpha_i^2) - \partial_{\alpha_i}\partial_{\alpha_{i+1}}(\alpha_i\alpha_{i+1})\right]f}$$

**分析**：
- $\partial^2(\alpha^2)$ 项 $\to$ 局部虚噪声（类似Boltz & Ihle）
- $\partial\partial(\alpha\alpha')$ 项 $\to$ 空间关联噪声

Langevin方程（线性化后Fourier空间）：
$$\partial_t \delta\tilde\rho_q = -D_{\text{eff}}q^2\delta\tilde\rho_q + i\tilde{c}\,q\sqrt{2[1-\cos(qa)]}\,\xi_q$$

Structure Factor：
$$N_S(k) = \rho_0 - \frac{\tilde{c}^2 a^2}{2D_{\text{eff}}^2 k^2}$$

**结果**：$k \to 0$ 时**发散**（趋向负无穷）！

**与Boltz & Ihle的对比**：

| 特征 | Boltz & Ihle | 晶格Chemotaxis |
|------|-------------|---------------|
| 场结构 | 双场（$\tilde\rho$, $\tilde m$） | 单场（$\tilde\rho$） |
| Streaming | $-v\nabla\tilde m$（磁化梯度驱动密度） | $-v_0\nabla\cdot(\tilde\rho\nabla\tilde\rho)$（自驱动） |
| 正则化 | 额外的过阻尼模式提供 $q^2$ 正则化 | 缺乏额外正则化 |
| S(k→0) | 有限（$\rho_0 - c^2/2D$） | 发散（$\sim -1/k^2$） |

**教训**：虽然有虚噪声，但由于缺乏Boltz & Ihle模型中的双场耦合结构，结果仍然发散。晶格化和Poisson表示本身不足以产生hyperuniformity——关键在于双场结构。

### 5.3 方案三：惯性动力学——引入独立速度场

**动机**：方案二的教训是单场结构不够。Boltz & Ihle有 $\tilde\rho$ 和 $\tilde m$ 两个独立场。如果在Chemotactic模型中引入速度 $\mathbf{v}$ 作为独立场，密度方程变为确定性的（无噪声），噪声只进入速度方程——这模仿了Boltz & Ihle的结构。

**实现**：
$$m\frac{d^2\mathbf{r}_i}{dt^2} = -\gamma\frac{d\mathbf{r}_i}{dt} + \alpha\nabla c(\mathbf{r}_i) + \sqrt{2\gamma k_BT}\,\boldsymbol{\xi}_i$$

场方程：
- 连续性方程（**无噪声**）：$\partial_t\rho + \nabla\cdot(\rho\mathbf{v}) = 0$
- 动量方程（**有噪声**）：$\partial_t\mathbf{v} = -\frac{\gamma}{m}\mathbf{v} + \frac{\alpha}{m}\nabla c + \boldsymbol{\eta}_v$

**结构对比**：
```
Boltz & Ihle:          Inertial Chemotaxis:
∂_t ρ̃ = -v∇m̃         ∂_t ρ + ∇·(ρv) = 0
∂_t m̃ = ... + noise    ∂_t v = -(γ/m)v + (α/m)∇c + η_v
(密度无噪声)            (密度无噪声)
```

**分析**：线性化后消去 $\delta\mathbf{v}$，过阻尼极限：
$$\partial_t\delta\rho_q = -D_{\text{eff}}(q)q^2\delta\rho_q + \frac{im\rho_0}{\gamma}\mathbf{q}\cdot\boldsymbol{\eta}_v$$

噪声：$\sim \mathbf{q}\cdot\boldsymbol{\eta}_v \sim q$（**单次梯度过滤**）

$$S(q) \sim \frac{q^2}{q^4} = \frac{1}{q^2} \to \infty$$

**仍然发散！**

**教训**：虽然有双场结构且密度方程无噪声，但噪声通过 $\nabla\cdot(\rho\mathbf{v})$ 只经过**一次梯度过滤**（velocity noise $\to$ density：$\sim q$）。Boltz & Ihle中是**双重过滤**（noise $\to$ $m$ $\to$ $\nabla m$ $\to$ $\rho$：$\sim k^2$）。区别在于耦合结构：Boltz & Ihle的 $\partial_t\rho = -v\nabla m$ 中密度由磁化的**梯度**驱动，而惯性Chemotaxis的 $\partial_t\rho = -\nabla\cdot(\rho\mathbf{v})$ 中密度由速度的**散度**直接驱动。

### 5.4 方案四：活性场论框架

**动机**：前三个方案从具体模型出发，但都没有成功。也许应该从更一般的理论框架来理解问题。Zheng, Klatt & Lowen (PRL 2024) 证明了一类广义的活性场论generic地导致hyperuniformity。

**通用活性场论结构**：
$$\partial_t \rho = -\nabla\cdot[\rho(v_0 \mathbf{P} + ...)] + D\nabla^2\rho + ...$$
$$\partial_t \mathbf{P} = -\frac{1}{\tau}\mathbf{P} + \frac{\lambda}{\rho}\nabla\rho + ... + \boldsymbol{\eta}_P$$

其中 $\lambda/\rho$ 项是**活性应力**，将极化场 $\mathbf{P}$（独立场）耦合到密度梯度。

**Chemotactic模型的改写**：
$$\mathbf{v} = \alpha\nabla c = \frac{\alpha k}{D_c}\nabla[(-\nabla^2 + \xi^{-2})^{-1}\rho]$$

Fourier空间：$\mathbf{v}_q \propto \frac{i\mathbf{q}\,\rho_q}{D_cq^2+|\mu|} = O(q)$。

**关键问题**：Chemotactic coupling是**非局部**的，且速度 $\mathbf{v}$ 是**完全由 $\rho$ 决定**的（slaved），不是独立场变量。这与通用活性场论有独立极化场 $\mathbf{P}$ 的结构根本不同。

**前三个方案的教训汇总**：

| 方案 | 尝试了什么 | 为什么失败 | 缺少什么 |
|------|-----------|-----------|---------|
| 化学场独立涨落 | 添加shot noise | 化学噪声 $O(q^4)$，被扩散噪声掩盖 | 噪声没有"保护"密度 |
| 晶格+Poisson | 虚噪声+非泊松统计 | 单场结构，无额外正则化 | 双场耦合结构 |
| 惯性动力学 | 速度作为独立场 | 单次梯度过滤（$\sim q$） | 双重梯度过滤（$\sim k^2$） |

**要实现hyperuniformity，需要**：
1. 一个**独立的辅助场**（接收噪声）
2. 辅助场通过**梯度耦合**影响密度：$\partial_t\rho \sim \nabla \times$（辅助场）
3. 形成**双重梯度过滤**：noise $\to$ aux $\to$ $\nabla$aux $\to$ $\rho$

---

## 6. 综合对比与结论

### 6.1 两种机制的根本物理差异

**Boltz & Ihle模型**——自限反馈：
```
密度过高 → 粒子更多碰撞 → 更多方向翻转 → 净流动被抑制 → 密度涨落被自发抑制
```

**Chemotactic模型**——正反馈不稳定：
```
密度过高 → 消耗更多化学 → 局部c降低 → 负趋化粒子向低c移动 → 向高密度区聚集 → 不稳定性
```

这两种机制在物理上是根本不同的：
- **反平行相互作用** → 自限、涨落抑制
- **负趋化+消耗** → 正反馈、聚集不稳定

### 6.2 实现Hyperuniformity的可能路径

**路径一：双场结构**（最接近Boltz & Ihle）
- 引入独立的辅助场（磁化、极化、或化学场的独立分量）
- 噪声只进入辅助场，辅助场通过梯度耦合影响密度
- 形成双重梯度过滤

**路径二：Manna-like守恒律**
- 引入质心守恒或类似约束
- 约束条件强制 $S(k \to 0) = 0$

**路径三：非线性饱和效应**
- 在非线性稳态中，chemotactic聚集可能形成有序但超均匀的图案

**路径四：多组分扩展**
- 多种粒子产生/消耗不同化学物质
- 交叉扩散效应可能抑制涨落

### 6.3 对Chemotactic模型"Hyperuniform"表述的理解

文档提到 $k < 0$ 情形"more relevant to the emergence of hyperuniformity"，可能的含义：

1. **临界点附近**：在 $D|\mu| \approx |\alpha|\rho_0$ 的临界区域，系统可能表现出特殊的标度行为
2. **非线性饱和态**：不稳定后的聚集态可能形成有序但超均匀的图案
3. **与活性场论的联系**：chemotactic coupling的结构与活性应力有形式相似性，可能通过适当修改实现hyperuniformity
4. **有效抑制**：虽然标准推导 $S \sim 1/q^2$，但在某些参数regime，chemotactic coupling可能有效降低涨落幅度

---

## 7. 附录

### 7.1 关键公式汇总

**Boltz & Ihle模型**：
- 弛豫率：$D = v^2/(r_0 + r_1\rho_0)$
- 噪声强度：$c = v\sqrt{r_1\rho_0/(r_0+r_1\rho_0)}$
- Structure Factor：$N_S(k) = \rho_0 - c^2[1-\cos(k)]/(Dk^2)$
- 极限值：$S_0 = \rho_0(1 - r_1\rho_0/[2(r_0+r_1\rho_0)])$

**Chemotactic模型（标准连续）**：
- 有效弛豫率：$\gamma(q) = Dq^2 - |\alpha||k|\rho_0 q^2/(D_cq^2+|\mu|)$
- 稳定性条件：$D|\mu| > |\alpha||k|\rho_0$
- Structure Factor：$S(q) = D\rho_0 q^2/\gamma(q)^2 \sim 1/q^2$（发散）

**晶格Chemotaxis + Poisson**：
- 噪声项：$J = -\lambda_0\tilde{\varepsilon}\sum_i[\partial_{\alpha_i}^2(\alpha_i^2) - \partial_{\alpha_i}\partial_{\alpha_{i+1}}(\alpha_i\alpha_{i+1})]f$
- 有效扩散：$D_{\text{eff}} = \lambda_0 a^2(1 - \tilde{\varepsilon}\rho_0)$
- Structure Factor：$N_S(k) = \rho_0 - \tilde{c}^2 a^2/(2D_{\text{eff}}^2 k^2)$（发散）

### 7.2 符号表

| 符号 | 含义 |
|------|------|
| $\alpha_i, \beta_i$ | Poisson representation中的复值参数 |
| $\tilde\rho = \alpha + \beta$ | 总密度Poisson参数 |
| $\tilde m = \alpha - \beta$ | 磁化强度Poisson参数 |
| $\rho_0$ | 平均粒子密度 |
| $D$ | 粒子扩散系数 |
| $D_c$ | 化学场扩散系数 |
| $\mu$ | 化学场衰变速率（$\mu < 0$） |
| $k$ | 粒子产生($k>0$)/消耗($k<0$)化学物的速率 |
| $\alpha$ | 趋化敏感度（$\alpha>0$吸引，$\alpha<0$排斥） |
| $r_0, r_1$ | 自发/相互作用驱动转化速率 |
| $v = \lambda a$ | 有效速度 |
| $\lambda_0$ | Hopping速率 |
| $\varepsilon, \tilde{\varepsilon}$ | Chemotactic偏置参数 |
| $c$ | 有效噪声强度 |
| $N_S(k)$ / $S(q)$ | 静态结构因子 |
| $\xi = \sqrt{D_c/|\mu|}$ | 化学场屏蔽长度 |

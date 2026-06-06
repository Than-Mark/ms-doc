# Chemotactic active matter 的 hyperuniform 性质：沿 Boltz--Ihle 方法论的理论分析

本文尝试把 Boltz 与 Ihle 在反对齐活性粒子模型中使用的分析路线，迁移到当前 `chemo_hyperuniform.pdf` 中的趋化模型。核心目标不是直接宣称系统是否 hyperuniform，而是问：

1. 从微观随机动力学出发，密度涨落的噪声结构是什么？
2. 线性化后小波数结构因子 \(S(q)\) 如何随 \(q\to 0\) 变化？
3. 与 Boltz--Ihle 模型相比，哪些结构保留了，哪些结构缺失了？

结论先写在前面：当前 PDF 中的标准 overdamped chemotactic 模型，一般不产生严格 hyperuniformity。在线性稳定的均匀相中，长波结构因子趋向一个有限常数。采用 Dean-Kawasaki 规范并只看 connected density fluctuation 时，

\[
S_{\rm DK}(q\to 0)
=
\frac{D\rho_0}{D+\alpha\rho_0/|\mu|}.
\]

因此 \(S(0)\neq 0\)。当 \(\alpha>0\) 且化学物被消耗时，趋化反馈提高有效扩散，可压低长波涨落，但只得到 finite suppression；当 \(\alpha<0\) 且化学物被消耗时，反馈降低有效扩散，稳定区中涨落反而增强，并在长波不稳定阈值处发散。

这与 Boltz--Ihle 的一维反对齐模型很相似的一点是：二者的严格解析结论都不是自动给出 \(S(0)=0\)。但差异也很关键：Boltz--Ihle 至少通过反对齐生成了 sub-Poissonian 的负涨落贡献，而标准趋化模型的普通扩散噪声直接进入密度方程，缺乏把 \(S(0)\) 推到零的保护机制。

---

## 1. Boltz--Ihle 可迁移的方法论骨架

Boltz--Ihle 的理论路线可以压缩为：

\[
\text{微观动力学}
\to
\text{主方程}
\to
\text{Poisson representation}
\to
\text{Fokker--Planck/Langevin}
\to
\text{线性涨落方程}
\to
S(k).
\]

真正有用的不是 Poisson representation 这个技术标签本身，而是它揭示了三件事：

1. 噪声必须从微观过程推导，而不是随手添加。
2. 判断 hyperuniformity 要看 \(S(k\to0)\)，不能只看确定性稳定性。
3. 小 \(k\) 行为由“弛豫率”和“噪声强度”共同决定。

对一个线性 Langevin 模式

\[
\partial_t \delta\rho_{\mathbf q}
=-\gamma(q)\delta\rho_{\mathbf q}
\zeta_{\mathbf q},
\]

若

\[
\langle \zeta_{\mathbf q}(t)\zeta_{-\mathbf q}(t')\rangle
=2A(q)\delta(t-t'),
\]

则稳态结构因子的基本尺度是

\[
S(q)=\frac{A(q)}{\gamma(q)}.
\]

所以：

- 若 \(A(q)\sim q^2\)、\(\gamma(q)\sim q^2\)，则 \(S(q\to0)\sim \text{const.}\)，普通非 hyperuniform。
- 若 \(A(q)\sim q^4\)、\(\gamma(q)\sim q^2\)，则 \(S(q)\sim q^2\)，严格 hyperuniform。
- 若 \(A(q)\sim q^2\)、\(\gamma(q)\to \gamma_0>0\)，则 \(S(q)\sim q^2\)，也可得到 hyperuniform-like 抑制。

Boltz--Ihle 的模型本质上通过极化场的间接噪声路径和 imaginary noise 产生了 sub-Poissonian 修正，但其一维解析结果仍是

\[
S(k)=S_0+A k^2+\cdots,
\qquad S_0>0.
\]

这正是本文分析 chemotactic 模型时要继承的警惕：有限小 \(k\) 上的下降趋势不等于 \(S(0)=0\)。

---

## 2. 当前 PDF 中的 chemotactic 模型

PDF 的粒子模型为

\[
\dot{\mathbf r}_i
=
\alpha\nabla c(\mathbf r_i,t)
+\sqrt{2D}\,\boldsymbol\xi_i(t),
\]

\[
\partial_t c
=D_c\nabla^2 c+\mu c
+\sum_j k\,\delta(\mathbf r-\mathbf r_j).
\]

后文重点讨论消耗情形。按 PDF 的无量纲写法，可写成

\[
\partial_t c
=D_c\nabla^2 c+\mu c-\rho,
\qquad \mu=-|\mu|<0,
\]

其中

\[
\rho(\mathbf r,t)=\sum_i\delta(\mathbf r-\mathbf r_i(t)).
\]

对应的连续平均场方程是

\[
\partial_t\rho
=-\alpha\nabla\cdot(\rho\nabla c)
+D\nabla^2\rho,
\]

\[
\partial_t c
=D_c\nabla^2 c+\mu c-\rho.
\]

若要分析 hyperuniformity，必须恢复密度噪声。由 overdamped Brownian 粒子的 Dean-Kawasaki 方程，线性噪声形式为

\[
\partial_t\rho
=-\alpha\nabla\cdot(\rho\nabla c)
+D\nabla^2\rho
+\nabla\cdot\left[\sqrt{2D\rho}\,\boldsymbol\eta\right].
\]

这一步是 Boltz--Ihle 方法论中“不能只看确定性方程”的对应物。

---

## 3. 均匀态的线性确定性结构

令

\[
\rho=\rho_0+\delta\rho,
\qquad
c=c_0+\delta c.
\]

线性化得到

\[
\partial_t\delta\rho
=
-\alpha\rho_0\nabla^2\delta c
+D\nabla^2\delta\rho
+\nabla\cdot\left[\sqrt{2D\rho_0}\,\boldsymbol\eta\right],
\]

\[
\partial_t\delta c
=
D_c\nabla^2\delta c
-|\mu|\delta c
-\delta\rho.
\]

在 Fourier 空间中：

\[
\partial_t
\begin{pmatrix}
\delta\rho_{\mathbf q}\\
\delta c_{\mathbf q}
\end{pmatrix}
=
\begin{pmatrix}
-Dq^2 & \alpha\rho_0 q^2\\
-1 & -|\mu|-D_cq^2
\end{pmatrix}
\begin{pmatrix}
\delta\rho_{\mathbf q}\\
\delta c_{\mathbf q}
\end{pmatrix}
+\zeta_{\rho,\mathbf q}.
\]

其中密度噪声满足

\[
\zeta_{\rho,\mathbf q}
=i\mathbf q\cdot \sqrt{2D\rho_0}\,\boldsymbol\eta_{\mathbf q},
\]

\[
\langle \zeta_{\rho,\mathbf q}(t)
\zeta_{\rho,-\mathbf q}(t')\rangle
=2D\rho_0 q^2\delta(t-t').
\]

注意这个噪声是守恒噪声，因此有一个 \(q^2\) 因子。但这只够保证普通流体的 \(S(0)\) 有限，并不足以保证 hyperuniform。

---

## 4. 消去化学场后的有效密度方程

若化学场比密度快，或者只关心长波低频极限，可以令化学场准静态响应：

\[
\delta c_{\mathbf q}
\simeq
-\frac{\delta\rho_{\mathbf q}}{|\mu|+D_cq^2}.
\]

代入密度方程：

\[
\partial_t\delta\rho_{\mathbf q}
=
-\gamma(q)\delta\rho_{\mathbf q}
+\zeta_{\rho,\mathbf q},
\]

其中

\[
\gamma(q)
=
q^2
\left[
D+\frac{\alpha\rho_0}{|\mu|+D_cq^2}
\right].
\]

小 \(q\) 极限为

\[
\gamma(q)
=D_{\rm eff}q^2+O(q^4),
\]

\[
D_{\rm eff}
=
D+\frac{\alpha\rho_0}{|\mu|}.
\]

因此均匀态的长波稳定条件是

\[
D_{\rm eff}>0.
\]

也就是

\[
D+\frac{\alpha\rho_0}{|\mu|}>0.
\]

若 \(\alpha<0\)，这给出 PDF 中得到的负趋化长波不稳定阈值：

\[
D|\mu|>|\alpha|\rho_0.
\]

---

## 5. 结构因子：标准模型不是严格 hyperuniform

线性 OU 过程给出

\[
S_{\rm DK}(q)
=
\frac{D\rho_0 q^2}{\gamma(q)}.
\]

代入 \(\gamma(q)\)：

\[
S_{\rm DK}(q)
=
\frac{D\rho_0}
{
D+\dfrac{\alpha\rho_0}{|\mu|+D_cq^2}
}.
\]

因此

\[
\boxed{
S_{\rm DK}(q\to0)
=
\frac{D\rho_0}
{D+\alpha\rho_0/|\mu|}
}
\]

只要系统在线性稳定区，这就是有限正数，而不是零。

这给出三个物理情形：

1. \(\alpha>0\)：消耗型化学场下，密度高处使 \(c\) 降低，正趋化粒子被吸向高 \(c\) 区，因此离开高密区。反馈提高有效扩散，

\[
D_{\rm eff}>D,
\]

于是

\[
S(0)<\rho_0.
\]

这是一种长波涨落抑制，但不是严格 hyperuniform。

2. \(\alpha<0\)：负趋化粒子偏向低 \(c\)，而高密区消耗更多化学物，恰好形成低 \(c\) 区。反馈降低有效扩散，

\[
D_{\rm eff}<D,
\]

于是

\[
S(0)>\rho_0.
\]

这不是涨落抑制，而是长波涨落增强。

3. 当 \(\alpha<0\) 且

\[
D|\mu|\to |\alpha|\rho_0
\]

时，

\[
D_{\rm eff}\to 0^+,
\]

线性结构因子发散。此时系统靠近长波聚集不稳定，而不是靠近 hyperuniform。

所以，如果当前 PDF 关注的是 \(\alpha<0,k<0\) 的 “chemorepulsion + consumption”，它的线性机制更像 positive feedback/clustering，而不是 Boltz--Ihle 式的 self-limiting fluctuation suppression。

---

## 6. 与 Boltz--Ihle 的结构性对比

Boltz--Ihle 的有效线性结构是

\[
\partial_t\delta\rho_k
=-D_B k^2\delta\rho_k
+i c\,\eta_k,
\]

其中

\[
\langle \eta_k\eta_{-k}\rangle
\sim 1-\cos k
\sim k^2.
\]

同时，由 Poisson representation 的 imaginary noise，物理结构因子中出现低于 Poisson 的负修正：

\[
S(k)
=
\rho_0
-\frac{c^2}{D_Bk^2}(1-\cos k)
=S_0+Ak^2+\cdots.
\]

chemotactic 模型的有效线性结构则是

\[
\partial_t\delta\rho_q
=-D_{\rm eff}q^2\delta\rho_q
+i\mathbf q\cdot\sqrt{2D\rho_0}\boldsymbol\eta_q.
\]

因此

\[
S(q\to0)=\frac{D\rho_0}{D_{\rm eff}}.
\]

两者的关键差异是：

| 结构 | Boltz--Ihle 反对齐模型 | 标准趋化模型 |
|---|---|---|
| 慢变量 | 总密度 \(\rho\) | 总密度 \(\rho\) |
| 辅助快变量 | 极化/磁化 \(m\) | 化学场 \(c\)，但通常被 \(\rho\) slaved |
| 噪声入口 | 主要进入 \(m\)，且是 imaginary noise | 普通 Brownian 噪声直接以守恒形式进入 \(\rho\) |
| 长波阻尼 | \(\sim k^2\) | \(\sim q^2\) |
| 噪声强度 | 有 sub-Poissonian 负修正 | \(+D\rho_0 q^2\)，普通正噪声 |
| 结果 | \(S(0)\) 有限但低于 Poisson | \(S(0)\) 有限；可低于或高于 Poisson，取决于反馈符号 |
| 严格 hyperuniform | 否，一维解析中 \(S_0>0\) | 否，标准稳定相中 \(S(0)>0\) |

---

## 7. 为什么标准 chemotaxis 很难直接给出 \(S(0)=0\)

要得到严格 hyperuniform，线性层面至少需要

\[
\frac{A(q)}{\gamma(q)}
\to0.
\]

对当前模型：

\[
A(q)\sim D\rho_0q^2,
\qquad
\gamma(q)\sim D_{\rm eff}q^2.
\]

二者同阶，所以比值是常数。换句话说，Dean-Kawasaki 守恒噪声只提供了一次梯度过滤，而扩散型密度弛豫本身也只在 \(q^2\) 阶发生。两者抵消后留下有限 \(S(0)\)。

若想使 \(S(q)\to0\)，需要下面任一结构：

1. 噪声更强地被长波过滤：

\[
A(q)\sim q^4,
\qquad
\gamma(q)\sim q^2.
\]

2. 密度模式有非扩散型长波恢复力：

\[
\gamma(q)\to \gamma_0>0,
\qquad
A(q)\sim q^2.
\]

3. 存在类似 Poisson representation 中 imaginary noise 的负涨落贡献，精确抵消 Poisson 基线：

\[
\rho_0-\Delta S(0)=0.
\]

标准 overdamped chemotaxis 三者都不满足。

---

## 8. 化学场自身噪声是否能改变结论？

若给化学场加入产生/消耗的 shot noise：

\[
\partial_t c
=D_c\nabla^2c-|\mu|c-\rho
+\sqrt{\sigma_c}\,\eta_c,
\]

那么该噪声通过 \(-\alpha\rho_0\nabla^2\delta c\) 作用到密度方程。准静态下，

\[
\delta c_{\mathbf q}^{\rm noise}
\sim
\frac{\sqrt{\sigma_c}}{|\mu|+D_cq^2}\eta_{c,\mathbf q},
\]

所以密度方程中的化学噪声幅度为

\[
\zeta_{c,\mathbf q}^{\rho}
\sim
\alpha\rho_0 q^2
\frac{\sqrt{\sigma_c}}{|\mu|+D_cq^2}
\eta_{c,\mathbf q}.
\]

对应强度

\[
A_c(q)\sim q^4.
\]

这部分单独看会给

\[
S_c(q)\sim \frac{q^4}{q^2}\sim q^2,
\]

确实是 hyperuniform-like 的。但它在小 \(q\) 下被 Brownian 扩散噪声

\[
A_\rho(q)\sim q^2
\]

压过。除非粒子扩散噪声被某种机制消除或显著压低，否则化学场 shot noise 不能决定 \(S(q\to0)\)。

---

## 9. 当前模型可能出现的“表观 hyperuniform”来源

虽然标准线性理论不给出严格 hyperuniform，但数值上仍可能看见 \(S(q)\) 在有限区间下降。可能来源包括：

1. 正趋化加消耗时 \(D_{\rm eff}>D\)，使 \(S(0)\) 小于 Poisson 基线；若 \(S_0\) 很小，有限 \(q\) 区间可伪装成幂律。
2. 化学场 shot noise 给出 \(q^2\) 型贡献，但只在 Brownian 噪声不主导的中间波数区间可见。
3. 非线性聚集或图案形成后，点阵结构可能产生额外的几何有序性；这已超出均匀态线性涨落理论，需要对非线性稳态重新展开。
4. 有限尺寸下最小 \(q=2\pi/L\) 不够小，无法分辨 \(S(q)=S_0+Aq^2\) 与 \(S(q)\sim q^\alpha\)。

这正呼应 Boltz--Ihle 的方法论提醒：不要把有限尺度幂律外观直接等同于 \(S(0)=0\)。

---

## 10. 若要让 chemotactic 模型真正趋向 hyperuniform，可能需要的结构修改

沿 Boltz--Ihle 的思路，最自然的修改不是单纯调参，而是改变噪声入口和场耦合结构。

### 10.1 抑制直接密度噪声

标准 Brownian 扩散给出

\[
\nabla\cdot\sqrt{2D\rho}\,\eta,
\]

这是 \(A(q)\sim q^2\) 的来源。若 \(D\to0\)，并且主要噪声来自化学场，则可能得到

\[
S(q)\sim q^2.
\]

但这不是普通 Brownian particle 模型，而是接近确定性运动加化学随机性的极限。

### 10.2 引入独立极化场

设

\[
\partial_t\rho=-v_0\nabla\cdot\mathbf p,
\]

\[
\partial_t\mathbf p
=-\Gamma\mathbf p
-\chi\nabla\rho
+\text{noise}.
\]

若噪声只进入 \(\mathbf p\)，且密度只通过 \(\nabla\cdot\mathbf p\) 感受它，则有机会形成类似 Boltz--Ihle 的间接噪声路径。但仍需检查最终 \(A(q)\) 是否达到 \(q^4\)，以及是否存在普通密度扩散噪声残留。

### 10.3 加入自限反馈而不是正反馈

Boltz--Ihle 的反对齐是“高密度/高极化 \(\to\) 更强阻尼极化 \(\to\) 抑制输运”。而 \(\alpha<0,k<0\) 的 chemotaxis 是“高密度 \(\to\) 更低化学浓度 \(\to\) 粒子继续进入高密区”。后者天然增强长波涨落。

如果希望趋化模型呈现涨落抑制，应考虑让反馈符号变成自限型，例如正趋化加消耗，或负趋化加生产。在线性层面这会提高 \(D_{\rm eff}\)，但仍只得到有限 \(S(0)\)；若再配合噪声保护机制，才可能接近严格 hyperuniform。

---

## 11. 最终判断

对当前 PDF 中的标准模型，线性涨落理论给出的判断是：

\[
\boxed{
\text{稳定均匀相中一般不是严格 hyperuniform，因为 }S(q\to0)>0.
}
\]

更细分地说：

- \(\alpha>0,k<0,\mu<0\)：长波涨落被压低，可能是 near-hyperuniform 或 sub-Poisson-like suppression，但 \(S(0)\) 有限。
- \(\alpha<0,k<0,\mu<0\)：长波涨落增强，靠近阈值时发散，并导致长波聚集不稳定。
- 化学场 shot noise 单独可以产生 \(q^2\) 型结构因子贡献，但通常被 Brownian 密度噪声主导。
- 若数值中观察到 \(S(q)\sim q^\alpha\)，应检查是否其实是 \(S_0+Aq^2\) 在有限区间的表观幂律。

因此，沿 Boltz--Ihle 的理论精神，当前 chemotactic 模型最重要的结论不是“它已经 hyperuniform”，而是：

\[
\text{它提供了一个可分析的反馈模型，}
\]

\[
\text{但标准 overdamped 版本缺少把长波密度涨落严格压到零的噪声保护结构。}
\]

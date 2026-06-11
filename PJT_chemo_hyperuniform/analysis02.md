可以分析，而且你的数值结果“正的小 $\mu$ 也能看到类似 hyperuniform 的涨落抑制”并不一定和线性理论矛盾。关键在于：

$$\mu>0$$

在**无限大系统、严格 $k\to0$** 下确实会带来长波不稳定；但在**有限尺寸数值模拟**中，只要 $\mu$ 小到低于最小非零波数对应的稳定阈值，所有可见的非零 Fourier 模式仍然可以稳定，并且可以表现出很强的涨落抑制。

下面详细推导。

你 PDF 里的 continuum model 是

$$\partial_t\rho = -\alpha\nabla\cdot(\rho\nabla c)+D\nabla^2\rho,$$

$$\partial_t c = D_c\nabla^2c+\mu c-\rho .$$

文档中也说明后续主要关注 consumption / negative-chemotaxis 更相关的情形，但数学上 $\mu>0$ 当然也可以分析。

---

# 1. 均匀态

设

$$\rho(\mathbf r,t)=\rho_0, \qquad c(\mathbf r,t)=c_0.$$

因为均匀态下

$$\nabla\rho_0=0, \qquad \nabla c_0=0, \qquad \nabla^2 c_0=0,$$

密度方程自动满足。

化学场方程给出

$$0=\mu c_0-\rho_0.$$

所以

$$\boxed{ c_0=\frac{\rho_0}{\mu}. }$$

如果 $\mu>0$ 很小，那么

$$c_0$$

会很大。但由于粒子运动只依赖

$$\nabla c,$$

均匀的 $c_0$ 本身不会直接产生趋化力。

不过这也提示一个数值细节：如果模拟中没有把 $c$ 的均值设在 $c_0=\rho_0/\mu$，那么 $c$ 的空间平均会按照

$$\frac{d\bar c}{dt} = \mu\bar c-\rho_0$$

演化；当 $\mu>0$ 时，这个均值方向是指数不稳定的。但由于它是空间均匀模式，不产生 $\nabla c$，对粒子动力学可能几乎不可见，除非数值上溢出或模型其他项依赖 $c$ 的绝对值。

---

# 2. 线性化方程

令

$$\rho=\rho_0+\delta\rho, \qquad c=c_0+\delta c.$$

加入 Dean 型守恒噪声后，线性化密度方程是

$$\partial_t\delta\rho = D\nabla^2\delta\rho + \alpha\rho_0\nabla^2\delta c + \eta,$$

其中

$$\eta = -\nabla\cdot \left[ \sqrt{2D\rho_0}\boldsymbol{\Lambda} \right].$$

化学场方程线性化为

$$\partial_t\delta c = D_c\nabla^2\delta c+\mu\delta c-\delta\rho.$$

Fourier 变换后：

$$\partial_t\delta\rho_{\mathbf k} = -Dk^2\delta\rho_{\mathbf k} + \alpha\rho_0 k^2\delta c_{\mathbf k} + \eta_{\mathbf k},$$

$$\partial_t\delta c_{\mathbf k} = -(D_ck^2-\mu)\delta c_{\mathbf k} - \delta\rho_{\mathbf k}.$$

写成矩阵形式：

$$\partial_t \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix} = A(k) \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix} + \begin{pmatrix} \eta_{\mathbf k}\\ 0 \end{pmatrix},$$

其中

$$A(k)= \begin{pmatrix} -Dk^2 & \alpha\rho_0k^2\\ -1 & \mu-D_ck^2 \end{pmatrix}.$$

---

# 3. 均匀态失稳条件

令线性增长率为 $s$，则

$$s$$

是矩阵 $A(k)$ 的本征值。

特征方程是

$$s^2-(\mathrm{Tr}A)s+\det A=0.$$

其中

$$\mathrm{Tr}A = \mu-(D+D_c)k^2,$$

$$\det A = (-Dk^2)(\mu-D_ck^2) + (\alpha\rho_0k^2)(-1).$$

所以

$$\det A = k^2 \left[ \alpha\rho_0-D\mu+DD_ck^2 \right].$$

对于二维线性系统，稳定条件是

$$\boxed{ \mathrm{Tr}A<0, }$$

$$\boxed{ \det A>0. }$$

因此：

$$\mu-(D+D_c)k^2<0,$$

即

$$\boxed{ \mu<(D+D_c)k^2. }$$

以及

$$\alpha\rho_0-D\mu+DD_ck^2>0,$$

即

$$\boxed{ \mu<D_ck^2+\frac{\alpha\rho_0}{D}. }$$

所以某个波数 $k$ 的稳定条件是

$$\boxed{ \mu< \min \left[ (D+D_c)k^2, D_ck^2+\frac{\alpha\rho_0}{D} \right]. }$$

---

# 4. 为什么 $\mu>0$ 在无限系统中一定有问题？

在无限大系统中，允许

$$k\to0.$$

对任意正的 $\mu$，总可以找到足够小的 $k$，使得

$$(D+D_c)k^2<\mu.$$

于是

$$\mathrm{Tr}A>0,$$

线性系统失稳。

所以在严格热力学极限下：

$$\boxed{ \mu>0 \quad\Rightarrow\quad k\to0 \text{ 长波模式不稳定。} }$$

这就是我之前说 $\mu>0$ 时不能直接谈稳态

$$S(k\to0)$$

的原因。

但这还不是数值模拟看到的全部故事。

---

# 5. 有限尺寸系统中的稳定范围

数值模拟中系统大小有限，周期边界下最小非零波数是

$$k_{\min}=\frac{2\pi}{L}.$$

你真正能看到的最小非零模式不是 $k=0$，而是

$$k_{\min}.$$

因此只要所有允许的非零模式都稳定，模拟中就可以看到一个稳定的涨落抑制态。

因为上面的两个稳定条件随 $k^2$ 增大而变得更容易满足，所以最危险的是

$$k=k_{\min}.$$

于是有限系统中的稳定条件是：

$$\boxed{ \mu< \mu_c(L) }$$

其中

$$\boxed{ \mu_c(L) = \min \left[ (D+D_c)k_{\min}^2, D_ck_{\min}^2+\frac{\alpha\rho_0}{D} \right]. }$$

也就是

$$\boxed{ \mu_c(L) = \min \left[ (D+D_c)\left(\frac{2\pi}{L}\right)^2, D_c\left(\frac{2\pi}{L}\right)^2+\frac{\alpha\rho_0}{D} \right]. }$$

如果 $\alpha>0$ 且 $\alpha\rho_0/D$ 不太小，通常第一项更严格：

$$\boxed{ \mu_c(L) \approx (D+D_c)\left(\frac{2\pi}{L}\right)^2. }$$

这说明：

$$\boxed{ \mu>0 \text{ 但很小，完全可能在有限系统中保持稳定。} }$$

但当你增大系统尺寸 $L$ 时，

$$k_{\min}\sim \frac1L,$$

所以

$$\mu_c(L)\sim \frac1{L^2}.$$

因此固定一个正的 $\mu$，只要 $L$ 足够大，最终一定会进入不稳定区。

---

# 6. 这如何解释你的数值结果？

你说：

> 对于正的小 $\mu$ 还是可以形成类似于 hyperuniform 的涨落抑制现象。

这很合理，可能有三种原因。

---

## 原因一：你的系统有限，最小波数还没有进入失稳区

如果

$$\mu<(D+D_c)k_{\min}^2,$$

那么所有非零模式的 trace 条件都满足：

$$\mathrm{Tr}A(k)<0.$$

这时虽然严格 $k\to0$ 下 $\mu>0$ 不稳定，但你的模拟根本没有那些更小的 $k$。

因此在可观测波数范围内，系统可以表现得稳定，并且有涨落抑制。

---

## 原因二：$k=0$ 化学场模式不影响粒子运动

对于 $k=0$：

$$A(0)= \begin{pmatrix} 0 & 0\\ -1 & \mu \end{pmatrix}.$$

本征值是

$$s_1=0, \qquad s_2=\mu.$$

所以 $\mu>0$ 时，$c$ 的空间均值模式会增长。

但是粒子运动由

$$\nabla c$$

决定。空间均匀的 $c$ 模式满足

$$\nabla c_{k=0}=0.$$

所以它不会直接驱动粒子聚集或扩散。

如果你的数值代码中对 $c$ 做了减均值、重标定，或者只通过梯度使用 $c$，那么 $k=0$ 的指数增长可能不会影响密度结构因子。

真正重要的是 $k>0$ 的模式稳定性。

---

## 原因三：你看到的是有限波数窗口中的 apparent hyperuniformity

定义

$$k_\mu=\sqrt{\frac{\mu}{D_c}}.$$

当

$$k\gg k_\mu,$$

有

$$D_ck^2\gg\mu.$$

此时

$$D_ck^2-\mu\approx D_ck^2,$$

所以系统看起来接近

$$\mu=0$$

的无屏蔽化学场情形。

如果在你的可观测波数范围内都有

$$k\gg k_\mu,$$

那么 $\mu>0$ 的影响几乎看不出来，系统可以表现出类似 hyperuniform 的下降趋势。

但如果继续往更小 $k$ 看，理论上会遇到 crossover 或失稳。

---

# 7. 对 $S(k)$ 的改进：不要直接用快化学场近似

之前我们用快场近似得到

$$S(k) \simeq \rho_0 \frac{D} { D+\dfrac{\alpha\rho_0}{D_ck^2-\mu} }.$$

这个式子在

$$D_ck^2-\mu>0$$

且化学场相对较快时有用。

但当 $\mu>0$ 且接近稳定边界时，

$$D_ck^2-\mu$$

可能很小，甚至为负。此时不能简单消去 $\delta c$。更稳妥的是使用完整二变量 Ornstein-Uhlenbeck 系统。

定义

$$a=Dk^2,$$

$$q=D_ck^2-\mu,$$

$$b=\alpha\rho_0k^2.$$

线性系统为

$$\partial_t \begin{pmatrix} \delta\rho_k\\ \delta c_k \end{pmatrix} = \begin{pmatrix} -a & b\\ -1 & -q \end{pmatrix} \begin{pmatrix} \delta\rho_k\\ \delta c_k \end{pmatrix} + \begin{pmatrix} \eta_k\\ 0 \end{pmatrix}.$$

噪声满足

$$\langle\eta_k(t)\eta_{-k}(t')\rangle = 2D\rho_0k^2\delta(t-t').$$

在稳定条件满足时，稳态结构因子可由 Lyapunov 方程得到：

$$\boxed{ S(k) = D\rho_0k^2 \frac{ b+q(a+q) } { (a+q)(aq+b) }. }$$

代回 $a,q,b$：

$$\boxed{ S(k) = D\rho_0k^2 \frac{ \alpha\rho_0k^2 + (D_ck^2-\mu)\left[(D+D_c)k^2-\mu\right] } { \left[(D+D_c)k^2-\mu\right] \left[ Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2 \right] }. }$$

这个表达式比快场近似更适合分析 $\mu>0$ 的小正值情况。

它要求分母中的两个稳定因子为正：

$$(a+q)>0,$$

$$(aq+b)>0.$$

也就是：

$$(D+D_c)k^2-\mu>0,$$

$$Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2>0.$$

这正对应前面的稳定条件。

---

# 8. 小正 $\mu$ 下的物理图像

对于某个可观测波数 $k$，如果

$$0<\mu\ll D_ck^2,$$

那么

$$q=D_ck^2-\mu>0.$$

化学场对该模式仍然是有效衰减的。

此时密度涨落会受到趋化反馈压制，机制仍然是：

$$\delta\rho_k \rightarrow \delta c_k \rightarrow \text{chemotactic flux} \rightarrow \text{opposes density fluctuation}.$$

所以小正 $\mu$ 不会立刻破坏涨落抑制。

真正的问题出现在

$$k\lesssim k_\mu=\sqrt{\frac{\mu}{D_c}}.$$

此时

$$D_ck^2-\mu<0,$$

化学场自身对该波数有增长倾向。若同时

$$\mu>(D+D_c)k^2,$$

则线性系统失稳。

因此对于正小 $\mu$，系统会有一个很重要的长度尺度：

$$\ell_\mu\sim \sqrt{\frac{D_c}{\mu}}.$$

只有当系统尺寸

$$L\ll \ell_\mu$$

时，你才比较容易在整个系统尺度上看到稳定的涨落抑制。

反过来，如果

$$L\gtrsim \ell_\mu,$$

长波模式会开始变得危险。

更准确地，有限系统稳定要求

$$\mu<(D+D_c)\left(\frac{2\pi}{L}\right)^2.$$

等价于

$$L<2\pi\sqrt{\frac{D+D_c}{\mu}}.$$

---

# 9. 你数值中看到 “hyperuniform-like” 的可能解释

如果你的参数满足：

$$0<\mu<\mu_c(L),$$

那么有限系统中所有可见非零模式都稳定。

如果进一步满足，在测量的波数窗口内：

$$D_ck^2\gg\mu,$$

那么系统会近似表现得像 $\mu=0$ 或弱屏蔽化学场。

在这种情况下，$S(k)$ 可能在有限窗口中明显下降，看起来像

$$S(k)\sim k^\alpha.$$

但严格地说，这种情况下要非常小心区分两种可能：

$$S(k)\sim k^\alpha$$

和

$$S(k)=S_0+A k^2+\cdots .$$

这和我们之前分析 antialignment 那篇文章时得到的警示类似：有限 $k$ 范围内的幂律外观，可能不是严格的 $k\to0$ 极限。

---

# 10. 对你后续数值验证的建议

你可以检查三个量。

第一，计算

$$k_{\min}=\frac{2\pi}{L}$$

并检查：

$$\mu<(D+D_c)k_{\min}^2.$$

如果成立，说明你的正小 $\mu$ 在当前有限系统里没有进入长波不稳定区。

第二，检查：

$$k_\mu=\sqrt{\frac{\mu}{D_c}}.$$

如果你测量 $S(k)$ 的波数大多满足

$$k\gg k_\mu,$$

那么你看到的涨落抑制很可能是“无屏蔽近似区间”的行为。

第三，做有限尺寸扫描。固定 $\mu>0$，逐渐增大 $L$。理论预测是：

$$\mu_c(L)\sim L^{-2}.$$

所以足够大 $L$ 后，应当出现长波失稳或至少最小波数处的涨落增强。

如果你看到随着 $L$ 增大，最小波数处的 $S(k_{\min})$ 仍持续下降并趋近 0，那说明模型中可能还有当前线性理论没有包含的稳定化机制。

---

# 11. 总结

对于 $\mu>0$，严格热力学极限下：

$$\boxed{ \text{任意正的 } \mu \text{ 都会使足够长波的化学场模式不稳定。} }$$

但对于有限系统，稳定条件是：

$$\boxed{ \mu< \mu_c(L) = \min \left[ (D+D_c)\left(\frac{2\pi}{L}\right)^2, D_c\left(\frac{2\pi}{L}\right)^2+\frac{\alpha\rho_0}{D} \right]. }$$

在这个范围内，正小 $\mu$ 仍然可以表现出稳定的密度涨落抑制。

更完整的结构因子估计应使用：

$$\boxed{ S(k) = D\rho_0k^2 \frac{ \alpha\rho_0k^2 + (D_ck^2-\mu)\left[(D+D_c)k^2-\mu\right] } { \left[(D+D_c)k^2-\mu\right] \left[ Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2 \right] }. }$$

这个公式只在对应 $k$ 模式稳定时有效。

所以你的数值现象可以这样解释：

$$\boxed{ \mu>0 \text{ 但很小，在有限系统和有限波数窗口中，仍可出现 apparent hyperuniform-like density suppression。} }$$

但它是否是真正的

$$S(k\to0)\to0$$

需要通过增大系统尺寸、降低可访问的最小 $k$，以及检查失稳阈值

$$\mu_c(L)\sim L^{-2}$$

来判断。

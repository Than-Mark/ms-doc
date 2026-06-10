你的理解里有一半是对的：**粒子噪声在平均场密度方程里确实会给出 $D\nabla^2\rho$**。但这只是**平均漂移项**，不是完整的涨落描述。

更完整地说：

$$D\nabla^2\rho$$

来自粒子布朗噪声的**平均效应**；而

$$-\nabla\cdot\left[\sqrt{2D\rho},\boldsymbol{\Lambda}(\mathbf r,t)\right]$$

来自同一个粒子布朗噪声的**样本间随机涨落效应**。

也就是说，二者不是互相替代关系，而是同一微观噪声在密度场层面的两个结果：

$$\text{粒子布朗噪声} \quad\Longrightarrow\quad D\nabla^2\rho + \text{守恒型密度噪声}.$$

你给的模型中，粒子运动方程本身包含

$$\dot{\mathbf r}_i = \alpha\nabla c(\mathbf r_i,t) + \sqrt{2D}\boldsymbol{\xi}_i(t),$$

其中 $D$ 是粒子扩散系数，$\xi$ 是零均值单位方差高斯白噪声；连续模型则写成平均场形式

$$\partial_t\rho = -\alpha\nabla\cdot(\rho\nabla c) + D\nabla^2\rho .$$

文档也说明这个 continuum model 是在 mean-field approximation 下得到的，并假设密度涨落可忽略；所以它保留了平均扩散项，但没有保留涨落噪声项。

---

## 1. 从粒子密度定义开始

定义 microscopic density：

$$\rho(\mathbf r,t) = \sum_{i=1}^N \delta(\mathbf r-\mathbf r_i(t)).$$

对时间求导：

$$\partial_t\rho(\mathbf r,t) = \sum_i \partial_t \delta(\mathbf r-\mathbf r_i(t)).$$

用链式法则：

$$\partial_t \delta(\mathbf r-\mathbf r_i(t)) = -\dot{\mathbf r}_i(t)\cdot \nabla \delta(\mathbf r-\mathbf r_i(t)).$$

所以

$$\partial_t\rho = -\nabla\cdot \sum_i \dot{\mathbf r}_i \delta(\mathbf r-\mathbf r_i).$$

这其实就是一个连续性方程：

$$\partial_t\rho = -\nabla\cdot\mathbf J.$$

粒子数守恒决定了噪声不能直接写成

$$+\eta(\mathbf r,t),$$

而必须写成某个随机流的散度：

$$-\nabla\cdot\mathbf J_{\rm noise}.$$

这就是“守恒噪声”的来源。

---

## 2. 代入粒子 Langevin 方程

粒子方程是

$$d\mathbf r_i = \alpha\nabla c(\mathbf r_i,t)dt + \sqrt{2D}d\mathbf W_i(t).$$

把它代入密度连续性方程，可以分成两部分。

确定性趋化漂移给出：

$$-\nabla\cdot \left[ \alpha \sum_i \nabla c(\mathbf r_i,t) \delta(\mathbf r-\mathbf r_i) \right].$$

在局域场近似下，

$$\sum_i \nabla c(\mathbf r_i,t) \delta(\mathbf r-\mathbf r_i) \approx \rho(\mathbf r,t)\nabla c(\mathbf r,t),$$

于是得到

$$-\alpha\nabla\cdot(\rho\nabla c).$$

这就是你 PDF 里的趋化项。

---

## 3. 为什么会出现 $D\nabla^2\rho$？

这里要用 Itô 公式。因为粒子运动含有布朗项，所以对

$$\delta(\mathbf r-\mathbf r_i(t))$$

求微分时，不只有一阶链式法则项，还会有二阶 Itô 修正项。

对单个粒子：

$$d\delta(\mathbf r-\mathbf r_i) = -\nabla\delta(\mathbf r-\mathbf r_i)\cdot d\mathbf r_i + D\nabla^2\delta(\mathbf r-\mathbf r_i)dt.$$

把所有粒子加起来：

$$\sum_i D\nabla^2\delta(\mathbf r-\mathbf r_i) = D\nabla^2 \sum_i\delta(\mathbf r-\mathbf r_i) = D\nabla^2\rho.$$

所以

$$D\nabla^2\rho$$

确实来自粒子布朗噪声，但它是 Itô 二阶项，也就是布朗运动对平均密度的扩散效应。

---

## 4. 那随机噪声项从哪里来？

同一个 Itô 展开里还有一阶随机项：

$$-\sum_i \nabla\delta(\mathbf r-\mathbf r_i) \cdot \sqrt{2D},d\mathbf W_i(t).$$

它可以写成散度形式：

$$-\nabla\cdot \left[ \sqrt{2D} \sum_i \delta(\mathbf r-\mathbf r_i) \boldsymbol{\xi}_i(t) \right].$$

定义 microscopic random current：

$$\mathbf J_{\rm noise}(\mathbf r,t) = \sqrt{2D} \sum_i \delta(\mathbf r-\mathbf r_i) \boldsymbol{\xi}_i(t).$$

于是密度方程变成：

$$\partial_t\rho = -\alpha\nabla\cdot(\rho\nabla c) + D\nabla^2\rho - \nabla\cdot\mathbf J_{\rm noise}.$$

这个式子仍然是精确的 microscopic density equation，只是 $\mathbf J_{\rm noise}$ 还写成粒子求和形式。

---

## 5. 为什么可以把它写成 $\sqrt{2D\rho}\boldsymbol{\Lambda}$？

我们计算这个随机流的相关函数。

设

$$\langle \xi_{i,a}(t)\xi_{j,b}(t')\rangle = \delta_{ij}\delta_{ab}\delta(t-t').$$

那么

$$\langle J_{{\rm noise},a}(\mathbf r,t) J_{{\rm noise},b}(\mathbf r',t')\rangle = 2D \delta_{ab} \delta(t-t') \sum_i \delta(\mathbf r-\mathbf r_i) \delta(\mathbf r'-\mathbf r_i).$$

利用

$$\delta(\mathbf r-\mathbf r_i)\delta(\mathbf r'-\mathbf r_i) = \delta(\mathbf r-\mathbf r') \delta(\mathbf r-\mathbf r_i),$$

得到

$$\langle J_{{\rm noise},a}(\mathbf r,t) J_{{\rm noise},b}(\mathbf r',t')\rangle = 2D\rho(\mathbf r,t) \delta_{ab} \delta(\mathbf r-\mathbf r') \delta(t-t').$$

这说明随机流可以等价写成：

$$\mathbf J_{\rm noise} = \sqrt{2D\rho(\mathbf r,t)} \boldsymbol{\Lambda}(\mathbf r,t),$$

其中

$$\langle \Lambda_a(\mathbf r,t)\Lambda_b(\mathbf r',t')\rangle = \delta_{ab} \delta(\mathbf r-\mathbf r') \delta(t-t').$$

于是得到 Dean 型密度方程：

$$\boxed{ \partial_t\rho = -\alpha\nabla\cdot(\rho\nabla c) + D\nabla^2\rho - \nabla\cdot \left[ \sqrt{2D\rho}\boldsymbol{\Lambda} \right]. }$$

所以答案是：**可以从粒子模型推导出来，而且它和 $D\nabla^2\rho$ 来自同一个微观布朗噪声，只是一个是平均扩散，一个是涨落流。**

---

## 6. 为什么平均场方程里只剩 $D\nabla^2\rho$？

因为平均场方程通常对噪声取平均。

由于

$$\langle \boldsymbol{\Lambda}\rangle=0,$$

所以

$$\left\langle -\nabla\cdot \left[ \sqrt{2D\rho}\boldsymbol{\Lambda} \right] \right\rangle \approx0.$$

因此平均密度满足：

$$\partial_t\langle\rho\rangle = -\alpha\nabla\cdot(\langle\rho\rangle\nabla c) + D\nabla^2\langle\rho\rangle.$$

这就是你 PDF 中的 continuum model。文档也明确说该连续模型通过 mean-field approximation 得到，并假设 fluctuations around the mean density are negligible。

但如果你要算

$$S(k)=\langle \delta\rho_k\delta\rho_{-k}\rangle,$$

就不能把这个零均值噪声丢掉。因为结构因子本身就是二阶涨落量。正如前面学习笔记总结的那样，研究 hyperuniform 不能只看平均场，必须研究噪声怎样进入粗粒化密度动力学。

---

## 7. 最容易混淆的点

你的直觉：

$$\text{粒子噪声} \Rightarrow D\nabla^2\rho$$

是对的，但更完整地应写成：

$$
\text{粒子噪声}\Rightarrow \begin{cases}
	D\nabla ^2\rho ,&		\text{平均扩散效应},\\
	-\nabla \cdot [\sqrt{2D\rho}\boldsymbol{\Lambda }],&		\text{密度涨落效应}.\\
\end{cases}
$$

如果只研究平均密度演化，可以只保留第一项。

如果研究结构因子、数目涨落、hyperuniformity，就必须保留第二项。

---

## 8. 在线性化时它变成什么？

围绕均匀态

$$\rho=\rho_0+\delta\rho$$

线性化时，

$$\sqrt{2D\rho} \approx \sqrt{2D\rho_0}.$$

因此噪声项变成加性守恒噪声：

$$-\nabla\cdot \left[ \sqrt{2D\rho_0} \boldsymbol{\Lambda} \right].$$

Fourier 空间中：

$$\eta_{\mathbf k} = -i\mathbf k\cdot \sqrt{2D\rho_0} \boldsymbol{\Lambda}_{\mathbf k}.$$

所以

$$\langle \eta_{\mathbf k}(t)\eta_{-\mathbf k}(t')\rangle = 2D\rho_0 k^2\delta(t-t').$$

这就是为什么我前面说：密度噪声在小 $k$ 下带 $k^2$ 强度。这个 $k^2$ 不是随便加的，而是粒子数守恒和布朗噪声共同决定的。

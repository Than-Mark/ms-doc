这个表达式是从二维线性 Ornstein-Uhlenbeck 系统的稳态协方差方程，也就是 Lyapunov 方程，直接解出来的。

我们从线性系统开始：

$$\frac{d}{dt} \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix} = \begin{pmatrix} -Dk^2 & \alpha\rho_0k^2\\ -1 & -(D_ck^2-\mu) \end{pmatrix} \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix} + \begin{pmatrix} \eta_{\mathbf k}\\ 0 \end{pmatrix}.$$

为了简化记号，定义

$$a=Dk^2,$$

$$b=\alpha\rho_0k^2,$$

$$q=D_ck^2-\mu.$$

于是矩阵写成

$$A= \begin{pmatrix} -a & b\\ -1 & -q \end{pmatrix}.$$

噪声满足

$$\langle \eta_{\mathbf k}(t)\eta_{-\mathbf k}(t')\rangle = 2D\rho_0k^2\delta(t-t').$$

记

$$\Gamma=D\rho_0k^2.$$

则噪声协方差矩阵是

$$Q= \begin{pmatrix} 2\Gamma & 0\\ 0 & 0 \end{pmatrix}.$$

---

## 1. 写出协方差矩阵

我们要算的是

$$S(k)= \langle \delta\rho_{\mathbf k}\delta\rho_{-\mathbf k}\rangle.$$

定义协方差矩阵

$$C= \left\langle \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix} \begin{pmatrix} \delta\rho_{-\mathbf k} & \delta c_{-\mathbf k} \end{pmatrix} \right\rangle.$$

写成分量形式：

$$C= \begin{pmatrix} x & y\\ y & z \end{pmatrix},$$

其中

$$x=\langle \delta\rho_{\mathbf k}\delta\rho_{-\mathbf k}\rangle=S(k),$$

$$y=\langle \delta\rho_{\mathbf k}\delta c_{-\mathbf k}\rangle,$$

$$z=\langle \delta c_{\mathbf k}\delta c_{-\mathbf k}\rangle.$$

我们要求的就是 $x$。

---

## 2. 稳态 Lyapunov 方程

对于线性随机系统

$$\dot{\mathbf u}=A\mathbf u+\boldsymbol{\eta},$$

如果系统稳定，稳态协方差满足

$$AC+CA^T+Q=0.$$

这里

$$\mathbf u= \begin{pmatrix} \delta\rho_{\mathbf k}\\ \delta c_{\mathbf k} \end{pmatrix}.$$

所以我们要解

$$\begin{pmatrix} -a & b\\ -1 & -q \end{pmatrix} \begin{pmatrix} x & y\\ y & z \end{pmatrix} + \begin{pmatrix} x & y\\ y & z \end{pmatrix} \begin{pmatrix} -a & -1\\ b & -q \end{pmatrix} + \begin{pmatrix} 2\Gamma & 0\\ 0 & 0 \end{pmatrix} =0.$$

---

## 3. 展开矩阵乘法

先算

$$AC= \begin{pmatrix} -a & b\\ -1 & -q \end{pmatrix} \begin{pmatrix} x & y\\ y & z \end{pmatrix} = \begin{pmatrix} -ax+by & -ay+bz\\ -x-qy & -y-qz \end{pmatrix}.$$

再算

$$CA^T= \begin{pmatrix} x & y\\ y & z \end{pmatrix} \begin{pmatrix} -a & -1\\ b & -q \end{pmatrix} = \begin{pmatrix} -ax+by & -x-qy\\ -ay+bz & -y-qz \end{pmatrix}.$$

所以

$$AC+CA^T+Q=0$$

给出三个独立方程。

---

## 4. 分量方程

### 第一行第一列

$$(-ax+by)+(-ax+by)+2\Gamma=0.$$

即

$$-2ax+2by+2\Gamma=0.$$

两边除以 $2$：

$$-ax+by+\Gamma=0.$$

所以

$$\boxed{ ax-by=\Gamma. }$$

---

### 第一行第二列

$$(-ay+bz)+(-x-qy)=0.$$

整理：

$$-x-(a+q)y+bz=0.$$

所以

$$\boxed{ x+(a+q)y-bz=0. }$$

---

### 第二行第二列

$$(-y-qz)+(-y-qz)=0.$$

即

$$-2y-2qz=0.$$

所以

$$\boxed{ y=-qz. }$$

---

## 5. 解这组三元线性方程

我们有：

$$ax-by=\Gamma,$$

$$x+(a+q)y-bz=0,$$

$$y=-qz.$$

由第三个式子：

$$z=-\frac{y}{q}.$$

代入第二个式子：

$$x+(a+q)y-b\left(-\frac{y}{q}\right)=0.$$

也就是

$$x+(a+q)y+\frac{b}{q}y=0.$$

所以

$$x+ \left( a+q+\frac{b}{q} \right)y=0.$$

得到

$$y= -\frac{x}{a+q+b/q}.$$

把分母通分：

$$a+q+\frac{b}{q} = \frac{q(a+q)+b}{q}.$$

所以

$$\boxed{ y= -\frac{xq}{q(a+q)+b}. }$$

现在代入第一个方程：

$$ax-by=\Gamma.$$

代入 $y$：

$$ax - b\left( -\frac{xq}{q(a+q)+b} \right) = \Gamma.$$

所以

$$x \left[ a+ \frac{bq}{q(a+q)+b} \right] = \Gamma.$$

因此

$$x = \Gamma \frac{q(a+q)+b} { a[q(a+q)+b]+bq }.$$

整理分母：

$$a[q(a+q)+b]+bq = aq(a+q)+ab+bq.$$

注意

$$ab+bq=b(a+q).$$

所以

$$aq(a+q)+b(a+q) = (a+q)(aq+b).$$

因此

$$\boxed{ x = \Gamma \frac{ q(a+q)+b } { (a+q)(aq+b) }. }$$

因为

$$x=S(k),$$

所以

$$\boxed{ S(k) = \Gamma \frac{ q(a+q)+b } { (a+q)(aq+b) }. }$$

---

## 6. 代回 $a,b,q,\Gamma$

现在代回：

$$a=Dk^2,$$

$$b=\alpha\rho_0k^2,$$

$$q=D_ck^2-\mu,$$

$$\Gamma=D\rho_0k^2.$$

于是

$$S(k) = D\rho_0k^2 \frac{ (D_ck^2-\mu)\left[Dk^2+(D_ck^2-\mu)\right] + \alpha\rho_0k^2 } { \left[Dk^2+(D_ck^2-\mu)\right] \left[Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2\right] }.$$

也可以写成

$$\boxed{ S(k) = D\rho_0k^2 \frac{ \alpha\rho_0k^2 + (D_ck^2-\mu)\left[(D+D_c)k^2-\mu\right] } { \left[(D+D_c)k^2-\mu\right] \left[ Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2 \right] }. }$$

如果把分母完全展开，就是：

$$\left[(D+D_c)k^2-\mu\right] \left[ Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2 \right].$$

展开后：

$$D^2k^4(D_ck^2-\mu) + D\alpha\rho_0k^4 + Dk^2(D_ck^2-\mu)^2 + \alpha\rho_0k^2(D_ck^2-\mu).$$

分子展开为：

$$Dk^2(D_ck^2-\mu) + (D_ck^2-\mu)^2 + \alpha\rho_0k^2.$$

所以得到之前写的形式：

$$\boxed{ S(k) = D\rho_0k^2 \frac{ Dk^2(D_ck^2-\mu) + \alpha\rho_0k^2 + (D_ck^2-\mu)^2 } { D^2k^4(D_ck^2-\mu) + D\alpha\rho_0k^4 + Dk^2(D_ck^2-\mu)^2 + \alpha\rho_0k^2(D_ck^2-\mu) }. }$$

---

## 7. 这个公式什么时候有效？

这个公式要求线性系统有稳定稳态。也就是矩阵 $A(k)$ 的两个本征值实部都为负。

对于

$$A= \begin{pmatrix} -a & b\\ -1 & -q \end{pmatrix},$$

稳定条件是：

$$a+q>0,$$

$$aq+b>0.$$

代回原变量：

$$\boxed{ (D+D_c)k^2-\mu>0, }$$

$$\boxed{ Dk^2(D_ck^2-\mu)+\alpha\rho_0k^2>0. }$$

也就是

$$\boxed{ \mu<(D+D_c)k^2, }$$

以及

$$\boxed{ \alpha\rho_0+D(D_ck^2-\mu)>0. }$$

如果这些条件不满足，Lyapunov 方程没有正定稳态解，结构因子作为稳态量就不成立。

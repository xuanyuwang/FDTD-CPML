
# FDTD-CPML 2D FDTD, TM, Hx, Hy & Ez.
FDTD with CPML (Convolutional Perfectly Matched Layer)

## Pre: Coefficients and The Set of PML

There are two types of form of coordinate stretched factor. The first is:
$$s_w=\kappa_w
+
\frac{\sigma_w}{j\omega\epsilon_0},$$
and the other is:
$$s_w=\kappa_w
+
\frac{\sigma_w}{a_w + j\omega\epsilon_0}.$$
Here, we adopt the second one.

$$\alpha = \frac{\sigma_w}{\epsilon_0 \kappa_w}$$

$$c_{w}=\frac{1}{\kappa_w}\left\{
    exp\left[
        -\frac{\sigma_w\Delta t}{\epsilon_0 \kappa_w}
        \right]-1
    \right\}$$

$$\sigma_{max}=\frac{m+1}{\sqrt{\epsilon_r}150\pi\delta}$$

$$\sigma_z(z) = \frac{\sigma_{max}|z-z_0|^m}{d^m}$$

$k_{max}=5 ~--~ 11$, There is 8. $d$ is the width of the PML layer. $m$ is a integer, we set $m=4$.

$$\kappa_z(z) = 1 + (\kappa_{max} - 1)\frac{|z-z_0|^m}{d^m}$$

$$CA(m) = \frac{
    1-\frac{\sigma(m)\Delta t}{2\epsilon(m)}
}{
1+\frac{\sigma(m)\Delta t}{2\epsilon(m)}    
}$$

$$CB(m)=\frac{
    \frac{\Delta t}{\epsilon(m)}
}{
 1+\frac{\sigma(m)\Delta t}{2\epsilon(m)}   
}$$

$$CP(m) = \frac{
    1-\frac{\sigma_m(m)\Delta t}{2\mu(m)}
}{
  1+\frac{\sigma_m(m)\Delta t}{2\mu(m)}  
}$$

$$CQ(m) = \frac{
    \frac{\Delta t}{\mu(m)}
}{
    1+\frac{\sigma_m(m)\Delta t}{2\mu(m)}
}$$

if $\sigma=\sigma_m=0$, 
so,$CA(m)=1, CB(m)=\frac{\Delta t}{\epsilon}, CP(m)=1, CQ(m)=\frac{\Delta t}{\mu}$

$$\kappa = \sqrt{\omega^2 \mu\epsilon}$$

## 3D equations, continual form

$$\frac{\partial D}{\partial t} 
=
\hat{x}\left(
\frac{1}{\kappa_y}\frac{\partial H_z}{\partial y}-\frac{1}{\kappa_z}\frac{\partial H_y}{\partial z}+\zeta_y(t)*\frac{\partial H_z}{\partial y}-\zeta_z(t)*\frac{\partial H_y}{\partial z}
\right)\\
+
\hat{y}\left(
\frac{1}{\kappa_z}\frac{\partial H_x}{\partial z}-\frac{1}{\kappa_x}\frac{\partial H_z}{\partial x}+\zeta_z(t)*\frac{\partial H_x}{\partial z}-\zeta_x(t)*\frac{\partial H_z}{\partial x}
\right)\\
+
\hat{z}\left(
\frac{1}{\kappa_x}\frac{\partial H_y}{\partial x}-\frac{1}{\kappa_y}\frac{\partial H_x}{\partial y}+\zeta_x(t)*\frac{\partial H_y}{\partial x}-\zeta_y(t)*\frac{\partial H_x}{\partial y}
\right)$$

<br/>

$$-\frac{\partial B}{\partial t} 
=
\hat{x}\left(
\frac{1}{\kappa_y}\frac{\partial E_z}{\partial y}-\frac{1}{\kappa_z}\frac{\partial E_y}{\partial z}+\zeta_y(t)*\frac{\partial E_z}{\partial y}-\zeta_z(t)*\frac{\partial E_y}{\partial z}
\right)\\
+
\hat{y}\left(
\frac{1}{\kappa_z}\frac{\partial E_x}{\partial z}-\frac{1}{\kappa_x}\frac{\partial E_z}{\partial x}+\zeta_z(t)*\frac{\partial E_x}{\partial z}-\zeta_x(t)*\frac{\partial E_z}{\partial x}
\right)\\
+
\hat{z}\left(
\frac{1}{\kappa_x}\frac{\partial E_y}{\partial x}-\frac{1}{\kappa_y}\frac{\partial E_x}{\partial y}+\zeta_x(t)*\frac{\partial E_y}{\partial x}-\zeta_y(t)*\frac{\partial E_x}{\partial y}
\right)$$

<br/>

$$\psi_{H_{wv}}(n) = c_w \frac{\partial E_v(n)}{\partial w} + \exp(-\alpha\Delta t)\psi_{H_{wv}}(n-1)$$

<br/>

$$\psi_{E_{wv}}(n) = c_w \frac{\partial H_v(n)}{\partial w} + \exp(-\alpha\Delta t)\psi_{E_{wv}}(n-1)$$

## 2D equations, continual form, TM($E_z, H_x, H_y$)

$$\frac{\partial}{\partial t}(\epsilon E_z)=
\frac{1}{\kappa_x}\frac{\partial H_y}{\partial x}-\frac{1}{\kappa_y}\frac{\partial H_x}{\partial y}+\zeta_x(t)*\frac{\partial H_y}{\partial x}-\zeta_y(t)*\frac{\partial H_x}{\partial y}
=
\frac{1}{\kappa_x}\frac{\partial H_y}{\partial x}-\frac{1}{\kappa_y}\frac{\partial H_x}{\partial y}+\psi_{E_{xy}}-\psi_{E_{yx}}$$

<br/>

$$-\frac{\partial}{\partial t}(\mu H_x)=
\frac{1}{\kappa_y}\frac{\partial E_z}{\partial y}-\frac{1}{\kappa_z}\frac{\partial E_y}{\partial z}+\zeta_y(t)*\frac{\partial E_z}{\partial y}-\zeta_z(t)*\frac{\partial E_y}{\partial z}
=
\frac{1}{\kappa_y}\frac{\partial E_z}{\partial y}-\frac{1}{\kappa_z}\frac{\partial E_y}{\partial z}+\psi_{H_{yz}}-\psi_{H_{zy}}$$

<br/>

$$-\frac{\partial}{\partial t}(\mu H_y)=
\frac{1}{\kappa_z}\frac{\partial E_x}{\partial z}-\frac{1}{\kappa_x}\frac{\partial E_z}{\partial x}+\zeta_z(t)*\frac{\partial E_x}{\partial z}-\zeta_x(t)*\frac{\partial E_z}{\partial x}
=
\frac{1}{\kappa_z}\frac{\partial E_x}{\partial z}-\frac{1}{\kappa_x}\frac{\partial E_z}{\partial x}+\psi_{H_{zx}}-\psi_{H_{xz}}$$

<br/>

## 2D equations, discrete form, TM($E_z, H_x, H_y$)

### 1. Normal FDTD

$$E_z^{n+1}(i,j)
=
E_z^n(i,j)
+\frac{\Delta t}{\epsilon}\left[
\frac{H_y^{n+1/2}(i+\frac{1}{2},j)-H_y^{n+1/2}(i-\frac{1}{2},j)}{\kappa_x \Delta x}-\frac{H_x^{n+1/2}(i,j+\frac{1}{2})-H_x^{n+1/2}(i,j-\frac{1}{2})}{\kappa_y \Delta y}
\right]
+\frac{\Delta t}{\epsilon}\left[
\psi_{E_{xy}}^{n+\frac{1}{2}}\left(i,j\right)-\psi_{E_{yx}}^{n+\frac{1}{2}}\left(i,j\right)
\right]$$

<br/>

$$H_x^{n+1/2}\left( i, j+\frac{1}{2} \right)
=
H_x^{n-1/2}\left( i, j+\frac{1}{2} \right)
-\frac{\Delta t}{\mu}\frac{E_z^n(i,j+1)-E_z^n(i,j)}{\kappa_y \Delta y}
+
\frac{\Delta t}{\mu}\left[
\psi^n_{H_{zy}}\left( i, j+\frac{1}{2} \right)-\psi^n_{H_{yz}}\left( i, j+\frac{1}{2} \right)
\right]$$

<br/>

$$H_y^{n+1/2}\left( i+\frac{1}{2},j \right)
=
H_y^{n-1/2}\left( i+\frac{1}{2},j \right)
+
\frac{\Delta t}{\mu}
\frac{E_z^n(i+1,j)-E_z^n(i,j)}{\kappa_x \Delta x}
+
\frac{\Delta t}{\mu}
\left[
\psi^n_{H_{xz}}\left( i+\frac{1}{2},j \right)-\psi^n_{H_{zx}}\left( i+\frac{1}{2},j \right)
\right]$$

<br>

### 2. Convolutions

$$\psi_{E_{xy}}^{n+1/2}(i,j) = c_x \frac{H_y^{n+1/2}(i+1/2,j)-H_y^{n+1/2}(i-1/2,j)}{\Delta x} + \exp(-\alpha\Delta t)\psi_{E_{xy}}^{n-1/2}(i,j)$$

$$\psi_{E_{yx}}^{n+1/2}(i,j) = c_y \frac{H_x^{n+1/2}(i,j+1/2)-H_x^{n+1/2}(i,j-1/2)}{\Delta y} + \exp(-\alpha\Delta t)\psi_{E_{yx}}^{n-1/2}(i,j)$$

$$\psi_{H_{wv}}(n) = c_w \frac{\partial E_v(n)}{\partial w} + \exp(-\alpha\Delta t)\psi_{H_{wv}}(n-1)$$

$$\psi_{H_{zy}}^{n}\left( i,j+\frac{1}{2} \right)
= 
0$$

$$\psi_{H_{yz}}^n \left( i,j+\frac{1}{2}\right)
=
c_y\frac{E_z^n(i,j+1)-E_z^n(i,j)}{\Delta y}
+
\exp(-\alpha\Delta t)\psi_{H_{yz}}^{n-1} \left( i,j+\frac{1}{2}\right)$$

$$\psi_{H_{xz}}^n \left( i+\frac{1}{2},j\right)
=
c_x\frac{E_z^n(i+1,j)-E_z^n(i,j)}{\Delta x}
+
\exp(-\alpha\Delta t)\psi_{H_{xz}}^{n-1} \left( i,j+\frac{1}{2}\right)$$

$$\psi_{H_{zx}}^{n}\left( i+\frac{1}{2},j \right)
= 
0$$

<br>

#### Computation steps

1. $E\rightarrow H_{zx} ~\& ~ H_{xz}$
2. $E, H_{zx}, H_{xz} \rightarrow H$
3. $H\rightarrow E_{yz}, E_{zy}$
4. $H, E_{zy}, E{yz} \rightarrow E$

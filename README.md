# FDTD-CPML
FDTD with CPML (Coordinate stretched Perfectly Matched Layer)

## 1D FDTD, Ex & Hy.

### Equations

#### Pre

$\psi_{E_{wv}}(n) = 
c_{w} \frac{\partial H_{v}(n)}{\partial w}+
\exp(-\alpha \Delta t)\psi_{E_{wv}}(n-1)$

$\psi_{H_{wv}}(n) = c_w \frac{\partial E_v(n)}{\partial w} + \exp(-\alpha\Delta t)\psi_{H_{wv}}(n-1)$

$\alpha = \frac{\sigma_w}{\epsilon_0 \kappa_w} + \frac{a_w}{\epsilon_0}$

$c_{w}=\frac{\sigma_w}{\sigma_w\kappa_w+\kappa_w^2 a_w}\left\{
    exp\left[
        -\left(
            \frac{\sigma_w}{\epsilon_0 \kappa_w}+\frac{a_w}{\epsilon_0}
            \right)\Delta t
        \right]-1
    \right\}$

$CA(m) = \frac{
    1-\frac{\sigma(m)\Delta t}{2\epsilon(m)}
}{
1+\frac{\sigma(m)\Delta t}{2\epsilon(m)}    
}$

$CB(m)=\frac{
    \frac{\Delta t}{\epsilon(m)}
}{
 1+\frac{\sigma(m)\Delta t}{2\epsilon(m)}   
}$

$CP(m) = \frac{
    1-\frac{\sigma_m(m)\Delta t}{2\mu(m)}
}{
  1+\frac{\sigma_m(m)\Delta t}{2\mu(m)}  
}$

$CQ(m) = \frac{
    \frac{\Delta t}{\mu(m)}
}{
    1+\frac{\sigma_m(m)\Delta t}{2\mu(m)}
}$

$\kappa = \sqrt{\omega^2 \mu\epsilon}$

#### 3D equations

$E_{x}^{n+1}\left( i + \frac{1}{2}, j, k \right) = CA(m)E_{x}^{n}\left( i + \frac{1}{2}, j, k \right)\\
+CB(m)\left[ 
    \frac{ H_z^{n+\frac{1}{2}}\left(i+\frac{1}{2},j+\frac{1}{2},k\right)-H_z^{n+\frac{1}{2}}\left(i+\frac{1}{2},j-\frac{1}{2},k\right) }{\kappa_y(m)\Delta y}\\ 
    -
    \frac{H_y^{n+\frac{1}{2}}\left(i+\frac{1}{2},j,k+\frac{1}{2}\right)-H_y^{n+\frac{1}{2}}\left(i+\frac{1}{2},j,k-\frac{1}{2}\right)}{\kappa_z(m)\Delta z}
\right]
+CB(m)\left[
\psi_{E_{yz}}^{n+\frac{1}{2}}\left(i+\frac{1}{2},j,k\right)-\psi_{E_{zy}}^{n+\frac{1}{2}}\left(i+\frac{1}{2},j,k\right)
\right]$

$H_y^{n+\frac{1}{2}}\left( i+\frac{1}{2},j,k+\frac{1}{2} \right) = CP(m)H_y^{n-\frac{1}{2}}\left( i+\frac{1}{2},j,k+\frac{1}{2} \right)\\
-CQ(m)\left[
\frac{E_x^n\left(i+\frac{1}{2}, j, k+1\right) - E_z^n\left(i+\frac{1}{2},j,k\right)}{\kappa_y(m) \Delta z}\\
-
\frac{E_z^n\left(i+1,j,k+\frac{1}{2}\right) - E_z^n\left(i,j,k+\frac{1}{2}\right)}{\kappa_z(m) \Delta x}\right]\\
-CQ(m)\left[
\psi_{H_{zx}}^n\left(i+\frac{1}{2},j,k+\frac{1}{2}\right)-\psi^n_{H_{xz}}\left(i+\frac{1}{2},j,k+\frac{1}{2}\right)
\right]$

$\psi_{E_{yz}}^{n+\frac{1}{2}}\left( i+\frac{1}{2},j,k \right)=\exp(-\alpha \Delta t)\psi_{E_{yz}}^{n-\frac{1}{2}}\left( i+\frac{1}{2},j,k \right)\\
+c_y(m)\frac{
H_z^{n+1/2}(i+1/2,j,k+1/2)-H_z^{n+1/2}(i+1/2,j,k-1/2)
}{
\Delta y
}$

#### 1D equations, TEM mode

$E_{x}^{n+1}\left( k \right) = CA(m)E_{x}^{n}\left( k \right)\\
+CB(m)\left[ -\frac{H_y^{n+\frac{1}{2}}\left(k+\frac{1}{2}\right)-H_y^{n+\frac{1}{2}}\left(k-\frac{1}{2}\right)}{\kappa_z(m)\Delta z}
\right]
+CB(m)\left[
\psi_{E_{yz}}^{n+\frac{1}{2}}\left(k\right)-\psi_{E_{zy}}^{n+\frac{1}{2}}\left(k\right)
\right]$

$H_y^{n+\frac{1}{2}}\left(k+\frac{1}{2} \right) = CP(m)H_y^{n-\frac{1}{2}}\left( k+\frac{1}{2} \right)\\
-CQ(m)\left[
\frac{E_x^n\left(k+1\right) - E_x^n\left(k\right)}{\kappa_z \Delta z}\right]
-CQ(m)\left[
\psi_{H_{zx}}^n\left(k+\frac{1}{2}\right)-\psi^n_{H_{xz}}\left(k+\frac{1}{2}\right)
\right]$

$\psi_{E_{yz}}^{n+\frac{1}{2}}\left( k \right)=\exp(-\alpha \Delta t)\psi_{E_{yz}}^{n-\frac{1}{2}}\left( k \right)
+c_y(m)\frac{
H_z^{n+1/2}(k+1/2)-H_z^{n+1/2}(k-1/2)
}{
\Delta y
}$

if $\sigma=\sigma_m=0$, $CA(m)=1$, $CB(m)=\frac{\Delta t}{\epsilon}$, $CP(m)=1$, $CQ(m)=\frac{\Delta t}{\mu}$

so,

$E_{x}^{n+1}\left( k \right) = E_{x}^{n}\left( k \right)\\
+\frac{\Delta t}{\epsilon}\left[ -\frac{H_y^{n+\frac{1}{2}}\left(k+\frac{1}{2}\right)-H_y^{n+\frac{1}{2}}\left(k-\frac{1}{2}\right)}{\kappa_z(m)\Delta z}
\right]
+\frac{\Delta t}{\epsilon}\left[
\psi_{E_{yz}}^{n+\frac{1}{2}}\left(k\right)-\psi_{E_{zy}}^{n+\frac{1}{2}}\left(k\right)
\right]$

$H_y^{n+\frac{1}{2}}\left(k+\frac{1}{2} \right) = H_y^{n-\frac{1}{2}}\left( k+\frac{1}{2} \right)\\
-\frac{\Delta t}{\mu}\left[
\frac{E_x^n\left(k+1\right) - E_x^n\left(k\right)}{\kappa_z \Delta z}\right]
-\frac{\Delta t}{\mu}\left[
\psi_{H_{zx}}^n\left(k+\frac{1}{2}\right)-\psi^n_{H_{xz}}\left(k+\frac{1}{2}\right)
\right]$
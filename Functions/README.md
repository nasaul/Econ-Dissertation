# Simulation Design

It is based on paper :

* Name    : **Adaptive GMM Shrinkage Estimation with Consistent Moment Selection**
* Author  : **Zhipeng Liao**

## Equations:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\begin{align*}&space;Y_i&space;&=&space;X_i&space;\theta_1&space;&plus;&space;W_{1,i}&space;\theta_2&space;&plus;&space;W_{2,i}&space;\theta_3&space;&plus;&space;u_i&space;\\&space;X_i&space;&=&space;Z_{1,i}&space;\pi_1&space;&plus;&space;W'_ii&space;\pi_2&space;&plus;&space;Z_{2,i}'&space;\pi_3&space;&plus;&space;vi&space;\\&space;W_i&space;&=&space;(W_{1,i},&space;W_{2,i})&space;\\&space;Z_{2,i}'&space;&=&space;(Z_{21,i},...,&space;Z_{25,i})&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\begin{align*}&space;Y_i&space;&=&space;X_i&space;\theta_1&space;&plus;&space;W_{1,i}&space;\theta_2&space;&plus;&space;W_{2,i}&space;\theta_3&space;&plus;&space;u_i&space;\\&space;X_i&space;&=&space;Z_{1,i}&space;\pi_1&space;&plus;&space;W'_ii&space;\pi_2&space;&plus;&space;Z_{2,i}'&space;\pi_3&space;&plus;&space;vi&space;\\&space;W_i&space;&=&space;(W_{1,i},&space;W_{2,i})&space;\\&space;Z_{2,i}'&space;&=&space;(Z_{21,i},...,&space;Z_{25,i})&space;\end{align*}" title="\large \begin{align*} Y_i &= X_i \theta_1 + W_{1,i} \theta_2 + W_{2,i} \theta_3 + u_i \\ X_i &= Z_{1,i} \pi_1 + W'_ii \pi_2 + Z_{2,i}' \pi_3 + vi \\ W_i &= (W_{1,i}, W_{2,i}) \\ Z_{2,i}' &= (Z_{21,i},..., Z_{25,i}) \end{align*}" /></a>

There are 20 invalids IVs F.

## Data generation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;(W_i',&space;Z_{1,i},&space;Z_{2,i}',&space;u_i,&space;v_i,&space;F^*_i)&space;\sim&space;N(\overrightarrow{0},&space;\Sigma)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;(W_i',&space;Z_{1,i},&space;Z_{2,i}',&space;u_i,&space;v_i,&space;F^*_i)&space;\sim&space;N(\overrightarrow{0},&space;\Sigma)" title="\large (W_i', Z_{1,i}, Z_{2,i}', u_i, v_i, F^*_i) \sim N(\overrightarrow{0}, \Sigma)" /></a>

with covariance:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\begin{align*}&space;&\Sigma&space;&=&&space;diag(\Sigma_z,&space;\Sigma_{u,v},&space;I_{20})&space;\\&space;&\Sigma_z(i,j)&space;&=&&space;0.2^{|i-j|}&space;\\&space;&\Sigma_{u,v}&space;&=&&space;\begin{bmatrix}&space;1&space;&&space;0.6\\&space;0.6&space;&&space;1&space;\end{bmatrix}&space;\\&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\begin{align*}&space;&\Sigma&space;&=&&space;diag(\Sigma_z,&space;\Sigma_{u,v},&space;I_{20})&space;\\&space;&\Sigma_z(i,j)&space;&=&&space;0.2^{|i-j|}&space;\\&space;&\Sigma_{u,v}&space;&=&&space;\begin{bmatrix}&space;1&space;&&space;0.6\\&space;0.6&space;&&space;1&space;\end{bmatrix}&space;\\&space;\end{align*}" title="\large \begin{align*} &\Sigma &=& diag(\Sigma_z, \Sigma_{u,v}, I_{20}) \\ &\Sigma_z(i,j) &=& 0.2^{|i-j|} \\ &\Sigma_{u,v} &=& \begin{bmatrix} 1 & 0.6\\ 0.6 & 1 \end{bmatrix} \\ \end{align*}" /></a>

The invalids IV's are generated:

<a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\begin{align*}&space;F_i&space;&=&space;F^*_i&space;&plus;&space;u_i&space;\times&space;l&space;\\&space;l[j]&space;&=&space;c_l&space;&plus;&space;\dfrac{(0.8&space;-&space;cl)&space;*&space;(j&space;-i)}{19}&space;\\&space;c_l&space;&\in&space;[0,&space;0.8]&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\begin{align*}&space;F_i&space;&=&space;F^*_i&space;&plus;&space;u_i&space;\times&space;l&space;\\&space;l[j]&space;&=&space;c_l&space;&plus;&space;\dfrac{(0.8&space;-&space;cl)&space;*&space;(j&space;-i)}{19}&space;\\&space;c_l&space;&\in&space;[0,&space;0.8]&space;\end{align*}" title="\large \begin{align*} F_i &= F^*_i + u_i \times l \\ l[j] &= c_l + \dfrac{(0.8 - cl) * (j -i)}{19} \\ c_l &\in [0, 0.8] \end{align*}" /></a>

so when cl is close to zero, the IV's in F with smaller index numbers behave more like valid IV's in finite samples.

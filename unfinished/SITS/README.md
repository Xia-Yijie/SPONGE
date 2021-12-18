## 选择性温度积分抽样(SITS)

### 简介

#### 经典SITS

​		部分稀有事件难以在有限的时间内出现，因此需要使用一些增强抽样算法。温度积分抽样（ITS）是一种以温度进行增强抽样的方式。

​		在NVT系综下，某构象的出现几率服从玻尔兹曼分布：
$$
p(\vec r)\propto{\rm exp}(-\beta_0 U(\vec{r}))\tag{SITS-1}\label{SITS-1}
$$
​		其中$\beta_0=1/{kT}$为系综温度，$U(\vec{r})$为坐标为$\vec{r}$时的能量。

​		通过构建新的有效势能，使得
$$
U_{\rm eff}(\vec r)=-\frac{1}{\beta_0}{\rm ln}\sum_k{n_k{\rm exp}\big(-\beta_kU(\vec r)\big)}\tag{SITS-2}\label{SITS-2}
$$
​		此时得到的分布服从下式
$$
p(\vec r)\propto\sum_k{n_k{\rm exp}\big(-\beta_kU(\vec r)\big)}\tag{SITS-3}\label{SITS-3}
$$
​		其中，$\beta_k$是一系列温度，$n_k$作为指前系数使得每个温度的贡献相等，也即
$$
n_1\sum_{\rm all\ samples}{{\rm exp}\big(-\beta_1U\big)}=n_2\sum_{\rm all\ samples}{{\rm exp}\big(-\beta_2U\big)}\\=\cdots=n_k\sum_{\rm all\ samples}{{\rm exp}\big(-\beta_kU\big)}\tag{SITS-4}\label{SITS-4}
$$


​		动力学过程中的力的变化则为
$$
\vec F_{\rm eff}=-\frac{\part U_{\rm eff}(\vec r)}{\part \vec r}=\frac{\part U_{\rm eff}}{\part U}\vec F=\frac{\sum_k{n_k\beta_k{\rm exp}\big(-\beta_kU\big)}}{\beta_0\sum_k{n_k{\rm exp}\big(-\beta_kU\big)}}\vec F=f_b\vec F\tag{SITS-5}\label{SITS-5}
$$
​		其中，
$$
f_b=\frac{\sum_k{n_k\beta_k{\rm exp}\big(-\beta_kU\big)}}{\beta_0\sum_k{n_k{\rm exp}\big(-\beta_kU\big)}}\tag{SITS-6}\label{SITS-6}
$$
​		最后的统计则需要还原每一份样品的概率，需要乘上$\frac{{\rm exp}(-\beta_0 U)}{\sum_k{n_k{\rm exp}(-\beta_kU)}}$的权重。

​		也可以只对体系的部分能量（$U_1+U_2=U$）进行选择性增强，即为SITS，同样地构建有效势能
$$
U_{\rm eff}(\vec r)=-\frac{1}{\beta_0}{\rm ln}\sum_k{n_k{\rm exp}\big(-\beta_kU_1(\vec r)\big)}+U_2(\vec r)\tag{SITS-7}\label{SITS-7}
$$

$$
p(\vec r)\propto\sum_k{n_k{\rm exp}\big(-\beta_kU_1(\vec r)\big)}{\rm exp}\big(-\beta_0 U_2(\vec r)\big)\tag{SITS-8}\label{SITS-8}
$$

$$
\vec F_{\rm eff}=-\frac{\part U_{\rm eff}(\vec r)}{\part \vec r}=-\bigg(\frac{\part U_{\rm eff}}{\part U_1}\frac{\part U_1}{\part \vec r}+\frac{\part U_{\rm eff}}{\part U_2}\frac{\part U_2}{\part \vec r} \bigg)\\ =\frac{\part U_{\rm eff}}{\part U_1}\vec F_1+\frac{\part U_{\rm eff}}{\part U_2}\vec F_2=\frac{\sum_k{n_k\beta_k{\rm exp}\big(-\beta_kU_1\big)}}{\beta_0\sum_k{n_k{\rm exp}\big(-\beta_kU_1\big)}}\vec F_1+\vec F_2\\=f_b\vec F_1+\vec F_2\tag{SITS-9}\label{SITS-9}
$$

​		还原时同样需要乘上$\frac{{\rm exp}(-\beta_0 U_1)}{\sum_k{n_k{\rm exp}(-\beta_kU_1)}}$的权重。通常能量划分时，将体系划分为AB两部分（AB两部分常常是蛋白质和水，也因此使用字母pw）$U_{\rm pp}$、$U_{\rm ww}$为两部分内部的能量，$U_{\rm pw}$为两部分之间的能量，则可以取$U_1=U_{\rm pp}+f_{\rm pw}U_{\rm pw}$，$U_2=U_{\rm ww}+(1-f_{\rm pw})U_{\rm pw}$。此时
$$
\vec F_{\rm eff}=f_b\vec F_{\rm pp}+(f_bf_{\rm pw}+1-f_{\rm pw})\vec F_{\rm pw}+\vec F_{\rm ww}\tag{SITS-10}\label{SITS-10}
$$


​		在实际操作中，有的时候增强效果不够理想或$n_k$很难收敛，并考虑到能量的平移性和叠加性，可以考虑对有效势能的一些部分进行线性变化，使之变为
$$
U_{\rm eff}(\vec r)=-\frac{1}{\beta_0}{\rm ln}\sum_k{n_k{\rm exp}\big(-\beta_k(aU_1+b)\big)}+U_2+cU_1\tag{SITS-11}\label{SITS-11}
$$

$$
p\propto\sum_k{n_k{\rm exp}\big(-\beta_k (aU_1+b)\big)}{\rm exp}\big(-\beta_0 U_2\big){\rm exp}\big(-\beta_0 cU_1\big)\tag{SITS-12}\label{SITS-12}
$$

$$
\vec F_{\rm eff}=\bigg(\frac{\sum_k{n_k\beta_k{\rm exp}\big(-\beta_k (aU_1+b)\big)}}{\beta_0\sum_k{n_k{\rm exp}\big(-\beta_k(aU_1+b)\big)}}+c\bigg)\vec F_1+\vec F_2=f_b'\vec F_1+\vec F_2\tag{SITS-13}\label{SITS-13}
$$

​		还原时需要乘上$\frac{{\rm exp}(-\beta_0 U_1)}{\sum_k{n_k{\rm exp}(-\beta_k(aU_1+b)}{\rm exp}\big(-\beta_0 cU_1\big)}$的权重。

#### 简单SITS

​		上述推导从$U_{\rm eff}$得到$f_b$，反之也可以从$f_b$逆推至$U_{\rm eff}$。令
$$
\vec F_{\rm eff}=f_b\vec F_1+\vec F_2\tag{SITS-14}\label{SITS-14}
$$
​		由此可推得在给定$f_b$下的有效势能和条件概率
$$
U_{\rm eff}=f_bU_1+U_2\tag{SITS-15}\label{SITS-15}
$$

$$
p(\vec r|f_b)\propto{\rm exp}\big(-\beta_0 f_bU_1\big){\rm exp}\big(-\beta_0 U_2\big)\tag{SITS-16}\label{SITS-16}
$$

​		因此可使用蒙特卡洛方法，使得$f_b$服从于某特定的分布$p(f_b)$。例如默认取均匀分布，
$$
p(f_b)=\left\{\begin{array}{**lr**}
0, & f_b<f_{b,\rm min}\\
1/(f_{b,\rm max}-f_{b,\rm min}), &f_{b,\rm max}<f_b<f_{b,\rm min} \\
0, &f_b>f_{b,\rm max}
\end{array}
\right.\tag{SITS-17}\label{SITS-17}
$$
​		也可根据一次体系初探的势能曲线进行设置一个分布曲线。

​		还原时，利用条件概率公式，有
$$
p(\vec r)=\sum_{f_b}p(f_b)\cdot p(\vec r|f_b)
$$

### IO-通用

#### sits_energy_record

记录能量和$f_b$的文件名。文本文件，默认值为`sits_energy_record.txt`。

### IO-经典SITS

#### sits_nk_traj_file

${\rm ln}(n_k)$的迭代轨迹文件名。二进制文件，存储${\rm ln}(n_k)$值，形状为(帧数,`sits_temperature_numbers`)。默认值为`sits_nk_traj_file.dat`。

#### sits_nk_rest_file

${\rm ln}(n_k)$的重开文件名。普通文本文件。默认值为`sits_nk_rest_file.txt`。

#### sits_nk_init_file

${\rm ln}(n_k)$的初始文件名。普通文本文件。无默认值。如果该值未设定，${\rm ln}(n_k)$会被初始化为0.0。如果该值设定，将会读取对应的文件中的值。

#### sits_norm_traj_file

${\rm ln}(W_k)$的迭代轨迹文件名。$W_k$是$n_k$迭代过程的归一化系数。二进制文件，存储${\rm ln}(W_k)$值，形状为(帧数,`sits_temperature_numbers`)。默认值为`sits_norm_traj_file.dat`。

#### sits_norm_rest_file

${\rm ln}(W_k)$的重开文件名。$W_k$是$n_k$迭代过程的归一化系数。普通文本文件。默认值为`sits_norm_rest_file.txt`。

#### sits_norm_init_file

${\rm ln}(W_k)$的初始文件名。$W_k$是$n_k$迭代过程的归一化系数。普通文本文件。无默认值。如果该值未设定，${\rm ln}(W_k)$会被初始化为`-FLT_MAX`。如果该值设定，将会读取对应的文件中的值。

### IO-简单SITS

#### sits_fcball_pdf

$p(f_b)$的离散化值（不用归一化）。如果未设定，则全部设为1.0。设定后会读取对应的文件中的$p(f_b)$值。

### 参数-通用

#### ==sits_mode==

使用的SITS方法。无默认值，必须设定。

- `未设定`：报错。
- `simple`：simple SITS
- `classical`：classical SITS
- `only_select_force`：只对力进行分块计算，不进行增强抽样
- `only_select_force_energy`：只对力和能量进行分块计算，不进行增强抽样

#### ==sits_atom_numbers==

设定前面`sits_atom_numbers`个原子为A部分。未设定会报错。

#### pwwp_enhance_factor

相互作用加强的比例，也即$\eqref{SITS-10}$的$f_{\rm pw}$参数。默认值为0.5。

### 参数-经典SITS

#### sits_temperature_numbers

积分的温度的离散化的份数。默认值为40。

#### sits_temperature_high

积分的最高温度。默认值为NVT系综的T的2倍。

#### sits_temperature_low

积分的最高温度。默认值为NVT系综的T的$1/1.2$。

#### sits_energy_multiple

$\eqref{SITS-11}$中的$a$值。默认值为`1.0`。

#### sits_energy_shift

$\eqref{SITS-11}$中的$b$值。默认值为`0.0`。

#### sits_fb_shift

$\eqref{SITS-11}$中的$c$值。默认值为`0.0`。

#### sits_constant_nk

$n_k$是否进行迭代计算。默认值为`0`

- `0`：进行迭代计算。
- 其他值：不进行迭代计算。

#### sits_record_interval

分子模拟每隔`sits_record_interval`步记录一次能量。默认值为`1`。

#### sits_update_interval

每记录`sits_update_interval`步能量迭代一次$n_k$。默认值为`100`。

### 参数-简单SITS

#### sits_fcball_pdf_grid_numbers

$f_b$参数的概率分布曲线离散化的格点数目。默认值为`1000`。

#### sits_constant_fcball

如果设定，将$f_b$固定为`sits_constant_fcball`值不进行蒙特卡洛模拟。默认未设定。

#### sits_fcball_max

$\eqref{SITS-17}$中的$f_{b,\rm max}$。默认值为`1.2`。

#### sits_fcball_min

$\eqref{SITS-17}$中的$f_{b,\rm min}$。默认值为`0.5`。

#### sits_fcball_move_length

每次$f_b$蒙特卡洛模拟过程中，尝试移动的可能的最长距离。

#### sits_fcball_random_seed

$f_b$蒙特卡洛模拟过程的随机数种子。
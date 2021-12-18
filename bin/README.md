# 北大高组分子动力学模拟程序SPONGE使用说明

[TOC]

SPONGE的控制命令均使用mdin中`Flag = Value`或命令行中`-Flag Value`的形式，两者完全等效，例如想设置模拟的总步长为10000步，其他命令均读取文件`mdin`，则可以在命令行中如下启动`SPONGE`:

```shell
SPONGE -i mdin -step_limit 10000
```

也可以在命令行中以`i`指定输入文件，然后在输入文件中填写`step_limit = 10000`

```shell
SPONGE -i mdin
```

如果某个`Flag`被两种方式重复赋值，将会报错。
如果某个`Flag`被设定了但是未被使用，将会报警。
下面将按照模块介绍各个`Flag`以及其对应的`Value`的功能。

==高亮==的`Flag`是被认为比较重要需要修改的参数。

## 控制(control)

### 简介

负责命令控制和输入输出的模块。

`SPONGE`的输入模块化，使用模块化的若干文本文件。

使用`AMBER`的输入（parm_file_name、rst7_file_name）跑normal MD，输入使用下列命令：

```shell
SPONGE -i mdin -amber_parm parm_file_name -c rst7_file_name 
```

使用`GROMACS`的输入跑normal MD，使用`tools`中的`input2sponge.py`转换

### IO

#### i

输入文件的文件名。默认值为`mdin`。
输入文件中，第一行是本次模拟的任务名，其余各行均以`Flag = Value`的形式输入命令。
可以使用英文半角逗号`,`或换行分隔各个命令。
每行行首或`,`后的英文半角`/`、`#`、`!`是注释符号，该行后面的内容将会被忽视
没有`=`的行将会被视为注释而被忽视。可以通过该方式书写多行注释。

#### c

模拟进行的初始坐标的文件名。没有默认值。
按照`SPONGE`的坐标文件格式读取坐标。
`SPONGE`坐标文件第一个数是系统原子数，后面是每个原子的三维坐标，最后是周期性盒子的长度和角度。
如果`amber_irest != -1`，则会以`AMBER`的ascii的rst7文件格式读取。
`GROMACS`的坐标文件请使用`tools/sponge2input.py`进行转换

#### v0

模拟进行的初始速度的文件名。没有默认值。
按照`SPONGE`的坐标文件格式读取速度。
`SPONGE`速度文件第一个数是系统原子数，后面是每个原子的三维速度。
如果`amber_irest != -1`，则该命令不会被使用。
如果`v0`没被定义，则所有原子的速度会被初始化为0。

#### r

模拟输出的坐标和速度的文件名。默认值为`restrt`。
按照`SPONGE`的坐标文件格式输出坐标到文件`r`\_coordinate.txt。
按照`SPONGE`的速度文件格式输出速度到文件`r`\_velocity.txt。
如果`amber_irest != -1`，则会以`AMBER`的ascii的rst7文件格式输出坐标和速度到文件`r`。

#### o

模拟输出信息的文件名。默认值为`mdout`。
第一行会显示每行输出信息的内容。
每`write_information_interval`步输出一次信息。

#### x

模拟输出的坐标轨迹文件名。默认值为`mdcrd`。
输出为3×原子个数×帧数的float二进制文件。
c/c++使用fread读取。
python使用np.fromfile(`x`, dtype=np.float32)读取。

#### vx

模拟输出的速度轨迹文件名。没有默认值。
如果没有设定，则不会输出速度轨迹。

#### box

模拟输出的盒子信息轨迹文件名。默认值为`mdbox`。

#### amber_parm

模拟需要读取的`AMBER`格式的参数文件名。
如果同时设定了`SPONGE`本身格式的参数文件，将会优先以`SPONGE`本身格式的参数文件初始化。

### 参数

#### amber_irest

用以兼容`AMBER`的坐标/速度输入输出命令。无默认值。

* 未设定：使用`SPONGE`的格式输入输出坐标、速度文件
* `0`：使用`AMBER`格式输入输出坐标、速度文件，且输入坐标文件不包含速度
* `1`：使用`AMBER`格式输入输出坐标、速度文件，且输入坐标文件包含速度

#### dont_check_input

通常来说，control会检查每个命令是否都被使用，以免输入错误或其他情况。可以设置dont_check_input = 1使得control模块不对此进行检查。

#### ==device==

运行的GPU设备编号。默认值为`0`。



## 核心(MD_core)

### 简介

md运行的核心模块，包括各个核心计算量和迭代。

### IO

#### mass_in_file

模拟需要读取的质量参数的文件名。
格式为第一个数为原子个数，然后是各个原子的相对原子质量。

#### charge_in_file

模拟需要读取的电荷量参数的文件名。
格式为第一个数为原子个数，然后是各个原子的电荷量。单位为$\frac{1}{18.2223}$元电荷。

#### residue_in_file

模拟需要读取的残基列表的文件名。
格式为第一个数为原子个数，第二个数为残基个数，然后是各个残基的原子个数。

### 参数

#### ==mode==

分子模拟的模式。默认值为`0`。

* `-1`: 梯度下降最小化。
  此时`dt`的平方为梯度下降的学习率。
* `0`: NVE系综的分子模拟
* `1`: NVT系综的分子模拟
* `2`: NPT系综的分子模拟

#### ==dt==

当`mode >= 0`时，为分子动力学模拟的每一步的时间长度；
当`mode == -1`时，为梯度下降的学习率的平方根。
单位为ps，默认值为`0.001`。支持`1e-3`之类的科学计数法表示。

#### ==step_limit==

分子模拟的总步数。默认值为`1000`。

#### ==write_information_interval==

模拟每隔`write_information_interval`步输出一次信息。默认值为`1000`。

#### calculate_energy_interval

模拟每隔`calculate_energy_interval`步额外计算一次能量。无默认值。
输出信息的时候会默认计算一次能量。

#### net_force

是否需要将合力手动归零。默认值为`0`。
因为存在计算误差，因此系统的合力可能不满足牛顿第三定律而导致合力不为0。
此时可以对所有原子添加一个力使得合力为0。

* `0`: 合力不归零。
* `1`: 归零。



## 键长项(bond)

### 简介

定义为：
$$
E=k(r-r_0)^{2}
$$
其中，

| 符号  |   物理量   |             单位             |
| :---: | :--------: | :--------------------------: |
|  $E$  |    能量    |     $\rm kcal·mol^{-1}$      |
|  $k$  |   力常数   | $\rm kcal·mol^{-1}·\AA^{-2}$ |
|  $r$  | 原子间距离 |          $\rm \AA$           |
| $r_0$ |  平衡距离  |          $\rm \AA$           |

### IO

#### bond_in_file

模拟需要读取的体系的所有的键长参数的文件名。
文件第一个数为键长项的个数`bond_numbers`。然后`bond_numbers`行，每行为计算键长的原子序号$a$、原子序号$b$、力常数$k$、平衡距离$r_0$。



## 键角项(angle)

### 简介

定义为：
$$
E=k(\theta-\theta_0)^{2}
$$
其中，

|    符号    |   物理量    |             单位             |
| :--------: | :---------: | :--------------------------: |
|    $E$     |    能量     |     $\rm kcal·mol^{-1}$      |
|    $k$     |   力常数    | $\rm kcal·mol^{-1}·rad^{-2}$ |
|  $\theta$  | $\ang{bac}$ |          $\rm rad$           |
| $\theta_0$ |  平衡角度   |          $\rm rad$           |

### IO

#### angle_in_file

模拟需要读取的体系的所有的键角参数的文件名。

文件第一个数为键长项的个数`angle_numbers`。然后`angle_numbers`行，每行为计算键长的原子序号$a$、原子序号$b$、原子序号$c$、力常数$k$、平衡角度$\theta_0$。

原子的连接顺序为a-b-c。

## 二面角项(dihedral)

### 简介

定义为：
$$
E=k(1+{\rm cos}(n\phi-\phi_0))
$$
其中，

|   符号   |            物理量            |        单位         |
| :------: | :--------------------------: | :-----------------: |
|   $E$    |             能量             | $\rm kcal·mol^{-1}$ |
|   $k$    |            力常数            | $\rm kcal·mol^{-1}$ |
|   $n$    |         二面角项周期         |          /          |
|  $\phi$  | 平面abc与平面bcd形成的二面角 |      $\rm rad$      |
| $\phi_0$ |           平衡角度           |      $\rm rad$      |

### IO

#### dihedral_in_file

模拟需要读取的体系的所有的二面角参数的文件名。

文件第一个数为键长项的个数`dihedral_numbers`。然后`dihedral_numbers`行，每行为计算键长的原子序号$a$、原子序号$b$、原子序号$c$、原子序号$d$、二面角项周期$n$（整数型）、力常数$k$、平衡角度$\phi_0$。

原子的连接顺序为a-b-c-d。

## 非键14项(nb14)

### 简介

定义为：
$$
E=k_{\rm LJ}(\frac{A_{ij}}{r^{12}}-\frac{B_{ij}}{r^{6}})+k_{\rm CF}\frac{q_{i}q_{j}}{r}
$$
其中，

|        符号        |                            物理量                            |        单位         |
| :----------------: | :----------------------------------------------------------: | :-----------------: |
|        $E$         |                             能量                             | $\rm kcal·mol^{-1}$ |
|    $k_{\rm LJ}$    |                 Lennard-Jones 14作用修正系数                 |          /          |
|    $k_{\rm CF}$    |                    Coulomb 14作用修正系数                    |          /          |
|        $r$         |                          原子间距离                          |        $\AA$        |
| $A_{ij}$、$B_{ij}$ | 原子$i$和原子$j$之间的Lennard-Jones作用参数。调用LJ模块读取的值 |          /          |
|        $q$         |                 电荷量。调用MD_core读取的值                  |          /          |

### IO

#### nb14_in_file

模拟需要读取的体系的所有的非键14参数的文件名。

文件第一个数为键长项的个数`nb14_numbers`。然后`nb14_numbers`行，每行为计算键长的原子序号$a$、原子序号$b$、Lennard-Jones 14作用修正系数$k_{\rm LJ}$，Coulomb 14作用修正系数$k_{\rm CF}$。

## 近邻表(Neighbor_List)

### 简介

对于非键作用，基于格点法构建近邻表，以此简化运算。

### IO

#### exclude_in_file

排除表的文件名。第一行两个参数，分别为原子个数`atom_number`和排除表总长。

下面`atom_number`行，每行第一个数为该原子的排除原子个数，然后接的是各个排除的原子。

注意，排除时a排除b和b排除a只写一遍。

### 参数

#### ==neighbor_list_refresh_interval==

每隔`neighbor_list_refresh_interval`步更新一次近邻表。默认值为`20`。

如果设置为0，则每当有原子移动超过了`skin`距离，则会立刻刷新近邻表。

#### ==cut==

非键计算的截断距离。单位为埃。默认值为`10.0`。不能大于周期性盒子长度的一半。

#### skin

两次刷新近邻表之间，一个原子不会移动超过的距离。单位为埃。默认值为`2.0`。

#### max_atom_in_grid_numbers

每个格点中最多有`max_atom_in_grid_numbers`个原子。默认值为``。当格点中的原子数超过该值时可能会出现错误。

#### max_neighbor_numbers

每个原子最多有`max_neighbor_numbers`个近邻。默认值为`800`。当一个原子的近邻超过该数目后可能会出现错误。

## LJ项(Lennard_Jones)

### 简介

定义为：
$$
E=\frac{A_{ij}}{r^{12}}-\frac{B_{ij}}{r^{6}}
$$
其中，

|   符号   |   物理量   |             单位              |
| :------: | :--------: | :---------------------------: |
|   $E$    |    能量    |      $\rm kcal·mol^{-1}$      |
| $A_{ij}$ | LJ排斥系数 | $\rm kcal·mol^{-1}·\AA^{-12}$ |
| $B_{ij}$ | LJ吸引系数 | $\rm kcal·mol^{-1}·\AA^{-6}$  |
|   $r$    | 原子间距离 |             $\AA$             |

### IO

#### lj_in_file

模拟需要读取的体系的所有的LJ参数的文件名。

文件第一个数为原子个数`atom_numbers`，第二个数为LJ类型数`atom_type_numbers`。

接下来的`atom_type_numbers`行给出$A_{ij}$的值。这`atom_type_numbers`行中，第$i$行的数字共$i$个。其中第$i$行第$j$列即为原子LJ类型$i$和原子LJ类型$j$之间的$A_{ij}$。

接下来的`atom_type_numbers`行给出$B_{ij}$的值。这`atom_type_numbers`行中，第$i$行的数字共$i$个。其中第$i$行第$j$列即为原子LJ类型$i$和原子LJ类型$j$之间的$B_{ij}$。

接下来的`atom_numbers`行给出每个原子的LJ类型。(注意类别从0开始计数)

## PME项(PME)

### 简介

使用Particle Mesh Ewald方法处理的静电能。静电能的定义为：
$$
E=\frac{q_{i}q_{j}}{r}
$$
其中，

| 符号 |   物理量   |        单位         |
| :--: | :--------: | :-----------------: |
| $E$  |    能量    | $\rm kcal·mol^{-1}$ |
| $q$  |   电荷量   | $\frac{1}{18.223}e$ |
| $r$  | 原子间距离 |        $\AA$        |

注意电荷量的单位为$\frac{1}{18.223}e$，也即18.223的SPONGE电荷量为一个电子的电荷量。使用该数值是为了使得点电荷静电能的公式能够直接等于电荷量之积除以距离，常数项系数为1。

### 参数

#### PME_Direct_Tolerance

直接部分的误差容忍度。默认值为`0.00001`。高斯电荷的$\beta$值迭代计算得${\rm erfc}(\beta r)/{r}=$`PME_Direct_Tolerance`。

#### fftx、ffty、fftz

傅里叶变换的三个维度的大小。

如果未设定，则会选取同时满足下列条件最小的整数：

- 大于等于周期性盒子对应维度的长度
- 4的倍数
- 除了2、3、5、7不能有其他质因数。

其中后两个条件是便于傅里叶变换加速

## 限制项(restrain)

### 简介

使用谐振子将原子限制在参考位置。被限制的原子会额外受到一个回复力。回复力的具体形式受`restrain_mode`调控。

### IO

#### ref

参考坐标的文件名。与`c`相同的读取方式。

#### reference_mask

被限制的原子的文件名。文件中是所有被限制的原子的原子序号。注意：原子序号是从0开始计数的。

可以使用`tools`中的`maskgen.py`调用`VMD`按照`VMD`的原子筛选方式生成对应的mask。

### 参数

#### ==restrain_mode==

限制的回复力形式。无默认值。

* 未设定：不启用限制。
* `0`：点限制。$\vec{F}=-k(\vec r-\vec r_{\rm ref})$
* `1`：面限制。$F_z=-k(z-z_{\rm ref})$

#### restrain_weight

限制回复力的力常数$k$的大小。单位$\rm kcal·mol^{-1}·\AA^{-2}$。默认值为`100.0`。

## 约束项(constrain)

### 简介

在迭代中加入约束方程，用于固定部分键长和键角。具体的约束算法由`constrain_mode`调控。

### 总参数

#### ==constrain_mode==

约束的算法。无默认值

* 未设定：不启用约束。

* `0`：简单约束算法。将所有包含有相对原子质量小于`constrain_mass`的原子的键长(bond)和键角(angle)进行固定。对键角的固定使用余弦定理转化为对键长的固定。在正常的迭代过后，对有约束键长的原子$i$坐标按照下列迭代式
  $$
  \vec{r}_i'=\vec{r}_i+\sum_{j}{k\frac{m_j}{m_i+m_j}\Delta t^2(1-\frac{r_0}{r_{ij}})\vec{r}_{ij,0}}
  $$


  迭代`simple_constrain_interation_numbers`次。迭代式中，$\vec{r}_i$是原子$i$的坐标，求和号针对所有与原子$i$存在约束键长的原子$j$，$m_{i}$是原子$i$的质量，$\vec{r}_{ij,0}$是正常迭代前的距离向量，$r_{ij}=|\vec{r}_i-\vec{r}_j|$，$r_0$是约束键长的长度，$\Delta t$是积分步长，$k$是简单约束算法迭代步长的系数，由`simple_constrain_step_length`控制。只能用于能量最小化、朗之万控温和朗之万-刘剑控温中。

### 总参数

#### constrain_mass

将对质量大于0小于`constrain_mass`的原子的键长和键角进行约束。当未定义了`constrain_in_file`时，默认值为`3.0`，否则默认值为`0.0`

#### constrain_in_file

描述哪些键需要进行约束的文件。文件每一行定义一个需要约束的键长，每一行三个数字，分别是需要约束的两个原子和约束到的长度。

### 简单约束参数

#### simple_constrain_step_length

简单约束算法的迭代步长的系数$k$。默认值为`1.0`。

`simple_constrain_step_length`值越小，越不容易发散，但是收敛需要的迭代次数相应会增加。

#### simple_constrain_iteration_numbers

简单约束算法的迭代次数。默认值为`50`。

## 热浴(thermostat)

### 简介

当`mode >= 1`时，会对系统的温度进行控制。具体的控温算法由使用的热浴`thermostat`决定。

### 总参数

#### ==thermostat==

使用的热浴类型。默认值为`0`。

- `0`：朗之万热浴
- `1`：朗之万-刘剑热浴

#### ==target_temperature==

目标温度。单位为开尔文。默认值为`300.0`

### 朗之万热浴参数

#### ==langevin_gamma==

朗之万动力学中的摩擦因子，又称为碰撞频率。单位为$\rm ps^{-1}$。默认值为`1.0`。

`langevin_gamma`越大，与热浴的耦合越强。

当`langevin_gamma = 0`时，不会控温。

#### langevin_seed

朗之万动力学中使用的伪随机数种子。如果未设定，默认值为以时间为种子伪随机选取的一个整数。

### 朗之万-刘剑热浴参数

#### ==langevin_gamma==

朗之万动力学中的摩擦因子，又称为碰撞频率。单位为$\rm ps^{-1}$。默认值为`1.0`。

`langevin_gamma`越大，与热浴的耦合越强。

当`langevin_gamma = 0`时，不会控温。

#### langevin_seed

朗之万动力学中使用的伪随机数种子。如果未设定，默认值为以时间为种子伪随机选取的一个整数。

#### velocity_max

迭代中允许出现的最大速度。默认值无。未设定时不会有影响，设定后每次使用迭代的速度超过该值将会强行使得速度等于该值。

## 压浴(barostat)

### 简介

当`mode = 2`时，会对系统的压强进行控制。具体的控温算法由使用的压浴`barostat`决定。

### 总参数

#### barostat

使用的压浴类型。默认值为`0`。

- `0`：蒙特卡洛控压
- `1`：贝伦德森控压

#### ==target_pressure==

目标压强。单位为bar。默认值为`1.0`。

### 蒙特卡洛控压参数

#### ==mc_baro_update_interval==

每隔`mc_baro_update_interval`步进行一次蒙特卡洛体积变换尝试。默认值为`100`。如果模拟中体积变化不明显，可以将该值调低，但会减慢模拟速度。

#### mc_baro_initial_ratio

初始的最大变化体积占总体积的比例。默认值为`0.01`，也即每次变化最大可能变化体积的1%。

#### mc_baro_check_interval

每隔`mc_baro_check_interval`次更新后对最大变化体积进行迭代。收敛的最大变化体积应使得蒙特卡洛接受率在一个合理的区间，该区间由下面两个参数确定。

#### mc_baro_accept_rate_low

使得收敛的最大变化体积应使得蒙特卡洛接受百分数的低值，默认值为`30.0`。

#### mc_baro_accept_rate_high

使得收敛的最大变化体积应使得蒙特卡洛接受百分数的高值，默认值为`40.0`。

### 贝伦德森控压参数

#### bd_baro_tau

弛豫时间。单位为皮秒。默认值为`1.0`。该值不会影响平衡，会影响部分体积变化速度。过大可能会导致较大的体积振荡，过小可能会影响控压效率。

#### bd_baro_compressibility

体系的压缩系数。单位为${\rm bar^{-1}}$。默认值为水的$4.5\times10^{-5}$。该值不会影响平衡，会影响部分体积变化速度。

## 虚原子(virtual atoms)

### 简介

部分模拟中，会构建出虚原子，虚原子的位置完全由其他原子（可以是实原子，也可以是虚原子）的位置产生，而且虚原子通常的质量也为0。为保持能量守恒，虚原子受到的力需要转移到实原子上。此时通过公式$\eqref{vatoms-1}$重新分配虚原子上的力：
$$
\begin{align}
\vec F_{\rm real}(\vec r_{\rm real},\vec r_{\rm virtual})
&=-\frac {\part U(\vec r_{\rm real},\vec r_{\rm virtual})}{\part \vec r_{\rm real}}
\\&=\vec F(\vec r_{\rm real}|\vec r_{\rm virtual})-\frac {\part U(\vec r_{\rm virtual}|\vec r_{\rm real})}{\part \vec r_{\rm virtual}} \frac {\part \vec r_{\rm virtual}}{\part \vec r_{\rm real}}
\\&=\vec F(\vec r_{\rm real}|\vec r_{\rm virtual})+\frac {\part \vec r_{\rm virtual}}{\part \vec r_{\rm real}}\vec F(\vec r_{\rm virtual}|\vec r_{\rm real})
\end{align}\tag{vatoms-1}\label{vatoms-1}
$$
其中，$\vec F_{\rm real}(\vec r_{\rm real},\vec r_{\rm virtual})$表示实原子受到的总的力，$U(\vec r_{\rm real},\vec r_{\rm virtual})$表示系统总能量，$\vec F(\vec r_{\rm virtual}|\vec r_{\rm real})$和$U(\vec r_{\rm virtual}|\vec r_{\rm real})$表示实原子坐标作为参数，虚原子受到的力和能量，$\vec F(\vec r_{\rm real}|\vec r_{\rm virtual})$和$U(\vec r_{\rm real}|\vec r_{\rm virtual})$表示虚原子坐标为参数，实原子受到的力和能量。

虚原子的构建可以参照于其他虚原子，此时原子有多个虚层级。本简介部分下文中的角标$i$表示第$i$层原子。第$i$层虚原子的坐标依赖于至少一个第$(i-1)$层原子和若干个小于第$(i-1)$层的原子。实原子定义为第0层原子。虚原子上的力从最高层递归地重新分配至实原子上：
$$
\begin{align}
\vec F_{0,\cdots,n-1}(\vec r_0,\vec r_1, \cdots,\vec r_n|\vec\sigma) &= -\frac {\part U(\vec r_0,\vec r_1, \cdots,\vec r_n|\vec\sigma)}{\part (\vec r_0,\vec r_1, \cdots,\vec r_{n-1})}
\\&=\vec F(\vec r_0, \cdots,\vec r_{n-1}|\vec r_n,\vec\sigma)-\frac {\part U(\vec r_n|\vec r_0, \cdots,\vec r_{n-1},\vec\sigma)}{\part \vec r_n}\frac {\part \vec r_n}{\part (\vec r_0,\vec r_1, \cdots,\vec r_{n-1})}
\\&=\vec F(\vec r_n|\vec r_0, \cdots,\vec r_{n-1},\vec\sigma)+\frac {\part \vec r_n}{\part (\vec r_0,\vec r_1, \cdots,\vec r_{n-1})}\vec F(\vec r_0, \cdots,\vec r_{n-1}|\vec r_n,\vec\sigma)
\end{align}\tag{vatoms-2}\label{vatoms-2}
$$
其中$\vec\sigma$表示其他参数，包含温度、电荷、LJ半径等参数，也包含$\vec r_{n+1}$等更高层虚原子的坐标。逗号连接的矢量为对应的直和空间的矢量。$|$前的量为变量，$|$后的量为参数。

### IO

#### ==virtual_atom_in_file==

定义虚原子信息的文件。

注意：如果某虚原子依赖于其他虚原子，那么该虚原子一定要在依赖的虚原子之后定义。

每行定义一个虚原子的信息。每行的前两个数分别是虚原子的函数类型和虚原子的编号。其余数字由虚原子类型决定。

下文中$\vec r_{ij}=\vec r_i-\vec r_j$

虚原子类型0：
$$
\vec r_{\rm v}=
\begin{bmatrix} 
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 2h-1\\
\end{bmatrix} \vec r_1
\label{vatoms-3}\tag{vatoms-3}
$$

$$
\vec F_1=
\begin{bmatrix} 
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & -1\\
\end{bmatrix} \vec F_{\rm v}
\label{vatoms-4}\tag{vatoms-4}
$$



​		后面跟1号父原子编号和$h$值。

​		输入例子：0 1 0 10.0表示1号原子是虚原子类型0，以$\eqref{vatoms-3}$的形式依赖于0号原子，$\eqref{vatoms-3}$中的$h=10.0$。

​        用于镜像电荷，一般可能会限制$z$坐标。

虚原子类型1：
$$
\vec r_{{\rm v}1}=
a
\vec r_{21}
\label{vatoms-5}\tag{vatoms-5}
$$

$$
\begin{align}
\vec F_1 &= (1-a)\vec F_{\rm v} \\
\vec F_2 &= a \vec F_{\rm v}
\end{align}
\label{vatoms-6}\tag{vatoms-6}
$$

​		后面跟1号父原子编号，2号父原子编号和$a$值。

​		输入例子：1 2 0 1 0.3表示2号原子是虚原子类型1，以$\eqref{vatoms-5}$的形式依赖于0号原子和1号原子，其中$a=0.3$。

​		虚原子和父原子共线。

虚原子类型2：
$$
\vec r_{{\rm v}1}=
a \vec r_{21} + b \vec r_{31}
\label{vatoms-7}\tag{vatoms-7}
$$

$$
\begin{align}
\vec F_1 &= (1-a-b)\vec F_{\rm v} \\
\vec F_2 &= a \vec F_{\rm v} \\
\vec F_3 &= b \vec F_{\rm v} \\
\end{align}
\label{vatoms-8}\tag{vatoms-8}
$$

​		后面跟1号父原子编号，2号父原子编号，3号父原子编号，$a$值和$b$值。

​		虚原子和父原子共平面。

虚原子类型3：
$$
\vec r_{{\rm v}1}=
d \frac {\vec r_{21} + k \vec r_{32}} {|\vec r_{21} + k \vec r_{32}|}
\label{vatoms-9}\tag{vatoms-9}
$$

$$
\begin{align}
\vec F_1 &= \vec F_{\rm v}-\vec F_{\rm p}\\
\vec F_2 &= (1-k) \vec F_{\rm p} \\
\vec F_3 &= k \vec F_{\rm p} \\
\vec F_{\rm p} &= \frac d {|\vec r_{21} + k \vec r_{32}|} \frac {\vec r_{{\rm v}1}·\vec F_{\rm v}} {\vec r_{{\rm v}1}·\vec r_{{\rm v}1}}\vec r_{{\rm v}1}
\end{align}
\label{vatoms-10}\tag{vatoms-10}
$$

​		后面跟1号父原子编号，2号父原子编号，3号父原子编号，$d$值和$k$值。

​		虚原子和父原子共平面，且到1号父原子的距离始终为$d$。

## Python接口(pyplugin)

### 简介

SPONGE的Python语言的接口模块，接入后可以使用Python脚本对模拟进行修改。除了依赖于Python的标准库以外，还使用了numpy库，因此需要使用

```shell
pip install numpy
```

进行安装。

### 模块安装

### Python环境

建议使用Python3.7及以上的版本

#### Linux

使用`sponge_builder.py`构建主程序和Makefile后，修改该`Makefile`中的下列部分

```makefile
PYTHON=python3.8
PYTHON_INCLUDE=/home/`whoami`/anaconda3/include/python3.8
NUMPY_INCLUDE=/home/`whoami`/anaconda3/lib/python3.8/site-packages/numpy/core/include
PYTHON_LIBRARY=/home/`whoami`/anaconda3/lib
```

改为你系统的Python动态库、include文件夹、numpy的include文件夹和Python的动态库文件夹。

推荐使用conda来管理python。

可使用下列python命令获得对应的文件夹。

```python
import sysconfig
print(sysconfig.get_config_var('LIBDIR')) #lib文件夹，不过因为不同python编译不同，该处可能不会给对应的结果
print(sysconfig.get_paths()["include"]) #include文件夹
import numpy
print(numpy.__path__[0]+"/core/include") #numpy的include文件夹
```

另外一种方式是使用shell命令查找

首先使用shell命令进入你的python交互模式

```python
import sysconfig
sysconfig.get_config_var('prefix')
```

该命令输出的即为你的shell命令使用的python的文件夹地址。在该文件夹地址中利用

```shell
find python的文件夹地址 -name Python.h
find python的文件夹地址 -name libpython*.so
```

 找到对应的python的include文件和动态库文件，然后将对应的文件填入上面的地址中

#### Windows

windows主要使用visual studio进行编译。

同样先进入python的交互模式获取python的文件夹地址。

```python
import sysconfig
sysconfig.get_config_var('prefix')
```

windows的include文件一般在python的文件夹地址\include中

lib文件在python的文件夹地址\libs中

在`项目`-`属性`-`配置属性`-`VC++`目录中的包含目录中加入

```
python的文件夹地址\include
python的文件夹地址\Lib\site-packages\numpy\core\include\numpy
```

在`项目`-`属性`-`配置属性`-`VC++`目录中的库目录中加入

```
python的文件夹地址\libs
```

### IO

#### ==py==

在命令行以`py = 文件名`的形式或终端中以`-py  文件名`的方式确定在主程序运行中的Python脚本名。

### Python端API

#### 简介

SPONGE的主程序类似一个Python主程序，在初始化时会import `py`指定的文件。

`src/quantum`文件夹中的`pyscf.py`、`g09.py`和`bdf.py`均调用了Python接口进行AIMD。可以查看这些文件作为示例。

在Python脚本中，先需要

```Python
import sponge
```

然后在函数定义时加入修饰符`@sponge.register(xxx)`注册

```Python
@sponge.register("Calculate_Force")
def funcA():
    crd = sponge.get_coordinate()
    sponge.set_coordinate(crd)
#Volume_Change的函数形式需要接受一个参数，该参数是坐标标度的系数，如所有坐标和盒子的边长都会乘0.9，则会执行funcWithAnyNameYouWant(0.9)
@sponge.register("Volume_Change")
def funcWithAnyNameYouWant(factor):
    frc = sponge.get_force()
    sponge.set_force(frc)
```

就会在对应的时候调用该函数。在函数内使用sponge.get_xxx获取对应的值，并以sponge.set_xxx将该值传回主程序。



#### 支持的修饰符

| 修饰符                                     | 作用                                                         |
| ------------------------------------------ | ------------------------------------------------------------ |
| @sponge.register("Iteration_1")            | 每步的初始迭代时调用。如速度Verlet算法会在每步开始时进行迭代。 |
| @sponge.register("Before_Calculate_Force") | 每步计算力之前调用。如将坐标转化为分数坐标、决定每步是否需要计算维里或原子能量。 |
| @sponge.register("Calculate_Force")        | 每步计算力的时候调用。                                       |
| @sponge.register("Calculate_Energy")       | 在需要计算能量的步骤计算每个分能量的时候调用。               |
| @sponge.register("After_Calculate_Energy") | 在需要计算能量的步骤计算每个分能量之后调用。主要是将能量从GPU下载，并将能量加和到总势能上。 |
| @sponge.register("After_Calculate_Force")  | 每步计算力之后的时候调用。如对力增强抽样、取消Before_Calculate_Force决定的计算维里或原子能量。 |
| @sponge.register("Before_Iteration")       | 坐标迭代之前调用。如简单约束会在此步获取坐标。               |
| @sponge.register("After_Iteration")        | 坐标迭代之后调用。如简单约束对坐标更改，压强调整坐标。       |
| @sponge.register("Iteration_2")            | 每步的迭代时调用。如速度Verlet算法后续的迭代，蛙跳法的迭代。 |
| @sponge.register("Volume_Change")          | 体积变化时调用。注意该函数必须接受一个浮点数作为坐标调整的系数。 |
| @sponge.register("Print")                  | 打印信息时调用。                                             |
| @sponge.register("Destroy")                | 程序结束时调用。如释放内存等。                               |

#### 支持的函数

- sponge.get_atom_numbers()
  参数：无
  返回值：int atom_numbers 分子模拟的原子数
- sponge.get_step()
  参数：无
  返回值：int step 当前的模拟步数
- sponge.get_potential_energy()
  参数：无
  返回值：float potential_energy 总势能
  注意：模拟并不是每一步都计算了总势能，确保在计算能量的步骤计算能量之后调用获得正确的值。
- sponge.set_potential_energy(energy)
  参数：float energy  设置的势能值
  返回值：无
- sponge.get_box_length()
  参数：无
  返回值：numpy.ndarray box_length 周期性盒子的三个边长组成的1维numpy数组
- sponge.get_mass(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray mass 第start号原子到第start + length - 1号原子的质量组成的1维numpy数组
- sponge.get_charge(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray charge 第start号原子到第start + length - 1号原子的电荷组成的1维numpy数组
- sponge.set_charge(charge, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray charge 设置的电荷
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.get_coordinate(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray coordinate 第start号原子到第start + length - 1号原子的坐标组成的2维numpy数组
- sponge.set_coordinate(coordinate, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray coordinate 设置的坐标
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.get_force(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray force第start号原子到第start + length - 1号原子的力组成的2维numpy数组
- sponge.set_force(force, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray force 设置的力
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.get_label(start = 0, length = sponge.get_atom_numbers())
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
  返回值：list label第start号原子到第start + length - 1号原子的标签组成的1维python list
- sponge.set_label(label, start = 0, length = sponge.get_atom_numbers())
  参数：
          list label 设置的标签
          int start 起始原子编号
          int length 需要获得的总原子数
  返回值：无
- sponge.need_atom_energy()
  参数：无
  返回值： bool need 是否需要计算原子能量
- sponge.get_atom_energy(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray atom_energy 第start号原子到第start + length - 1号原子的能量组成的1维numpy数组
- sponge.set_atom_energy(atom_energy, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray atom_energy 设置的原子能量
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.need_virial()
  参数：无
  返回值： bool need 是否需要计算维里
- sponge.get_atom_virial(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：numpy.ndarray atom_virial第start号原子到第start + length - 1号原子的维里组成的1维numpy数组
- sponge.set_atom_virial(atom_virial, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray atom_virial设置的原子维里
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.get_virial(start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否从gpu上下载
  返回值：float virial 系统的总标量维里
- sponge.set_virial(virial, start = 0, length = sponge.get_atom_numbers(), cuda_memcpy = False)
  参数：
          numpy.ndarray atom_virial 设置的系统的总标量维里
          int start 起始原子编号
          int length 需要获得的总原子数
          bool cuda_memcpy 是否上传至gpu上
  返回值：无
- sponge.command_exist(command)
  参数：str command 命令
  返回值：bool existance 命令存在与否
- sponge.command(command)
  参数：str command 命令
  返回值：str value 命令对应的除去前后空格的第一个值。
  注意：
          如果命令不存在也不会报错，可能会给出错误的值。
          例如mdin中设置，“mycommand = what a nice day ”，sponge.command(mycommand)得到的返回值是“what”。
          例如命令行中输入，“-test haha hahaha ”，sponge.command(test)得到的返回值是"haha"
- sponge.origin_command(command)
  参数：str command 命令
  返回值：str value 命令对应的未除去前后空格的原本值。
  注意：
          如果命令不存在也不会报错，可能会给出错误的值。
          例如mdin中设置，“mycommand = what a nice day  ,”sponge.origin_command(mycommand)得到的返回值是“ what a nice day  ”。
          例如命令行中输入，“-test haha hahaha ”，sponge.origin_command(test)得到的返回值是" haha hahaha"
- sponge.add_print_head(head)
  参数：str head 需要增加的表头。调用后，每步print的表头将额外打出该内容。只需用于初始化。
  返回值：无
- sponge.add_print(content)
  参数：str content 需要打印的内容。调用后，该步print内容中会增加该内容，也会存至mdout文件中。在`Print`步骤调用。

## 量子力学(quantum mechanics)

### 简介

可以使用python脚本调用外部的量子化学程序，以实现AIMD或QM/MM方法。程序内置了一些简单的样例，可以实现一些基础的功能。

### 参数

#### ==qm_program==

使用的外部量化程序。目前支持pyscf、Gaussian和Beijing Density Function(BDF)。无默认值。

- 未设定或`0`：不启用
- `1`：使用pyscf
- `2`：使用Gaussian
- `3`：使用BDF

#### ==qm_command==

使用的外部量化程序的命令行命令。

对于pyscf无效，因为pyscf是直接调用Python。

对于Gaussian，默认值为g09。

对于BDF，无默认值。参考值为$BDFHOME/sbin/run.sh。

注意：因为是用Python调用shell，不能以~/.bashrc中使用alias等效的命令作为命令行命令。

#### ==qm_atoms==

参与量化计算的原子数目。程序会默认将编号为0至qm_atoms - 1的原子进行量化计算。无默认值。必须设定。

#### ==qm_method==

量化计算的方法。无默认值。大小写不敏感。必须设定。

对于Gaussian，会直接将该命令的值填入对应的关键字位置，例如“b3lyp”会和基组(此处以"6-31g"为例)处理为关键字"b3lyp/6-31g"。

对于pyscf和BDF，优先比对该命令是否为"RHF"、“UHF”和"ROHF"之一，如果是将启用对应的Hatree-Folk方法。对于其他值，如"b3lyp"，则会认为使用dft方法，且对应的xc函数为该处填的值。同时，还需通过`qm_ks_method`给定对应的Kohn-Sham计算的方法。

#### qm_ks_method

Kohn-Sham计算的方法。只对pyscf和BDF有用。可选值为"RKS"、"UKS"和"ROKS"。

#### ==qm_basis==

量化计算的基组。默认值为`6-31g`。

#### qm_charge

量化计算部分的电荷量，单位为元电荷。无默认值。未设定时，将会通过分子力学中力场电荷计算得出。

#### qm_spin

量化计算部分的自旋。对于pyscf，该值为2S，默认值为`0`；对于Gaussian和BDF，该值为2S+1，默认值为`1`

#### qm_max_memory

量化部分最大使用的内存。单位为MB。默认值为`1000`。对于BDF无效。

### Gaussian特有参数

#### qm_extra_key_words

额外的关键字。将原封不动地加入至关键字行中。

#### qm_temp_file_name

临时生成的gjf的文件名。默认值为`qm_tmp_file.gjf`。

### BDF特有参数

#### qm_temp_prefix

临时生成的各文件的文件前缀。默认值为`qm_tmp_file`。

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
p(\vec r)=\sum_{f_b}p(f_b)\cdot p(\vec r|f_b)=\sum_{f_b}p(f_b)·\frac{{\rm exp}(-\beta f_b U(\vec r))}{\sum_{\vec r} {\rm exp}(-\beta f_b U(\vec r))p(\vec r){\rm exp}(\beta U(\vec r))}
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

积分的最低温度。默认值为NVT系综的T的$1/1.2$。

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

每次$f_b$蒙特卡洛模拟过程中，尝试移动的可能的最长距离。默认值为0.01。

#### sits_fcball_random_seed

$f_b$蒙特卡洛模拟过程的随机数种子。默认值为0。


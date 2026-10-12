# 高频电场模块（Maxwell）架构与工作流程

本文说明 `miniapps/electromagnetics` 里的时域全波求解器。它是这组电磁小程序中的高频模块：电场和磁场互相耦合，随时间传播。静电场在 `volta`，静磁场在 `tesla`，准静态涡流和焦耳热在 `joule`。

频域谐波形式（`curl curl E - ω²ε E = f`）不在这个小程序里，对应的示例是 `examples/ex22.cpp` 和 `examples/ex25.cpp`。

## 1. 它在库里的位置

MFEM 自下而上分成四层。Maxwell 小程序只负责物理方程、材料、边界和时间循环，离散和线性代数都交给库。

```mermaid
flowchart TB
  subgraph app ["应用层 miniapps/electromagnetics"]
    main["maxwell.cpp<br/>网格、材料、源、时间循环"]
    solver["MaxwellSolver<br/>maxwell_solver.cpp / .hpp"]
  end

  subgraph common ["小程序公共层 miniapps/common"]
    spaces["ND_ParFESpace / RT_ParFESpace"]
    curlop["ParDiscreteCurlOperator"]
  end

  subgraph fem ["有限元层 fem/"]
    forms["ParBilinearForm / ParMixedBilinearForm"]
    integ["质量、弱 curl 积分器"]
    gf["ParGridFunction"]
  end

  subgraph linalg ["线性代数与时间积分 linalg/"]
    hypre["HypreParMatrix / HyprePCG"]
    sia["SIAVSolver<br/>辛积分，1–4 阶"]
  end

  main --> solver
  solver --> spaces
  solver --> curlop
  solver --> forms
  solver --> sia
  forms --> integ
  forms --> hypre
  sia --> solver
```

并行前提是 MPI 和 hypre。`MaxwellSolver` 整段包在 `MFEM_USE_MPI` 里，单精度构建会直接跳过。

## 2. 连续方程

求解的是一阶 Maxwell 方程组，未知量是电场 `E` 和磁通密度 `B`：

```text
ε ∂E/∂t = curl(B / μ) - σ E - J
  ∂B/∂t = - curl E
```

| 符号 | 代码 | 缺省 |
| --- | --- | --- |
| `ε` | `epsCoef_` | 真空 `epsilon0_`，可用介质球覆盖 |
| `1/μ` | `muInvCoef_` | 真空 `1/mu0_`，可用磁壳覆盖 |
| `σ` | `sigmaCoef_` | 不导电；导电球内部为给定电导率 |
| `J` | `jCoef_` / `jd_` | 圆柱电流脉冲，时间包络是高斯导数 |
| `Y` | `etaInvCoef_` | 吸收边界导纳 `sqrt(ε0/μ0)` |

无耗时（没有电导、没有吸收边界）时间格式保持离散能量。有耗时，损耗项用梯形平均放到左端，能量只因物理损耗下降。

电磁能量为

```text
W = 1/2 ∫ ( ε |E|² + |B|² / μ ) dV
```

对应 `GetEnergy()`：`0.5 * (Eᵀ Mε E + Bᵀ M_{1/μ} B)`。

## 3. 离散空间

两个场放在兼容的 de Rham 序列上，多项式阶数来自命令行 `-so`。

| 场 | 空间 | 单元 | 构造 |
| --- | --- | --- | --- |
| `E`、`∂E/∂t`、`J` | H(curl) | Nédélec 棱元 | `ND_ParFESpace(mesh, order)`，内部是 `ND_FECollection(order)` |
| `B` | H(div) | Raviart–Thomas 面元 | `RT_ParFESpace(mesh, order)`，内部是 `RT_FECollection(order-1)` |

MFEM 的 RT 阶数比 ND 少 1，这样 `curl` 把 ND 自由度精确映到 RT 自由度。因此磁场更新不需要解线性方程组，只乘离散 curl。电场更新在弱形式下进行，左端是质量矩阵，必须求解。

```mermaid
flowchart LR
  E["E 在 H(curl)<br/>棱自由度"] -->|"离散 curl<br/>NegCurl_ = -C"| B["B 在 H(div)<br/>面通量自由度"]
  B -->|"弱 curl<br/>(B/μ, curl w)"| RHS["电场方程右端"]
  RHS -->|"PCG 解质量矩阵"| dEdt["∂E/∂t"]
  dEdt -->|"E += Δτ ∂E/∂t"| E
```

## 4. 类怎么分工

`maxwell.cpp` 的 `main` 是驱动，不组装矩阵。

`MaxwellSolver` 继承 `TimeDependentOperator`。对辛积分器来说它是力算子 `F`：给定 `B`，算出 `∂E/∂t`。另一个算子 `P` 是 `GetNegCurl()`，给定 `E`，算出 `∂B/∂t = -curl E`。

`SIAVSolver`（`linalg/ode.cpp`）只知道这两个算子。它交替更新 `q = B` 和 `p = E`，不关心有限元。

```text
siaSolver.Init(Maxwell.GetNegCurl(), Maxwell);          // P = -curl, F = Maxwell
siaSolver.Run(Maxwell.GetBField(), Maxwell.GetEField(), t, dt, t_snapshot);
```

自由度有两份拷贝：

- `E_`、`B_` 是 hypre 真自由度，时间积分只改它们。
- `e_`、`b_` 是带共享自由度的网格函数，给局部矩阵乘、投影和可视化用。
- 快照前调用 `SyncGridFuncs()`，把真自由度散布回网格函数。

## 5. 主程序流程

入口是 `maxwell.cpp` 的 `main`。命令行时间用纳秒，进入求解器前乘 `tScale_ = 1e-9` 变成秒。

```mermaid
flowchart TD
  start["MPI / hypre 初始化"] --> parse["解析网格、阶数、材料、边界、起止时间"]
  parse --> mesh["读串行网格"]
  mesh --> nurbs{"NURBS?"}
  nurbs -->|是| quad["均匀加密一次并设成二次曲面"]
  nurbs -->|否| sref
  quad --> sref["串行均匀加密 -rs"]
  sref --> part["按 MPI 进程划分成 ParMesh"]
  part --> pref["并行均匀加密 -rp"]
  pref --> ctor["构造 MaxwellSolver<br/>空间、材料、矩阵、CFL"]
  ctor --> ic["投影初始 E、B<br/>缺省全零"]
  ic --> snap["按 0.95 CFL 取整时间步<br/>使总步数是 10 和 5 的整齐倍数"]
  snap --> ode["SIAVSolver 阶数 -to，默认 1"]
  ode --> vis0["写出 0 时刻 VisIt / GLVis"]
  vis0 --> loop{"t < tf ?"}
  loop -->|是| run["推进到下一个快照时刻 -ts"]
  run --> energy["打印能量"]
  energy --> sync["SyncGridFuncs"]
  sync --> vis["写 VisIt，刷新 GLVis"]
  vis --> loop
  loop -->|否| done["结束"]
```

构造函数里一次做完、之后复用的东西：

1. 建 H(curl)、H(div) 并行空间。有电导或吸收边界时把算子标成 `IMPLICIT`，否则 `EXPLICIT`。
2. 材料系数：介电常数、磁导率倒数；按需再建电导率和吸收边界导纳。
3. Dirichlet 面（`-dbcs`）上取 H(curl) 真自由度 `dbc_dofs_`。边界值约束的是切向 `∂E/∂t`，缺省为 0。初始切向电场也是 0 时，这就是理想导体壁。
4. 组装 H(div) 质量矩阵 `M2MuInv_`，权是 `1/μ`。
5. 组装混合弱 curl：试验函数在 H(curl)，试探函数在 H(div)，积分 `(B/μ, curl w)`。
6. 若有耗，再组装损耗双线性型：体积分 `σ`，吸收面上边界质量 `Y`。
7. 组装离散 curl，整体乘 `-1` 得到 `NegCurl_`。
8. 电流源做成 `ParLinearForm`，每个子步按当前时间重装配。
9. 用幂迭代估计最大稳定时间步。

## 6. 一个时间步

`SIAVSolver::Step` 把一个整步拆成 `order` 个子步。系数 `a_[i]` 推进一步磁场，`b_[i]` 推进一步电场，子步之后时间增加 `a_[i] * dt`。实现依据 Candy 和 Rozmus 1991 年的可分哈密顿辛格式。

| 阶 | `a`（更新 B） | `b`（更新 E） | 结构 |
| --- | --- | --- | --- |
| 1 | `1` | `1` | 先电场，再磁场 |
| 2 | `1/2, 1/2` | `0, 1` | 半步磁场、整步电场、半步磁场 |
| 3 | `2/3, -2/3, 1` | `7/24, 3/4, -1/24` | 三阶辛格式 |
| 4 | Forest–Ruth 系数 | 同左，首项为 0 | 四阶辛格式 |

无耗时 `F` 是显式的，调用 `Mult`，它内部仍走 `implicitSolve(0, ...)`，只是左端没有损耗。有耗时调用 `ImplicitSolve(b_i * dt, B, dEdt)`，把这半个电步的步长送进左端矩阵。

```mermaid
flowchart TD
  step["SIAVSolver::Step"] --> stage["子步 i = 0 .. order-1"]
  stage --> bchk{"b_i == 0 ?"}
  bchk -->|是| skip["本子步不改 E<br/>二阶格式的第一子步就是这样"]
  bchk -->|否| kind{"有耗 ?"}
  kind -->|否| expl["Mult(B) → ∂E/∂t<br/>左端只有 ε 质量矩阵"]
  kind -->|是| impl["ImplicitSolve(b_i·dt, B)<br/>左端含 ε 和半步损耗"]
  expl --> kick["E += b_i·dt · ∂E/∂t"]
  impl --> kick
  skip --> drift
  kick --> drift["∂B/∂t = -curl E<br/>B += a_i·dt · ∂B/∂t"]
  drift --> adv["t += a_i·dt"]
  adv --> more{"还有子步?"}
  more -->|是| stage
  more -->|否| out["返回主循环"]
```

磁场那一步没有求解。`NegCurl_` 是已组装的 hypre 矩阵，直接作用在电场真自由度上。

## 7. 电场子步里解什么

`implicitSolve` 是高频模块的核心。给定当前 `B` 和时间，求 `∂E/∂t`。

```mermaid
flowchart TD
  dist["把 B 的真自由度散布到 b_"] --> weak["rhs = (B/μ, curl w)"]
  weak --> loss{"有损耗矩阵?"}
  loss -->|是| sub["rhs -= (σ E, w) + 吸收边界项"]
  loss -->|否| src
  sub --> src{"有电流源?"}
  src -->|是| jterm["按当前时间重装 jd_，rhs -= J"]
  src -->|否| bc
  jterm --> bc{"有 Dirichlet 面?"}
  bc -->|是| tan["把切向 ∂E/∂t 投影到边界"]
  bc -->|否| sys
  tan --> sys["按这个子步步长取或创建左端矩阵 A"]
  sys --> form["FormLinearSystem<br/>代入本质边界，取出真自由度"]
  form --> pcg["Jacobi 预处理 PCG<br/>容差 1e-12，最多 200 步"]
  pcg --> rec["RecoverFEMSolution<br/>写回 dedt_"]
```

无耗左端是

```text
(ε ∂E/∂t, w) = (B/μ, curl w) - (J, w)
```

有耗时，电导和吸收边界用梯形法则。记本子步电场增量为 `Δτ = b_i dt`，则

```text
(ε ∂E/∂t, w) + (Δτ/2) (σ ∂E/∂t, w) + (Δτ/2) <Y ∂E/∂t, w>_Γ
    = (B/μ, curl w) - (σ E, w) - <Y E, w>_Γ - (J, w)
```

更新 `E ← E + Δτ ∂E/∂t` 之后，损耗项等价于用了 `(E_old + E_new) / 2`。代码里 `0.5 * dt` 就是这个 `Δτ/2`，见 `setupSolver`。

高阶格式的各个 `b_i dt` 不同，左端矩阵也不同。`setupSolver` 用

```text
idt = (int)(1e6 * dt / dtMax)
```

作键，把 `A1_`、对角预处理器和 `HyprePCG` 放进 `std::map`。无耗时键恒为 0，全程只有一个质量矩阵。`idt == 0` 时不把损耗加进左端，幂迭代和能量都使用这份纯介电质量矩阵。

质量矩阵条件数好，所以预处理器只是 hypre 对角缩放，没有 AMG。

三种边界的离散含义：

| 类型 | 命令行 | 离散处理 |
| --- | --- | --- |
| 自然边界 | 未列入 `-dbcs`、`-abcs` 的面 | 弱 curl 不产生边界积分，切向电流为 0 |
| 理想导体 | `-dbcs`，缺省 `dEdtBCFunc = 0` | 约束切向 `∂E/∂t`。初始切向 `E` 为 0 时等效 `n × E = 0` |
| 一阶吸收边界 | `-abcs` | Sommerfeld 条件，边界质量权为自由空间导纳 |

## 8. 时间步长

`GetMaximumTimeStep` 对下面这个算子做 20 次幂迭代，或直到相对变化小于 `0.001`：

```text
v  →  Mε⁻¹  Cᵀ  M_{1/μ}  C  v
```

最大特征值 `λ` 对应最高角频率的平方，稳定极限是 `dtMax = 2 / sqrt(λ)`。主程序再乘安全因子 `-sf`（默认 0.95），并用 `SnapTimeStep` 把步长稍微改小，使 `(tf - ti) / dt` 落在 1、2、5 乘 10 的幂这种步数上。

## 9. 一次运行里数据怎么走

以示例命令里的电流脉冲、介质球和吸收球为例：

```text
mpirun -np 4 maxwell -m ../../data/ball-nurbs.mesh -rs 2 \
  -abcs '-1' \
  -dp '-0.3 0.0 0.0 0.3 0.0 0.0 0.1 1 .5 .5'
```

```mermaid
sequenceDiagram
  participant Main as maxwell.cpp
  participant M as MaxwellSolver
  participant S as SIAVSolver
  participant P as PCG

  Main->>M: 构造（空间、矩阵、CFL）
  Main->>M: 投影 E=0, B=0
  Main->>S: Init(-curl, Maxwell)
  loop 直到 tf
    S->>M: 给定 B，求 ∂E/∂t
    M->>M: 弱 curl、减去 J 和损耗
    M->>P: A ∂E/∂t = rhs
    P-->>M: ∂E/∂t
    S->>S: E += Δτ ∂E/∂t
    S->>S: B += Δτ (-curl E)
    Main->>M: 能量、SyncGridFuncs、VisIt/GLVis
  end
```

电流只在圆柱棒内非零，方向沿棒轴。时间因子是

```text
a (t - b) exp( -0.5 ((t-b)/c)² ) / c²
```

`a、b、c` 在输入时是纳秒，函数内部乘 `tScale_`。这个脉冲的频谱决定了它是一个宽带高频源；网格加密和时间步必须跟得上脉冲宽度。

## 10. 建议的阅读顺序

1. `maxwell.cpp` 文件头的方程注释，再读 `main` 的时间循环。
2. `maxwell_solver.hpp` 的成员，对照第 4 节的矩阵。
3. `MaxwellSolver::MaxwellSolver` 的组装顺序。
4. `implicitSolve` 和 `setupSolver`，这是电场更新。
5. `linalg/ode.cpp` 里 `SIAVSolver::Step`，对照第 6 节的子步。
6. `GetMaximumTimeStep` 和 `GetEnergy`，确认 CFL 与能量用的是同一套矩阵。

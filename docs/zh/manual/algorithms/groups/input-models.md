# 数据加载与输入模型（Data Loading and Input Models）

<a href="../../../../en/index.html">English</a> · **简体中文**

这一门类覆盖算法的输入前提：把经典数据、布尔真值表与算子装进 Oracle、块编码和量子数据结构。多数条目来自 `oracq.algorithms.input_model` 子包（`prepare-select` 来自 `common`）；在使用任何估计、求解或模拟类算法之前，通常需要先在这里选择并构造合适的访问模型。

```{toctree}
:maxdepth: 1

../state-preparation
../alias-preparation
../xor-database
../adjacency-oracle
../select-swap
../qrom-lookup
../block-encoding-algebra
../sparse-access
../sparse-block-encoding
../prepare-select
../purification
../gibbs-state
../double-factorization
../thc
```

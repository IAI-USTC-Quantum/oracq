# 运行与修改算法展示目录

展示目录把不同类别的算法放在统一的运行流程中，便于比较它们的输入、寄存器和读出方式。

```bash
uv run python examples/algorithm_gallery.py
```

每个子目录包含 `closed.rir.json`、`modular.originir` 和 `toffoli_u3_cz.originir`。Bernstein–Vazirani 还保存绑定前的开放描述。`index.json` 记录类别、模块数量和读出说明。

## 选择要检查的结果

不要对所有算法使用同一种读出方式。BV 读取秘密字符串；QPE 和振幅估计读取 phase；Hadamard test 估计 probe 的 Z 期望；重复码恢复则检查逻辑态，并保留 syndrome。

用以下代码查找案例：

```{testcode}
from pyqecclang.applications.gallery import algorithm_gallery

cases = {case.name: case for case in algorithm_gallery()}
assert "qaoa_maxcut" in cases
assert cases["order_finding"].family == "number_theory"
assert len(cases) == 22
```

需要添加应用例子时，可以在自己的脚本里调用对应算法文件。若要扩展公开展示目录，则为 `GalleryCase` 提供操作和具体读出说明，并增加独立的数学见证。

## 真实后端对拍

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

脚本逐项比较参考执行器、PySparQ 和 OriginIR 的复幅度，超出容差立即退出。它验证同一线路在各后端中的语义一致；某个参数是否适合真实问题，仍需要应用层分析。

# 命令行

命令行处理已经保存的 RIR，也提供数学函数编译入口。

```bash
pyqecclang validate program.rir.json
pyqecclang inspect program.rir.json
pyqecclang requirements program.rir.json
pyqecclang emit program.rir.json -o program.originir
pyqecclang emit program.rir.json --basis toffoli-u3-cz -o basis.originir
```

`validate` 检查结构并报告开放状态。`inspect` 输出入口的寄存器、资源和能力描述。`requirements` 列出入口可达的未实现 oracle 及调用路径。

## 绑定实现

```bash
pyqecclang bind open.rir.json --bindings bindings.json --report binding-report.json -o closed.rir.json
```

绑定清单将槽名映射到实现的 RIR 文件与资源名称：

```json
{
  "Function": {
    "program": "lookup.rir.json",
    "resources": {"table": "values"}
  }
}
```

`program` 相对于绑定清单所在目录解析。实现的布局和 alpha 必须与槽位一致。

## 开放与闭合资源分析

```bash
pyqecclang estimate open.rir.json --allow-open -o open-cost.json
pyqecclang estimate closed.rir.json -o closed-cost.json
```

开放报告保留 oracle 调用数与未知工作区；已知成本不能当作最终总成本。
绑定失败时仍会写出 `--report` 指定的诊断，并以非零状态退出。

## 执行

```bash
pyqecclang run closed.rir.json --memory memory.json
pyqecclang run closed.rir.json --memory memory.json --backend pysparq
pyqecclang run closed.rir.json --memory memory.json --backend originir
```

默认使用参考执行器。内存文件按入口资源名提供列表或地址字典，缺省单元为零。执行失败会返回非零退出码，不会静默省略未完成模块。

## 数学函数与 QHAM

```bash
pyqecclang compile-function examples/math_functions.py --function pressure \
  --width 12 --fraction 6 --mir-output out/pressure.mir.json -o out/pressure.rir.json

python -m pyqecclang.applications.qham --example burgers --order 2 --eta=-0.4
```

QHAM 的旧入口 `python -m pyqecclang.qham` 保留兼容。它和新入口调用同一实现。

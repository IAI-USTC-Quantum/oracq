# 命令行

命令行处理已经保存的 RIR，也提供[数学函数编译](math-functions.md#命令行与案例)入口。RIR 输入自动识别 YAML 与 JSON 两种文本（约定扩展名为 `.rir.yaml`）；写出 RIR 文本时以 `--format yaml|json` 选择，默认 `yaml`。

```bash
oracq validate program.rir.yaml
oracq inspect program.rir.yaml
oracq requirements program.rir.yaml
oracq emit program.rir.yaml -o program.originir
oracq emit program.rir.yaml --basis toffoli-u3-cz -o basis.originir
```

{obj}`validate <oracq.infrastructure.validation.validate>` 检查结构并报告开放状态。`inspect` 输出入口的寄存器、资源和能力描述。`requirements` 列出入口可达的未实现 oracle 及调用路径。

## 绑定实现

```bash
oracq bind open.rir.yaml --bindings bindings.json --report binding-report.json -o closed.rir.yaml
```

绑定清单是 JSON 文件，将槽名映射到实现的 RIR 文件与资源名称：

```json
{
  "Function": {
    "program": "lookup.rir.yaml",
    "resources": {"table": "values"}
  }
}
```

`program` 相对于绑定清单所在目录解析，指向的 RIR 文件可为 YAML 或 JSON 文本。实现的布局和 alpha 必须与槽位一致。`--format` 作用于 {obj}`bind <oracq.infrastructure.linking.bind>`、`canonicalize` 与 `compile-function` 写出的 RIR 文本；两种格式解析得到的程序逐字段一致。绑定语义与可运行示例见[教程：替换 oracle](../tutorials/oracle-binding.md)。

## 开放与闭合资源分析

```bash
oracq estimate open.rir.yaml --allow-open -o open-cost.json
oracq estimate closed.rir.yaml -o closed-cost.json
```

开放报告保留 oracle 调用数与未知工作区；已知成本不能当作最终总成本。
计数口径与已知差距见[资源估计](resource-estimation.md)。
绑定失败时仍会写出 `--report` 指定的诊断，并以非零状态退出。

## 执行

```bash
oracq run closed.rir.yaml --memory memory.qram.yaml
oracq run closed.rir.yaml --memory memory.qram.yaml --backend pysparq
oracq run closed.rir.yaml --memory memory.qram.yaml --backend originir
```

默认使用[参考执行器](backends.md#参考执行器)。内存文件为 qram YAML（格式见规范参考的 [QRAM 内存定义](../reference/qram-memory.md)），按入口资源名列出段，缺省单元为零。执行失败会返回非零退出码，不会静默省略未完成模块。

## 数学函数与 QHAM

```bash
oracq compile-function examples/math_functions.py --function pressure \
  --width 12 --fraction 6 --mir-output out/pressure.mir.json -o out/pressure.rir.yaml

python -m oracq.applications.qham --example burgers --order 2 --eta=-0.4
```

QHAM 的旧入口 `python -m oracq.qham` 保留兼容。它和新入口调用同一实现，完整用法见[一般 QHAM 自动生成](qham.md)。

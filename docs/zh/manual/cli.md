# 命令行

<a href="../../manual/cli.html">English</a> · **简体中文**

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

### 独立编译工作台

```bash
oracq gui
oracq gui --no-browser --port 8000
```

安装包自带独立本地界面，不需要另取 demo-suite，也不访问外部网页资源。
编辑有界的纯函数，选择定点位宽、小数位、数学核阶数和检查点大小，即可将真实
编译结果与保留全部中间量的基线比较。两套程序使用完全相同的数值配置；检查点
可能减少空间并增加门数，界面同时报告两者。

峰值比特数来自实际严格门集后端布局，包含额外工作比特。可以下载标准压缩
RIR JSON、严格模块化门表、数学 IR 和来源报告，供下游软件消费。这些是旋转
综合、QRAM 时间建模和容错调度之前的逻辑成本；仅完成编译不等于认证数值精度。

每次编译使用独立进程，限时60秒。编辑器支持不含循环、至多128个 AST 节点、
位宽2..16、数学核阶数1..8的示例；下列 CLI 支持更大的源程序及静态循环。
源码由受限编译前端解释，不作为 Python 脚本执行；服务只绑定本机回环地址。

### 函数及 PDE 命令行

```bash
oracq compile-function examples/math_functions.py --function pressure \
  --width 12 --fraction 6 --checkpoint-size 16 \
  --mir-output out/pressure.mir.json -o out/pressure.rir.yaml

python -m oracq.applications.qham --example burgers --order 2 --eta=-0.4
```

QHAM 的旧入口 `python -m oracq.qham` 保留兼容。它和新入口调用同一实现，完整用法见[一般 QHAM 自动生成](qham.md)。

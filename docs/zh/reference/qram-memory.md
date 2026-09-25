# QRAM 内存定义（qram YAML）

<a href="../../en/reference/qram-memory.html">English</a> · **简体中文**

内存数据不写入程序序列化文本（见 [RIR](rir.md)），而是作为执行输入另行绑定。本文定义绑定文件 `*.qram.yaml` 的格式：一个文件描述一组 QRAM 段（`qram_segments` 列表），执行时按入口资源名与程序交叉校验。

```yaml
qram_segments:
  - name: values            # 存储名字，代码中按此索引
    address_length: 2       # 地址位宽
    word_length: 3          # 数据字位宽
    type: uint              # 数据类型：uint / sint / fixedpoint
    data: [1, 2, 4, 7]      # 具体数据，稠密数组，下标即地址
```

## 字段

| 字段 | 类型 | 约束 | 说明 |
|------|------|------|------|
| `qram_segments` | list | 必填，顶层唯一键 | 一组 QRAM 定义 |
| `name` | string | 非空，文件内唯一 | 存储名字，执行时与入口资源名精确匹配 |
| `address_length` | integer | 1..64 | 地址位宽，对应 RIR QRAM 的 `address_width` |
| `word_length` | integer | 1..64 | 数据字位宽，对应 RIR QRAM 的 `data_width` |
| `type` | string | `uint` / `sint` / `fixedpoint` | data 的源数据类型，加载时按其编码为无符号字 |
| `data` | 数组 | 见下表 | 稠密源数据数组，下标即地址 |

三种数据类型的取值范围与编码：

| type | data 元素范围 | 加载时编码为无符号字 |
|------|--------------|---------------------|
| `uint` | 无符号整数，0 .. 2^`word_length` − 1 | 原样 |
| `sint` | 有符号整数，−2^(`word_length`−1) .. 2^(`word_length`−1) − 1 | 二补码位模式 |
| `fixedpoint` | 定点数，[0, 1) | 乘以 2^`word_length` 后向零截断（与 `FixedFormat.encode` 同口径） |

## 语义规则

- 加载期按 `word_length` 本地校验每个源数据（范围见上表），越界即报错；`data` 长度不能超过 2^`address_length`，短于它时高位单元按零补齐（padding）。
- 段名在文件内必须唯一。
- 执行期以程序内的 QRAM 声明为准：{obj}`check_memory <oracq.infrastructure.execution.check_memory>` 要求段的资源名集合与入口声明完全一致、不多不少；编码后的字必须落在程序声明位宽内，`address_length`/`word_length` 与程序声明不一致不单独报错。
- {obj}`dump_qram_yaml <oracq.infrastructure.qram_schema.dump_qram_yaml>` 产出的段类型恒为 `uint`（执行器侧的字已是位宽内无符号整数）。

## Python API

```python
from oracq import dump_qram_yaml, load_qram_yaml, simulate

memory = load_qram_yaml("memory.qram.yaml")  # -> {"values": [1, 2, 4, 7]}
state = simulate(program, memory)             # 与既有 memory 参数完全兼容

text = dump_qram_yaml(program, memory)        # 规整校验后输出文本，由调用方写盘
```

{obj}`load_qram_yaml <oracq.infrastructure.qram_schema.load_qram_yaml>` 返回资源名到稠密字数组的映射（已按 type 编码并补零到全长），与 {obj}`simulate <oracq.infrastructure.execution.simulate>`、{obj}`run_pysparq <oracq.infrastructure.backends.pysparq.run_pysparq>`、{obj}`run_originir <oracq.infrastructure.backends.originir.run_originir>` 的 `memory` 参数同形。`dump_qram_yaml` 先经 `check_memory` 交叉校验并把稀疏字典稠密化，再按入口资源声明顺序生成段；入口无 QRAM 资源时输出空 `qram_segments`。

## CLI

```bash
oracq run closed.rir.yaml --memory memory.qram.yaml
```

`--memory` 只接受本格式。结构或数据不合法时以退出码 2 报中文错误，不会静默截断。全部子命令见[命令行](../manual/cli.md)。

## Schema

结构约束见 [qram-memory.schema.json](schemas/qram-memory.schema.json)，含按 `type` 分支的 `data` 元素类型约束。跨字段约束（数组长度对 `address_length`、字值范围对 `word_length` 与 `type`、段名唯一）无法用 JSON Schema 表达，由 Python 加载器强制。

pyqecclang 包总览
============================

根包 ``pyqecclang`` 汇总导出公开 API（共 92 个名字）。
名字按定义模块分组；模块标题链接到对应 API 页，成员链接到模块页内的完整说明。

基础设施 API
----------------

:doc:`Toffoli / U3 / CZ 降低 <infrastructure/backends/basis>`
    ``pyqecclang.infrastructure.backends.basis`` —— :obj:`export_toffoli_u3_cz <pyqecclang.infrastructure.backends.basis.export_toffoli_u3_cz>`

:doc:`OriginIR-ext 后端 <infrastructure/backends/originir>`
    ``pyqecclang.infrastructure.backends.originir`` —— :obj:`OriginIRArtifact <pyqecclang.infrastructure.backends.originir.OriginIRArtifact>`、:obj:`export_originir <pyqecclang.infrastructure.backends.originir.export_originir>`、:obj:`run_originir <pyqecclang.infrastructure.backends.originir.run_originir>`

:doc:`PySparQ 后端 <infrastructure/backends/pysparq>`
    ``pyqecclang.infrastructure.backends.pysparq`` —— :obj:`run_pysparq <pyqecclang.infrastructure.backends.pysparq.run_pysparq>`、:obj:`run_pysparq_rir <pyqecclang.infrastructure.backends.pysparq.run_pysparq_rir>`

:doc:`Quantikz 线路导出 <infrastructure/backends/quantikz>`
    ``pyqecclang.infrastructure.backends.quantikz`` —— :obj:`quantikz <pyqecclang.infrastructure.backends.quantikz.quantikz>`

:doc:`严格网表导出 <infrastructure/backends/strict>`
    ``pyqecclang.infrastructure.backends.strict`` —— :obj:`StrictArtifact <pyqecclang.infrastructure.backends.strict.StrictArtifact>`、:obj:`export_strict <pyqecclang.infrastructure.backends.strict.export_strict>`

:doc:`模块构造器 <infrastructure/builder>`
    ``pyqecclang.infrastructure.builder`` —— :obj:`Builder <pyqecclang.infrastructure.builder.Builder>`、:obj:`Operation <pyqecclang.infrastructure.builder.Operation>`

:doc:`资源估计 <infrastructure/estimate>`
    ``pyqecclang.infrastructure.estimate`` —— :obj:`ResourceEstimate <pyqecclang.infrastructure.estimate.ResourceEstimate>`、:obj:`estimate_resources <pyqecclang.infrastructure.estimate.estimate_resources>`

:doc:`寄存器参考执行器 <infrastructure/execution>`
    ``pyqecclang.infrastructure.execution`` —— :obj:`RegisterState <pyqecclang.infrastructure.execution.RegisterState>`、:obj:`simulate <pyqecclang.infrastructure.execution.simulate>`

:doc:`RIR 对象 <infrastructure/ir>`
    ``pyqecclang.infrastructure.ir`` —— :obj:`VERSION <pyqecclang.infrastructure.ir.VERSION>`、:obj:`Adjoint <pyqecclang.infrastructure.ir.Adjoint>`、:obj:`Bits <pyqecclang.infrastructure.ir.Bits>`、:obj:`Call <pyqecclang.infrastructure.ir.Call>`、:obj:`Control <pyqecclang.infrastructure.ir.Control>`、:obj:`Load <pyqecclang.infrastructure.ir.Load>`、:obj:`Module <pyqecclang.infrastructure.ir.Module>`、:obj:`Primitive <pyqecclang.infrastructure.ir.Primitive>`、:obj:`Program <pyqecclang.infrastructure.ir.Program>`、:obj:`QRAM <pyqecclang.infrastructure.ir.QRAM>`、:obj:`Rational <pyqecclang.infrastructure.ir.Rational>`、:obj:`Ref <pyqecclang.infrastructure.ir.Ref>`、:obj:`Register <pyqecclang.infrastructure.ir.Register>`、:obj:`RegType <pyqecclang.infrastructure.ir.RegType>`、:obj:`Repeat <pyqecclang.infrastructure.ir.Repeat>`、:obj:`Resource <pyqecclang.infrastructure.ir.Resource>`、:obj:`SInt <pyqecclang.infrastructure.ir.SInt>`、:obj:`Span <pyqecclang.infrastructure.ir.Span>`、:obj:`Store <pyqecclang.infrastructure.ir.Store>`、:obj:`UInt <pyqecclang.infrastructure.ir.UInt>`、:obj:`ValidationError <pyqecclang.infrastructure.ir.ValidationError>`、:obj:`fuse <pyqecclang.infrastructure.ir.fuse>`

:doc:`绑定与能力分析 <infrastructure/linking>`
    ``pyqecclang.infrastructure.linking`` —— :obj:`Binding <pyqecclang.infrastructure.linking.Binding>`、:obj:`OracleRequirement <pyqecclang.infrastructure.linking.OracleRequirement>`、:obj:`bind <pyqecclang.infrastructure.linking.bind>`、:obj:`capabilities <pyqecclang.infrastructure.linking.capabilities>`、:obj:`unresolved <pyqecclang.infrastructure.linking.unresolved>`

:doc:`数学函数编译入口 <infrastructure/mathfunc>`
    ``pyqecclang.infrastructure.mathfunc`` —— :obj:`compile_function <pyqecclang.infrastructure.mathfunc.compile_function>`、:obj:`lower_math_ir <pyqecclang.infrastructure.mathfunc.lower_math_ir>`

:doc:`Python 数学函数前端 <infrastructure/mathfunc/frontend>`
    ``pyqecclang.infrastructure.mathfunc.frontend`` —— :obj:`FunctionCompileError <pyqecclang.infrastructure.mathfunc.frontend.FunctionCompileError>`

:doc:`MIR 对象 <infrastructure/mathfunc/graph>`
    ``pyqecclang.infrastructure.mathfunc.graph`` —— :obj:`Index <pyqecclang.infrastructure.mathfunc.graph.Index>`、:obj:`MathProgram <pyqecclang.infrastructure.mathfunc.graph.MathProgram>`

:doc:`数学函数降低 <infrastructure/mathfunc/lowering>`
    ``pyqecclang.infrastructure.mathfunc.lowering`` —— :obj:`CompiledFunction <pyqecclang.infrastructure.mathfunc.lowering.CompiledFunction>`

:doc:`数学核生成 <infrastructure/mathfunc/numeric>`
    ``pyqecclang.infrastructure.mathfunc.numeric`` —— :obj:`MathConfig <pyqecclang.infrastructure.mathfunc.numeric.MathConfig>`

:doc:`原生实现注册 <infrastructure/native>`
    ``pyqecclang.infrastructure.native`` —— :obj:`DynamicCppFactory <pyqecclang.infrastructure.native.DynamicCppFactory>`、:obj:`NativeRegistry <pyqecclang.infrastructure.native.NativeRegistry>`

:doc:`QRAM 指针式读写 <infrastructure/qmem>`
    ``pyqecclang.infrastructure.qmem`` —— :obj:`QMem <pyqecclang.infrastructure.qmem.QMem>`、:obj:`QPtr <pyqecclang.infrastructure.qmem.QPtr>`

:doc:`RIR 序列化 <infrastructure/serialization>`
    ``pyqecclang.infrastructure.serialization`` —— :obj:`dumps <pyqecclang.infrastructure.serialization.dumps>`、:obj:`loads <pyqecclang.infrastructure.serialization.loads>`

:doc:`结构验证 <infrastructure/validation>`
    ``pyqecclang.infrastructure.validation`` —— :obj:`validate <pyqecclang.infrastructure.validation.validate>`

算法 API
------------

:doc:`可逆算术 <algorithms/arithmetic>`
    ``pyqecclang.algorithms.arithmetic`` —— :obj:`FixedFormat <pyqecclang.algorithms.arithmetic.FixedFormat>`、:obj:`arithmetic_native_registry <pyqecclang.algorithms.arithmetic.arithmetic_native_registry>`、:obj:`fixed_arithmetic <pyqecclang.algorithms.arithmetic.fixed_arithmetic>`

:doc:`算法契约与报告 <algorithms/contracts>`
    ``pyqecclang.algorithms.contracts`` —— :obj:`ContractError <pyqecclang.algorithms.contracts.ContractError>`、:obj:`ContractIssue <pyqecclang.algorithms.contracts.ContractIssue>`、:obj:`ContractReport <pyqecclang.algorithms.contracts.ContractReport>`、:obj:`InputRequirement <pyqecclang.algorithms.contracts.InputRequirement>`、:obj:`OracleCapabilities <pyqecclang.algorithms.contracts.OracleCapabilities>`、:obj:`OracleSpec <pyqecclang.algorithms.contracts.OracleSpec>`、:obj:`ProtocolContract <pyqecclang.algorithms.contracts.ProtocolContract>`、:obj:`describe_oracle <pyqecclang.algorithms.contracts.describe_oracle>`、:obj:`requires <pyqecclang.algorithms.contracts.requires>`

:doc:`QODE 组装接口 <algorithms/ode>`
    ``pyqecclang.algorithms.ode`` —— :obj:`QODEProblem <pyqecclang.algorithms.ode.QODEProblem>`、:obj:`QODEProtocol <pyqecclang.algorithms.ode.QODEProtocol>`

:doc:`算子包装与基本组合 <algorithms/operators>`
    ``pyqecclang.algorithms.operators`` —— :obj:`BlockEncoding <pyqecclang.algorithms.operators.BlockEncoding>`、:obj:`Generator <pyqecclang.algorithms.operators.Generator>`、:obj:`block_encoding <pyqecclang.algorithms.operators.block_encoding>`、:obj:`identity <pyqecclang.algorithms.operators.identity>`、:obj:`linear_combination <pyqecclang.algorithms.operators.linear_combination>`、:obj:`pauli_x <pyqecclang.algorithms.operators.pauli_x>`、:obj:`product <pyqecclang.algorithms.operators.product>`、:obj:`scale <pyqecclang.algorithms.operators.scale>`、:obj:`zero <pyqecclang.algorithms.operators.zero>`

:doc:`Oracle 声明与实现 <algorithms/oracles>`
    ``pyqecclang.algorithms.oracles`` —— :obj:`declare <pyqecclang.algorithms.oracles.declare>`

:doc:`量子数据结构（qsample 与 sample-and-query） <algorithms/qdata>`
    ``pyqecclang.algorithms.qdata`` —— :obj:`QMatrix <pyqecclang.algorithms.qdata.QMatrix>`、:obj:`QVector <pyqecclang.algorithms.qdata.QVector>`

:doc:`量子线性系统 <algorithms/qlss>`
    ``pyqecclang.algorithms.qlss`` —— :obj:`BlockSystem <pyqecclang.algorithms.qlss.BlockSystem>`、:obj:`LinearSystem <pyqecclang.algorithms.qlss.LinearSystem>`、:obj:`QLSSProtocol <pyqecclang.algorithms.qlss.QLSSProtocol>`、:obj:`SolveResult <pyqecclang.algorithms.qlss.SolveResult>`、:obj:`SparseSystem <pyqecclang.algorithms.qlss.SparseSystem>`、:obj:`SpectralPromise <pyqecclang.algorithms.qlss.SpectralPromise>`

:doc:`KP 量子推荐系统 <algorithms/recommendation>`
    ``pyqecclang.algorithms.recommendation`` —— :obj:`KPRecommendationConfig <pyqecclang.algorithms.recommendation.KPRecommendationConfig>`、:obj:`RecommendationResult <pyqecclang.algorithms.recommendation.RecommendationResult>`、:obj:`kp_recommendation <pyqecclang.algorithms.recommendation.kp_recommendation>`、:obj:`sigma_from_phase <pyqecclang.algorithms.recommendation.sigma_from_phase>`

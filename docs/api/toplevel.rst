oracq package overview
============================================

The root package ``oracq`` re-exports the public API (104 names).
Names are grouped by their defining module; module titles link to the corresponding API page and members link to their full descriptions there.

Infrastructure API
------------------------------------

:doc:`Toffoli / U3 / CZ lowering <infrastructure/backends/basis>`
    ``oracq.infrastructure.backends.basis`` — :obj:`export_toffoli_u3_cz <oracq.infrastructure.backends.basis.export_toffoli_u3_cz>`

:doc:`OriginIR-ext backend <infrastructure/backends/originir>`
    ``oracq.infrastructure.backends.originir`` — :obj:`OriginIRArtifact <oracq.infrastructure.backends.originir.OriginIRArtifact>`, :obj:`export_originir <oracq.infrastructure.backends.originir.export_originir>`, :obj:`run_originir <oracq.infrastructure.backends.originir.run_originir>`

:doc:`PySparQ backend <infrastructure/backends/pysparq>`
    ``oracq.infrastructure.backends.pysparq`` — :obj:`run_pysparq <oracq.infrastructure.backends.pysparq.run_pysparq>`, :obj:`run_pysparq_rir <oracq.infrastructure.backends.pysparq.run_pysparq_rir>`

:doc:`Quantikz circuit export <infrastructure/backends/quantikz>`
    ``oracq.infrastructure.backends.quantikz`` — :obj:`quantikz <oracq.infrastructure.backends.quantikz.quantikz>`

:doc:`Strict netlist export <infrastructure/backends/strict>`
    ``oracq.infrastructure.backends.strict`` — :obj:`StrictArtifact <oracq.infrastructure.backends.strict.StrictArtifact>`, :obj:`export_strict <oracq.infrastructure.backends.strict.export_strict>`

:doc:`Module builder <infrastructure/builder>`
    ``oracq.infrastructure.builder`` — :obj:`Builder <oracq.infrastructure.builder.Builder>`, :obj:`Operation <oracq.infrastructure.builder.Operation>`

:doc:`Resource estimation <infrastructure/estimate>`
    ``oracq.infrastructure.estimate`` — :obj:`ResourceEstimate <oracq.infrastructure.estimate.ResourceEstimate>`, :obj:`estimate_resources <oracq.infrastructure.estimate.estimate_resources>`, :obj:`OracleCall <oracq.infrastructure.estimate.OracleCall>`, :obj:`RotationCounts <oracq.infrastructure.estimate.RotationCounts>`

:doc:`Register reference executor <infrastructure/execution>`
    ``oracq.infrastructure.execution`` — :obj:`RegisterState <oracq.infrastructure.execution.RegisterState>`, :obj:`simulate <oracq.infrastructure.execution.simulate>`

:doc:`RIR objects <infrastructure/ir>`
    ``oracq.infrastructure.ir`` — :obj:`VERSION <oracq.infrastructure.ir.VERSION>`, :obj:`Adjoint <oracq.infrastructure.ir.Adjoint>`, :obj:`Bits <oracq.infrastructure.ir.Bits>`, :obj:`Call <oracq.infrastructure.ir.Call>`, :obj:`Control <oracq.infrastructure.ir.Control>`, :obj:`Load <oracq.infrastructure.ir.Load>`, :obj:`Module <oracq.infrastructure.ir.Module>`, :obj:`Primitive <oracq.infrastructure.ir.Primitive>`, :obj:`Program <oracq.infrastructure.ir.Program>`, :obj:`QRAM <oracq.infrastructure.ir.QRAM>`, :obj:`Rational <oracq.infrastructure.ir.Rational>`, :obj:`Ref <oracq.infrastructure.ir.Ref>`, :obj:`Register <oracq.infrastructure.ir.Register>`, :obj:`RegType <oracq.infrastructure.ir.RegType>`, :obj:`Repeat <oracq.infrastructure.ir.Repeat>`, :obj:`Resource <oracq.infrastructure.ir.Resource>`, :obj:`SInt <oracq.infrastructure.ir.SInt>`, :obj:`Span <oracq.infrastructure.ir.Span>`, :obj:`Store <oracq.infrastructure.ir.Store>`, :obj:`UInt <oracq.infrastructure.ir.UInt>`, :obj:`ValidationError <oracq.infrastructure.ir.ValidationError>`, :obj:`fuse <oracq.infrastructure.ir.fuse>`

:doc:`Binding and capability analysis <infrastructure/linking>`
    ``oracq.infrastructure.linking`` — :obj:`Binding <oracq.infrastructure.linking.Binding>`, :obj:`OracleRequirement <oracq.infrastructure.linking.OracleRequirement>`, :obj:`bind <oracq.infrastructure.linking.bind>`, :obj:`capabilities <oracq.infrastructure.linking.capabilities>`, :obj:`unresolved <oracq.infrastructure.linking.unresolved>`, :obj:`BindingError <oracq.infrastructure.linking.BindingError>`, :obj:`BindingReport <oracq.infrastructure.linking.BindingReport>`, :obj:`BindingResult <oracq.infrastructure.linking.BindingResult>`, :obj:`bind_with_report <oracq.infrastructure.linking.bind_with_report>`

:doc:`Math-function compilation entry <infrastructure/mathfunc>`
    ``oracq.infrastructure.mathfunc`` — :obj:`compile_function <oracq.infrastructure.mathfunc.compile_function>`, :obj:`lower_math_ir <oracq.infrastructure.mathfunc.lower_math_ir>`

:doc:`Python math-function frontend <infrastructure/mathfunc/frontend>`
    ``oracq.infrastructure.mathfunc.frontend`` — :obj:`FunctionCompileError <oracq.infrastructure.mathfunc.frontend.FunctionCompileError>`

:doc:`MIR objects <infrastructure/mathfunc/graph>`
    ``oracq.infrastructure.mathfunc.graph`` — :obj:`Index <oracq.infrastructure.mathfunc.graph.Index>`, :obj:`MathProgram <oracq.infrastructure.mathfunc.graph.MathProgram>`

:doc:`Math-function lowering <infrastructure/mathfunc/lowering>`
    ``oracq.infrastructure.mathfunc.lowering`` — :obj:`CompiledFunction <oracq.infrastructure.mathfunc.lowering.CompiledFunction>`

:doc:`Math kernel generation <infrastructure/mathfunc/numeric>`
    ``oracq.infrastructure.mathfunc.numeric`` — :obj:`MathConfig <oracq.infrastructure.mathfunc.numeric.MathConfig>`

:doc:`Native implementation registry <infrastructure/native>`
    ``oracq.infrastructure.native`` — :obj:`DynamicCppFactory <oracq.infrastructure.native.DynamicCppFactory>`, :obj:`NativeRegistry <oracq.infrastructure.native.NativeRegistry>`

:doc:`QRAM pointer-style read/write <infrastructure/qmem>`
    ``oracq.infrastructure.qmem`` — :obj:`QMem <oracq.infrastructure.qmem.QMem>`, :obj:`QPtr <oracq.infrastructure.qmem.QPtr>`

:doc:`QRAM YAML memory format <infrastructure/qram_schema>`
    ``oracq.infrastructure.qram_schema`` — :obj:`dump_qram_yaml <oracq.infrastructure.qram_schema.dump_qram_yaml>`, :obj:`load_qram_yaml <oracq.infrastructure.qram_schema.load_qram_yaml>`

:doc:`RIR serialization <infrastructure/serialization>`
    ``oracq.infrastructure.serialization`` — :obj:`dumps <oracq.infrastructure.serialization.dumps>`, :obj:`loads <oracq.infrastructure.serialization.loads>`

:doc:`Structural validation <infrastructure/validation>`
    ``oracq.infrastructure.validation`` — :obj:`validate <oracq.infrastructure.validation.validate>`

Algorithms API
----------------------------

:doc:`Reversible arithmetic <algorithms/common/arithmetic>`
    ``oracq.algorithms.common.arithmetic`` — :obj:`FixedFormat <oracq.algorithms.common.arithmetic.FixedFormat>`, :obj:`arithmetic_native_registry <oracq.algorithms.common.arithmetic.arithmetic_native_registry>`, :obj:`fixed_arithmetic <oracq.algorithms.common.arithmetic.fixed_arithmetic>`

:doc:`Algorithm contracts and reports <algorithms/input_model/contracts>`
    ``oracq.algorithms.input_model.contracts`` — :obj:`AlgorithmContract <oracq.algorithms.input_model.contracts.AlgorithmContract>`, :obj:`ResolvedInputs <oracq.algorithms.input_model.contracts.ResolvedInputs>`, :obj:`ContractError <oracq.algorithms.input_model.contracts.ContractError>`, :obj:`ContractIssue <oracq.algorithms.input_model.contracts.ContractIssue>`, :obj:`ContractReport <oracq.algorithms.input_model.contracts.ContractReport>`, :obj:`InputRequirement <oracq.algorithms.input_model.contracts.InputRequirement>`, :obj:`OracleCapabilities <oracq.algorithms.input_model.contracts.OracleCapabilities>`, :obj:`OracleSpec <oracq.algorithms.input_model.contracts.OracleSpec>`, :obj:`ProtocolContract <oracq.algorithms.input_model.contracts.ProtocolContract>`, :obj:`describe_oracle <oracq.algorithms.input_model.contracts.describe_oracle>`, :obj:`requires <oracq.algorithms.input_model.contracts.requires>`

:doc:`Operator wrappers and basic composition <algorithms/input_model/operators>`
    ``oracq.algorithms.input_model.operators`` — :obj:`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`, :obj:`Generator <oracq.algorithms.input_model.operators.Generator>`, :obj:`block_encoding <oracq.algorithms.input_model.operators.block_encoding>`, :obj:`identity <oracq.algorithms.input_model.operators.identity>`, :obj:`linear_combination <oracq.algorithms.input_model.operators.linear_combination>`, :obj:`pauli_x <oracq.algorithms.input_model.operators.pauli_x>`, :obj:`product <oracq.algorithms.input_model.operators.product>`, :obj:`scale <oracq.algorithms.input_model.operators.scale>`, :obj:`zero <oracq.algorithms.input_model.operators.zero>`

:doc:`Oracle declarations and implementations <algorithms/input_model/oracles>`
    ``oracq.algorithms.input_model.oracles`` — :obj:`declare <oracq.algorithms.input_model.oracles.declare>`

:doc:`Quantum data structures (qsample and sample-and-query) <algorithms/input_model/qdata>`
    ``oracq.algorithms.input_model.qdata`` — :obj:`QMatrix <oracq.algorithms.input_model.qdata.QMatrix>`, :obj:`QVector <oracq.algorithms.input_model.qdata.QVector>`

:doc:`Quantum linear systems <algorithms/qlss/qlss>`
    ``oracq.algorithms.qlss.qlss`` — :obj:`BlockSystem <oracq.algorithms.qlss.qlss.BlockSystem>`, :obj:`LinearSystem <oracq.algorithms.qlss.qlss.LinearSystem>`, :obj:`QLSSProtocol <oracq.algorithms.qlss.qlss.QLSSProtocol>`, :obj:`SolveResult <oracq.algorithms.qlss.qlss.SolveResult>`, :obj:`SparseSystem <oracq.algorithms.qlss.qlss.SparseSystem>`, :obj:`SpectralPromise <oracq.algorithms.qlss.qlss.SpectralPromise>`, :obj:`QLSSSolver <oracq.algorithms.qlss.qlss.QLSSSolver>`

:doc:`KP quantum recommendation system <algorithms/qml/recommendation>`
    ``oracq.algorithms.qml.recommendation`` — :obj:`KPRecommendationConfig <oracq.algorithms.qml.recommendation.KPRecommendationConfig>`, :obj:`RecommendationResult <oracq.algorithms.qml.recommendation.RecommendationResult>`, :obj:`kp_recommendation <oracq.algorithms.qml.recommendation.kp_recommendation>`, :obj:`sigma_from_phase <oracq.algorithms.qml.recommendation.sigma_from_phase>`

:doc:`QODE assembly interface <algorithms/qode/ode>`
    ``oracq.algorithms.qode.ode`` — :obj:`QODESolver <oracq.algorithms.qode.ode.QODESolver>`, :obj:`QODEProblem <oracq.algorithms.qode.ode.QODEProblem>`, :obj:`QODEProtocol <oracq.algorithms.qode.ode.QODEProtocol>`

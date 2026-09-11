# 验证覆盖矩阵

本文档记录"每个算法的证据在哪里"，与 `validation-plan.md` 的"如何证明正确"和 `algorithm-coverage.md` 的"实现了什么"行级对应（见 `validation-plan.md` §6 的分工说明）。`validation.md` / `validation.json` 继续承担版本验收快照的职责；本表描述的是当前开发版本的实时覆盖状态，新增算法或新见证必须同步本表（维护规则见末节）。逐算法的接口与验证详情见各算法页（`../manual/algorithms/index.md`）。

类别与判定准则沿用 `validation-plan.md` §2（C1 精确离散、C2 近似连续、C3 概率分布、C4 启发式/优化、C5 数据访问、C6 组合骨架）；阶段（V1–V4）见同文档 §5。三层证据指结构合法、小规模数值见证、可绑定性（计划总原则 §1）。

## 见证原语（`tests/core/witness.py`）

`tests/core/witness.py` 提供四个跨算法共享的断言原语，所有原语都以 `unittest.TestCase` 实例为首参数，失败抛 `AssertionError`。文件名不带 `test_` 前缀，不会被 unittest discover 收集；正确性自证测试位于 `tests/core/test_witness.py`（`UnitaryWitnessTests`、`UncomputationWitnessTests`、`BindInvariantWitnessTests`、`BlockEqualsWitnessTests` 四个测试类，共 11 个用例）。

| 原语 | 签名 | 用途 | 约定 |
|---|---|---|---|
| `assert_unitary` | `(case, program, *, places=9, samples=None)` | W†W = I 的抽样见证：各基态列归一（Σ\|a\|²=1）且两两正交 | `samples` 显式给出基态下标列表；缺省时小寄存器空间（≤16 个基态）取全部，大空间用 `random.Random(0)` 固定种子抽 16 个。失败信息含违规列与 \|a\|² 或内积 |
| `assert_uncomputation` | `(case, program, *, initial=None, work_registers=None)` | 复净见证：模拟器在 LocalExit 处强制复净，原语把 `ValidationError` 转为 `AssertionError`；`work_registers` 列出的根寄存器在所有非零幅度基态中取值为 0 | 返回 `simulate` 的状态便于调用方继续断言。判定零用阈值 `_ZERO_AMPLITUDE = 1e-9`（模拟器丢 \|a\| < 1e-15 的幅度，复净检查需更宽松） |
| `assert_bind_invariant` | `(case, abstract_program, bindings, *, tolerance=0.0)` | 同一抽象槽位的各候选实现须给出一致的可观察分布 | `abstract_program` 须恰有一个未绑定槽位（用 `unresolved` 解析）；`bindings` 字典 `{标签: Operation}`。`tolerance=0` 时精确对拍；QRAM 等有量化误差的绑定传误差界（如 `0.02`），与 `integration.py` 三层一致性的 1e-12 精度分开 |
| `assert_block_equals` | `(case, be, matrix, *, places=9)` | BE 的 (0,0) 块逐列与稠密矩阵对拍（乘回 `alpha`） | 仅依赖 `block_column` helper：先 simulate 初始基态再把 `(row, 0)` 分支的幅度乘 `be.alpha` 复原第 column 列 |

辅助函数 `block_column(be, column)` 单独导出以便调用方在自定义断言中复用（如 `test_lowrank.py` 的 `test_matches_pauli_encoding_block` 自行构造另一侧 BE 后逐行对拍）。

## 覆盖矩阵

### infrastructure / 数据访问层（C5）

| 模块 | 算法 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `oracles.py` | [XorDatabase](../manual/algorithms/xor-database.md)、[StatePreparation](../manual/algorithms/state-preparation.md)、[SparseAccess](../manual/algorithms/sparse-access.md) | C5 | `tests/core/test_contracts.py:OracleContractTests` | `tests/core/test_contracts.py:OracleContractTests.test_dirty_work_is_rejected_by_algorithm` | 同类 `test_capability_restrictions_survive_annotation_and_wrapping` | 三绑定一致性参数化（V4 铺开） | V4 |
| `data_loading.py` | [Select-Swap QROM](../manual/algorithms/select-swap.md)、[QROM Lookup](../manual/algorithms/qrom-lookup.md) | C5 | `tests/core/test_data_loading.py:DataLoadingTests.test_generated_operations_carry_cost_attributes` | 同类 `test_select_swap_matches_gate_database_per_address` / `test_cost_model_matches_formulas_and_tradeoff_curve` | `gate/qram` 两绑定由 `test_result_invariant_across_partitions_and_qrom_lookup_baseline` 间接覆盖；三绑定一致性参数化待 V4 | 三绑定参数化（V4） | V4 |
| `mathfunc/` | [经典函数→可逆线路](../manual/algorithms/mathfunc.md) | C5 | `tests/core/test_mathfunc.py:MathFunctionTests.test_helper_module_reuse_and_roundtrip` | `test_all_elementary_families_construct` / `test_roe_is_a_compiled_pure_function` | — | 定点量化误差界的显式见证未补 | V2 |
| `prepare_select.py` | [PREPARE/SELECT](../manual/algorithms/prepare-select.md)、[alias](../manual/algorithms/alias-preparation.md) | C5/C1 | `tests/core/test_prepare_select.py:PrepareSelectTests.test_gate_and_qram_prepare_bindings_agree` / `test_pauli_hamiltonian_input_and_single_term_shortcut` | 同类 `test_block_encoding_recovers_hamiltonian_over_alpha` / `test_alias_block_encoding_approaches_hamiltonian` / `test_qram_prepare_block_encoding_approaches_hamiltonian` | `test_gate_and_qram_prepare_bindings_agree` 已覆盖 gate/qram 两绑定；`abstract_prepare → gate_prepare` 由 `test_abstract_prepare_binds_inside_block_encoding` 覆盖 | alias 三绑定一致性（gate/qram 当前缺独立场景） | V4 |

### 精确离散算法（C1）

| 模块 | 算法 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `oracle_algorithms.py` | [DJ](../manual/algorithms/deutsch-jozsa.md)/[BV](../manual/algorithms/bernstein-vazirani.md)/[Simon](../manual/algorithms/simon.md) | C1 | `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests.test_legacy_imports_reference_canonical_objects` | 同类 `test_bernstein_vazirani_recovers_secret_with_affine_bias` / `test_simon_sampling_and_rank_aware_postprocess` | — | — | V1 |
| `fourier.py` | [QFT](../manual/algorithms/qft.md)、[QFT 加法](../manual/algorithms/fourier-addition.md) | C1 | `AlgorithmExpansionTests` | `test_qft_matches_positive_fourier_matrix` / `test_fourier_add_all_basis_inputs` | — | — | V1 |
| `arithmetic.py` | [定点加减乘](../manual/algorithms/fixed-point-arithmetic.md)、[布尔网络](../manual/algorithms/boolean-networks.md) | C1 | `tests/core/test_stage2.py:Stage2StructureTests` | `test_arithmetic_construction_and_basis` / `test_gate_network_is_executable_small_example` | — | — | V1 |
| `number_theory.py` | [模乘](../manual/algorithms/modular-multiplication.md)、[order finding](../manual/algorithms/order-finding.md) | C1 | `AlgorithmExpansionTests` | `test_modular_multiplication_total_permutation_and_order`（覆盖 factors_from_phase 解码路径） | — | 已注明不可外推为完整 Shor（保留原有口径） | — |
| `error_correction.py` | [重复码](../manual/algorithms/repetition-codes.md) | C1 | `AlgorithmExpansionTests` | `test_repetition_codes_preserve_arbitrary_logical_amplitudes` | — | — | V1 |
| `walks.py` | [coined cycle walk](../manual/algorithms/coined-cycle-walk.md) | C1 | `AlgorithmExpansionTests` | `test_ansatz_zero_angles_and_walk_one_step` | — | — | V1 |
| `graph_walks.py` | [邻接 oracle](../manual/algorithms/adjacency-oracle.md)、[Szegedy](../manual/algorithms/szegedy-walk.md)、[MNRS](../manual/algorithms/mnrs-search.md) | C1 | `tests/core/test_graph_walks.py:AdjacencyOracleTests` / `SzegedyWalkTests` / `QuantumWalkSearchTests` / `HittingTimeTests` | `test_gate_adjacency_implements_neighbor_table` / `test_stationary_state_is_fixed_point` / `test_search_amplifies_marked_vertex_on_hypercube` / `test_cycle_hitting_times_match_closed_form` | `test_abstract_handles_bind_to_gate_implementations`（邻接 oracle gate/qram 双绑定） | Johnson 图后需新见证（element distinctness） | V3 |

### 近似连续算法（C2）

| 模块 | 算法 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `transforms.py` | [qubitization walk](../manual/algorithms/qubitization-walk.md)、[QSVT 序列](../manual/algorithms/qsvt-sequence.md)、[OAA](../manual/algorithms/oblivious-amplification.md) | C2 | `tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` / `test_trotter_only_input_does_not_need_block_encoding` | `tests/core/test_qsvt.py:PhaseSynthesisTests.test_convention_matches_qsvt_sequence`（逐点一致 1e-10） | — | — | V1 |
| `qsvt.py` | [相位合成](../manual/algorithms/qsp-phase-synthesis.md)、[求逆](../manual/algorithms/qsvt-matrix-inversion.md)、[过滤](../manual/algorithms/eigenstate-filtering.md)、[HS 模拟](../manual/algorithms/qsvt-hamiltonian-simulation.md)、[定点搜索](../manual/algorithms/fixed-point-search.md) | C2 | `PhaseSynthesisTests.test_input_validation` / `TransformWitnessTests.test_transform_input_validation` | `test_roundtrip_chebyshev` / `test_roundtrip_with_explicit_imaginary_part` / `test_negated_phases_conjugate_polynomial` / `test_matrix_inversion_block` / `test_eigenstate_filter_isolates_eigenvalue` / `test_hamiltonian_simulation_block` / `test_fixed_point_search_phase_properties` / `test_fixed_point_search_circuit_amplifies` | — | 收敛性扫描缺失；病态负例（度数 ≳ 16 或精度过低）的参数化测试已通过 input_validation 局部覆盖，单独的"病态输入必须抛错"批量参数化待补 | V2 |
| `hamiltonian.py` | [Trotter](../manual/algorithms/trotter.md)、[Taylor](../manual/algorithms/taylor-block-encoding.md)、[hamsim](../manual/algorithms/hamiltonian-simulation.md) | C2 | `AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` / `test_trotter_validates_each_term` | 由 `test_trotter_protocol_keeps_phase_and_repeat` 间接覆盖阶段顺序 | — | Trotter 阶数误差率扫描 | V2 |
| `sparse.py` / `block_encoding.py` | [稀疏 BE](../manual/algorithms/sparse-block-encoding.md)、[BE 代数](../manual/algorithms/block-encoding-algebra.md) | C2 | `tests/core/test_language.py:BlockEncodingTests` | 同类 `test_unequal_normalization_and_signed_sum` / `test_complex_coefficients` / `test_product_order_and_independent_signal_spaces` | `test_alpha_survives_ir_serialization` | — | V1 |
| `qlss.py` | [CKS](../manual/algorithms/cks.md)、[Costa](../manual/algorithms/costa-walk.md) | C2 | `tests/core/test_qlss_input_models.py:QLSSInputTests` | `test_signed_sparse_encoding_and_chebyshev` / `test_protocols_consume_different_models_and_scale_kappa` / `test_costa_rhs_reflection_is_independent_of_unitary_extension` / `test_matrix_probe_recovers_scalar_system_norm` | `test_costa_rhs_reflection_is_independent_of_unitary_extension`（U 扩展独立性） | 与 HHL 论文参考值的端到端对拍（可入 catalog 目录） | V3 |
| `ode.py` / `lchs.py` / `cbmd.py` / `schrodingerization.py` / `carleman.py` | [QODE 问题](../manual/algorithms/qode-problem.md)与求解器 [LCHS](../manual/algorithms/lchs.md)、[CBMD](../manual/algorithms/cbmd.md)、[Schrödingerization](../manual/algorithms/schrodingerization.md)、[Carleman](../manual/algorithms/carleman.md) | C2 | `tests/core/test_differential.py:DifferentialStructureTests` | `test_four_methods_keep_input_oracles` / `test_cbmd_plan_records_omitted_terms` / `test_qfvm_uses_raw_flow_and_modular_arithmetic` / `test_local_riemann_updates_only_neighbor_faces_and_tree` | 同类注册与契约检查 | 统一的"解析可解 ODE 族"收敛基准 | V2 |
| `sde.py` | [Fokker–Planck](../manual/algorithms/fokker-planck.md) | C2 | `tests/core/test_sde.py:GeneratorDiscretizationTests` / `StatePreparationTests` / `SolverContractTests` | `test_ou_moments_match_closed_form` / `test_stationary_approaches_boltzmann` / `test_column_sums_vanish` | `SolverContractTests.test_qode_problem_accepted_by_lchs` | SDE 与 LCHS 端到端 catalog 案例 | V3 |
| `density.py` | [纯化](../manual/algorithms/purification.md)、[Gibbs](../manual/algorithms/gibbs-state.md) | C2 | `tests/core/test_density.py:PurificationTests` | 同类 `test_gate_purification_recovers_rho` / `test_maximally_mixed_purification` / `test_classical_tools`；`GibbsTests.test_diagonal_hamiltonian_gibbs` / `test_non_diagonal_hamiltonian_gibbs` | `PurificationTests.test_abstract_purification_binds` | —（V1 新增 `test_error_convergence_decreases` 单调不增扫描 + `test_error_bound_uniform_in_beta` β 网格，每档实测距离 5.6e-4 / 8.7e-5 / 8.7e-5 与 2.7e-6 / 1.5e-5 / 1.9e-4） | V1 |
| `lowrank.py` | [DF](../manual/algorithms/double-factorization.md)/[THC](../manual/algorithms/thc.md) → BE | C2 | `tests/core/test_lowrank.py:DiagonalizeSymmetricTests` / `DoubleFactorizationTests` / `ThcTests` | `test_single_rank_block_equals_hamiltonian` / `test_thc_block_equals_hamiltonian` / `test_feeds_qubitization_walk` | 由 `block_column` helper 在 `test_matches_pauli_encoding_block` 内跨 BE 绑定 | α 紧性见证：V1 新增 `test_alpha_matches_closed_form_eigenvalues`（2×2 闭式独立对拍 places=10）+ `test_df_alpha_tighter_than_pauli`（DF α=1.4 ≤ Pauli α=1.5）+ `test_thc_alpha_matches_hand_computed_bound`（THC α=1.996 独立手算） | V1 |

### 概率分布算法（C3）

| 模块 | 算法 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `estimation.py` | [QPE](../manual/algorithms/qpe.md)、[QAE](../manual/algorithms/qae.md)、[量子计数](../manual/algorithms/quantum-counting.md)、[Hadamard test](../manual/algorithms/hadamard-test.md)/[SWAP test](../manual/algorithms/swap-test.md) | C3 | `tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` | `test_amplitude_estimation_half_probability` / `test_hadamard_test_real_and_imaginary` / `test_swap_test_overlap` | — | QAE 置信区间声明未见证（只验证了点估计） | V2 |
| `search.py` | [Grover](../manual/algorithms/grover.md)、[振幅放大](../manual/algorithms/amplitude-amplification.md) | C3 | `AlgorithmExpansionTests` | `test_amplification_of_known_quarter_probability` | — | — | V1 |
| `integration.py` | [Heinrich 求和](../manual/algorithms/heinrich-summation.md)/[积分](../manual/algorithms/heinrich-integration.md) | C3 | `tests/core/test_integration.py:SumPreparationTests` / `QuantumSumTests` / `RateTests` | `test_mean_on_qae_grid_is_exact` / `test_ramp_mean_within_qae_resolution` / `test_quantum_integral_trapezoid_scale` / `test_heinrich_rate_known_values` | `SumPreparationTests.test_qram_binding_matches_gate` / `test_abstract_loader_binds`（三层一致性） | — | V1 |
| `gradient.py` | [Jordan 梯度](../manual/algorithms/jordan-gradient.md) | C3 | `tests/core/test_gradient.py:GradientTests` | `test_linear_function_exact_two_dimensions` / `test_function_phase_oracle_from_mathfunc` | `test_abstract_oracle_binds_to_gate_implementation` | —（V1 新增 `test_perturbed_linear_concentrates_with_grid_bits` 失败概率衰减率断言 q_{m+1} ≤ 0.34·q_m，实测比率 0.292 / 0.268） | V1 |
| `qpca.py` | [QPCA](../manual/algorithms/qpca.md)、[密度矩阵指数化](../manual/algorithms/density-matrix-exponentiation.md) | C3 | `tests/core/test_qpca.py:DensityMatrixExponentiationTests` / `QpcaTests` | `test_pure_state_eigenvalues` / `test_mixed_state_eigenvalue_via_purification` / `test_small_step_first_order_accurate` / `test_error_halves_with_copies` | 由 `gate_purification.as_state_preparation()` 适配路径覆盖 | — | V1 |

### 启发式/优化算法（C4）

| 模块 | 算法 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `variational.py` | [QAOA-MaxCut](../manual/algorithms/qaoa-maxcut.md)、[VQE](../manual/algorithms/vqe.md)、[ansatz](../manual/algorithms/variational-ansatz.md) | C4 | `AlgorithmExpansionTests` | `test_qaoa_single_edge_optimal_layer` / `test_vqe_measurement_of_y_eigenstate` / `test_ansatz_zero_angles_and_walk_one_step` | — | "达到小实例最优割"的强见证 | V3 |
| `dqi.py` | [DQI](../manual/algorithms/dqi.md) | C4 | `tests/core/test_dqi.py:DqiTests` | `test_planted_instance_beats_random_guessing` / `test_identity_instance_with_weight_two` / `test_dicke_state_weight_and_uniformity` | `test_abstract_decoder_binds_to_witness`（分布逐点对拍 Krawtchouk） | identity 边界已按 ≥/≤ 容差修复 | V1 |

### 应用层

| 模块 | 内容 | 类别 | 结构证据 | 数值证据 | 绑定证据 | 缺口 | 阶段 |
|---|---|---|---|---|---|---|---|
| `qfvm.py` / `qham/` | [QFVM](../manual/algorithms/qfvm.md)、[QHAM](../manual/algorithms/qham.md) | C6 | `tests/core/test_qham_general.py:QhamGeneralTests` | `test_rules_roundtrip_and_unsupported_nonlinearity` / `test_lazy_closure_and_rank_unrank` / `test_quadratic_dimension_and_nontrivial_eta_weight` / `test_multidimensional_coupled_identity` / `test_m1_reduces_to_previous_special_case` | `test_open_qode_inputs_keep_operator_and_initial_oracles` | — | V1 |
| `catalog.py` | 33 案例对拍目录 | C6 | `tests/core/test_workloads.py:WorkloadTests` | `test_every_catalog_case_has_a_closed_description` | `test_alpha_change_requires_regeneration` | 新算法（qsvt 求逆、MNRS 搜索、integration、qpca）未登记 | V3 |

### 横切机制

| 机制 | 实现位置 | 自证测试 |
|---|---|---|
| 见证原语 | `tests/core/witness.py` | `tests/core/test_witness.py:UnitaryWitnessTests` / `UncomputationWitnessTests` / `BindInvariantWitnessTests` / `BlockEqualsWitnessTests`（共 11 个用例，含故意错误的反例） |
| 核心静态与单元测试 | `tests/core/test_*.py` | 见上各模块结构/数值/绑定证据列 |

## 维护规则

- 新增算法或新见证时，必须同步本表的模块/算法/缺口/阶段列以及见证原语一节（若引入新原语）。
- 新增算法必须同步新增 `docs/manual/algorithms/` 页面（模板见 `docs/development/writing-docs.md` 的“算法页面”小节）。
- 表格写法与 `algorithm-coverage.md` 一致：模块名按 `src/pyqecclang/` 下的实际文件；算法名沿用对应模块导出的公开符号。
- 类别（C1–C6）按 `validation-plan.md` §2 的判定准则填写；同一算法在不同见证下属于不同类别时，取主导类别并在缺口列说明。
- 三层证据的位置统一写成 `tests/core/<file>.py:<TestClass>` 形式（必要时附 `.<test_method>` 子串），便于读者直接跳转到测试类。
- 阶段（V1–V4）按 `validation-plan.md` §5 填写；不属于任何阶段的注明"—"。

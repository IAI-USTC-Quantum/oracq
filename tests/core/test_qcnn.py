"""QCNN（arXiv:1911.01117）构件的语义与端到端验证。

验证分四层：经典基底（im2col/前向/池化/QRAM 树）逐元素精确；QRAM 行
制备与角度量化镜像一致；Hadamard 内积电路的测量概率满足 Eq. (20)；
采样驱动的输出满足 Eq. (34)-(39) 的语义（max 池化恢复区域最大值、
eta 阈值置零、average 池化收敛到 f^2 加权期望）。
"""

import math
import random
import unittest

from pyqecclang import Builder, simulate
from pyqecclang.algorithms.qcnn import (
    QCNNQRAM,
    ConvSpec,
    cap_relu,
    convolution_forward,
    im2col,
    kernel_columns,
)
from pyqecclang.algorithms.qcnn_layer import (
    qcnn_inner_product,
    qcnn_sampled_layer,
    qcnn_vector_prep,
    quantized_prepared_state,
    vector_angle_tables,
)


def _probe(operation, memory, *, initial=None):
    b = Builder("probe", {r.name: r.type for r in operation.module.registers},
                {r.name: r.type for r in operation.module.resources})
    for name, value in (initial or {}).items():
        for bit in range(b[name].width):
            if (value >> bit) & 1:
                b.x(b[name][bit])
    b.call(operation, **{r.name: b[r.name] for r in operation.module.registers},
           resources={r.name: r.name for r in operation.module.resources})
    return simulate(b.finish().program(), memory=memory, max_steps=2_000_000).amplitudes


def _make_case(seed=7, h=5, w=5):
    spec = ConvSpec(input_shape=(h, w, 1), kernel_shape=(2, 2, 1, 2), cap=3.0, pool=2)
    rng = random.Random(seed)
    x = [rng.uniform(-1, 1) for _ in range(h * w)]
    kernel = [rng.uniform(-1, 1) for _ in range(2 * 2 * 1 * 2)]
    return spec, x, kernel


class ClassicalSubstrateTests(unittest.TestCase):
    def test_im2col_matches_naive_convolution(self):
        spec, x, kernel = _make_case()
        a = im2col(x, spec)
        f = kernel_columns(kernel, spec)
        oh, ow, dp = spec.output_shape
        for i in range(oh):
            for j in range(ow):
                for q in range(dp):
                    naive = sum(
                        x[(i + ki) * 5 + (j + kj)] * kernel[((ki * 2 + kj) * 1) * dp + q]
                        for ki in range(2)
                        for kj in range(2)
                    )
                    got = sum(ar * fc for ar, fc in zip(a[i * ow + j], f[q], strict=True))
                    self.assertAlmostEqual(naive, got, places=12)

    def test_forward_with_pooling_matches_naive(self):
        spec, x, kernel = _make_case()
        out = convolution_forward(x, kernel, spec)
        ph, pw, dp = spec.pooled_shape
        for pi in range(ph):
            for pj in range(pw):
                for q in range(dp):
                    region = [
                        cap_relu(
                            sum(
                                x[(pi * 2 + di + ki) * 5 + (pj * 2 + dj + kj)]
                                * kernel[((ki * 2 + kj) * 1) * dp + q]
                                for ki in range(2)
                                for kj in range(2)
                            ),
                            spec.cap,
                        )
                        for di in range(2)
                        for dj in range(2)
                    ]
                    self.assertAlmostEqual(max(region), out[pi * pw + pj][q], places=12)

    def test_qram_trees_and_pooling_overwrite(self):
        spec, x, kernel = _make_case()
        qram = QCNNQRAM(im2col(x, spec))
        for p, row in enumerate(im2col(x, spec)):
            self.assertEqual(qram.row(p), row)
            self.assertAlmostEqual(
                qram.norm(p), math.sqrt(sum(v * v for v in row)), places=12
            )
        small = QCNNQRAM([[0.0, 0.0]])
        small.update_with_pooling(0, 0, 2.0, "max")
        small.update_with_pooling(0, 0, 1.0, "max")
        self.assertEqual(small.row(0)[0], 2.0)
        avg = QCNNQRAM([[0.0]])
        avg.update_with_pooling(0, 0, 2.0, "average", state=0)
        avg.update_with_pooling(0, 0, 4.0, "average", state=1)
        self.assertAlmostEqual(avg.row(0)[0], 3.0, places=12)


class VectorPrepTests(unittest.TestCase):
    ROWS = [[1.0, 2.0, -1.0, 0.5], [0.5, 0.5, 0.5, 0.5]]

    def test_matches_normalized_rows(self):
        angles, signs = vector_angle_tables(self.ROWS, 10)
        prep = qcnn_vector_prep(2, 2, 10)
        for index, row in enumerate(self.ROWS):
            state = _probe(prep, {"angles__table": angles, "signs__table": signs},
                           initial={"index": index})
            norm = math.sqrt(sum(v * v for v in row))
            for key, amp in state.items():
                self.assertAlmostEqual(amp, row[key[1]] / norm, places=3)

    def test_quantized_mirror_agrees_with_circuit(self):
        angles, signs = vector_angle_tables(self.ROWS, 10)
        prep = qcnn_vector_prep(2, 2, 10)
        for index, row in enumerate(self.ROWS):
            state = _probe(prep, {"angles__table": angles, "signs__table": signs},
                           initial={"index": index})
            mirror = quantized_prepared_state(row, 10)
            for key, amp in state.items():
                self.assertAlmostEqual(amp, mirror[key[1]], places=9)


class InnerProductTests(unittest.TestCase):
    def test_measurement_probabilities_satisfy_equation_20(self):
        rows = [[1.0, 2.0, -1.0, 0.5], [0.5, 0.5, 0.5, 0.5]]
        cols = [[0.3, -0.7, 0.2, 0.4], [1.0, 0.0, 0.0, 0.0]]
        operation = qcnn_inner_product(2, 2, 2, 10)
        memory = {}
        for prefix, vectors in (("row", rows), ("col", cols)):
            angles, signs = vector_angle_tables(vectors, 10)
            memory[prefix + "__angles__table"] = angles
            memory[prefix + "__signs__table"] = signs
        state = _probe(operation, memory)
        prob = {}
        for key, amp in state.items():
            triplet = (key[0], key[1], key[2])
            prob[triplet] = prob.get(triplet, 0.0) + abs(amp) ** 2

        def normalized(a):
            norm = math.sqrt(sum(v * v for v in a))
            return [v / norm for v in a]

        for p in range(2):
            for q in range(2):
                overlap = sum(
                    u * v for u, v in zip(normalized(rows[p]), normalized(cols[q]), strict=True)
                )
                expected = (1 + overlap) / (2 * 2 * 2)
                self.assertAlmostEqual(
                    prob.get((p, q, 0), 0.0), expected, places=4, msg=f"p={p} q={q}"
                )


class SampledLayerTests(unittest.TestCase):
    def test_max_pooling_recovers_region_maxima(self):
        spec, x, kernel = _make_case()
        classical = convolution_forward(x, kernel, spec)
        pooled, _ = qcnn_sampled_layer(x, kernel, spec, samples=40000, eta=0.0, seed=1)
        for row_c, row_q in zip(classical, pooled, strict=True):
            for expected, got in zip(row_c, row_q, strict=True):
                self.assertLessEqual(abs(expected - got), 5e-3)

    def test_eta_threshold_zeroes_everything(self):
        spec, x, kernel = _make_case()
        pooled, stats = qcnn_sampled_layer(x, kernel, spec, samples=1000, eta=1e9, seed=1)
        self.assertTrue(all(v == 0 for row in pooled for v in row))
        self.assertEqual(stats["kept"], 0)

    def test_average_pooling_converges_to_quadratic_mean(self):
        spec = ConvSpec(
            input_shape=(5, 5, 1), kernel_shape=(2, 2, 1, 1), cap=3.0, pool=2,
            pool_kind="average",
        )
        rng = random.Random(7)
        x = [rng.uniform(-1, 1) for _ in range(25)]
        kernel = [rng.uniform(-1, 1) for _ in range(4)]
        pooled, _ = qcnn_sampled_layer(x, kernel, spec, samples=60000, eta=0.0, seed=2)
        rows, cols = im2col(x, spec), kernel_columns(kernel, spec)
        oh, ow = spec.output_shape[:2]
        ph, pw, d = spec.pooled_shape
        for pi in range(ph):
            for pj in range(pw):
                for q in range(d):
                    num = den = 0.0
                    for di in range(2):
                        for dj in range(2):
                            p = (pi * 2 + di) * ow + (pj * 2 + dj)
                            y = cap_relu(
                                sum(ar * fc for ar, fc in zip(rows[p], cols[q], strict=True)),
                                spec.cap,
                            )
                            num += y**3
                            den += y * y
                    expectation = num / den if den else 0.0
                    self.assertLessEqual(
                        abs(pooled[pi * pw + pj][q] - expectation), 2e-2,
                        msg=f"({pi},{pj},{q})",
                    )


if __name__ == "__main__":
    unittest.main()

"算法库共享的 Python 结构协议；允许应用在自己的模块中继续定义协议。"

from __future__ import annotations

from inspect import signature
from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

from oracq.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractError,
    InputRequirement,
    fail,
    require_instance,
)

if TYPE_CHECKING:
    from oracq.algorithms.input_model.operators import BlockEncoding
    from oracq.algorithms.input_model.oracles import (
        SparseAccess,
        StateOracle,
        StatePreparation,
        XorDatabase,
    )
    from oracq.infrastructure.builder import Operation


@runtime_checkable
class StatePreparationProtocol(Protocol):
    """要求 ``state_preparation()`` 方法；满足即可作为算法的初态输入。"""

    def state_preparation(self) -> StatePreparation:
        """返回该输入对应的 ``StatePreparation`` 访问视图。

        Returns:
            StatePreparation: 该输入在零输入寄存器上制备目标态的视图。
        """
        ...


@runtime_checkable
class UnitaryProtocol(StatePreparationProtocol, Protocol):
    """在态制备协议之上增加 ``unitary()`` 的完整酉接口。"""

    def unitary(self) -> Operation:
        """返回完整寄存器空间上的酉 ``Operation``。

        Returns:
            Operation: 作用在该输入全部寄存器上的酉操作。
        """
        ...


@runtime_checkable
class BlockEncodingProtocol(Protocol):
    """要求 ``block_encoding()`` 方法；满足即可作为块编码输入。"""

    def block_encoding(self) -> BlockEncoding:
        """返回该输入对应的 ``BlockEncoding`` 访问视图。

        Returns:
            BlockEncoding: 把该输入编码为酉的正规化子块的访问视图。
        """
        ...


@runtime_checkable
class CKSSparseProtocol(Protocol):
    """要求 ``sparse_access()`` 方法；满足即可作为 CKS 稀疏矩阵输入。"""

    def sparse_access(self) -> SparseAccess:
        """返回 CKS 位置与元素操作组成的 ``SparseAccess``。

        Returns:
            SparseAccess: 由位置 oracle 与元素 oracle 组成的稀疏访问视图。
        """
        ...


@runtime_checkable
class XorDatabaseProtocol(Protocol):
    """要求 ``xor_database()`` 方法；满足即可作为 XOR 数据库输入。"""

    def xor_database(self) -> XorDatabase:
        """返回该输入对应的 ``XorDatabase`` 访问视图。

        Returns:
            XorDatabase: 以 XOR 方式应答查询的数据库访问视图。
        """
        ...


@runtime_checkable
class StateOracleProtocol(Protocol):
    """要求 ``state_oracle()`` 方法；满足即可作为输出态 oracle。"""

    def state_oracle(self) -> StateOracle:
        """返回含成功信号的 ``StateOracle``。

        Returns:
            StateOracle: 在输出态中标注成功子空间的 oracle。
        """
        ...


_View = TypeVar("_View")


def _view(value: object, method_name: str, cls: type[_View]) -> _View:
    """只取得一次访问视图；调用边界错误与提供方内部异常分别处理。"""
    if isinstance(value, cls):
        return value
    method = getattr(value, method_name, None)
    if method is None:
        fail("INPUT_PROTOCOL", method_name, "provider method", type(value).__name__, "缺少访问方法")
    if not callable(method):
        fail("INPUT_CALLABLE", method_name, "callable", type(method).__name__, "访问属性不可调用")
    try:
        parameters = signature(method)
    except (TypeError, ValueError):
        # 原生可调用对象可能不暴露签名；实际调用的异常保留原始 traceback。
        parameters = None
    if parameters is not None:
        try:
            parameters.bind()
        except TypeError as exc:
            from oracq.algorithms.input_model.contracts import ContractIssue

            raise ContractError((ContractIssue(
                "INPUT_SIGNATURE", method_name, "zero-argument call", str(parameters),
                "访问方法必须支持无参数调用",
            ),)) from exc
    return require_instance(method(), cls, method_name + ".output")


def as_state_preparation(value: StatePreparationProtocol) -> StatePreparation:
    """把满足协议的输入转换成 ``StatePreparation``。

    Args:
        value: 满足 ``StatePreparationProtocol`` 的对象。

    Returns:
        StatePreparation: 调用 ``state_preparation()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``StatePreparation`` 时以 ``INPUT_TYPE`` 上报。
    """
    from oracq.algorithms.input_model.oracles import StatePreparation

    return _view(value, "state_preparation", StatePreparation)


def checked_state_preparation(
    value: StatePreparationProtocol,
    *,
    adjoint: bool = False,
    controlled: bool = False,
    path: str = "preparation",
) -> StatePreparation:
    """取得零输入、干净工作区的制备；由调用算法选择需要的变换能力。

    Args:
        value: 满足 ``StatePreparationProtocol`` 的输入对象。
        adjoint: 制备是否必须提供厄米共轭版本。
        controlled: 制备是否必须支持受控版本。
        path: 契约报告中该输入的路径名。

    Returns:
        StatePreparation: 通过零输入与干净工作区检查的制备视图。
    """
    result = as_state_preparation(value)
    requirement = InputRequirement(
        path,
        (StatePreparationProtocol,),
        zero_input=True,
        clean_work=True,
        adjoint=adjoint,
        controlled=controlled,
    )
    _, issues = requirement.inspect(result)
    if issues:
        raise ContractError(issues)
    return result


def as_block_encoding(value: BlockEncodingProtocol) -> BlockEncoding:
    """把满足协议的输入转换成 ``BlockEncoding``。

    Args:
        value: 满足 ``BlockEncodingProtocol`` 的对象。

    Returns:
        BlockEncoding: 调用 ``block_encoding()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``BlockEncoding`` 时以 ``INPUT_TYPE`` 上报。
    """
    from oracq.algorithms.input_model.operators import BlockEncoding

    return _view(value, "block_encoding", BlockEncoding)


def as_sparse_access(value: CKSSparseProtocol) -> SparseAccess:
    """把满足协议的输入转换成 ``SparseAccess``。

    Args:
        value: 满足 ``CKSSparseProtocol`` 的对象。

    Returns:
        SparseAccess: 调用 ``sparse_access()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``SparseAccess`` 时以 ``INPUT_TYPE`` 上报。
    """
    from oracq.algorithms.input_model.oracles import SparseAccess

    return _view(value, "sparse_access", SparseAccess)


def as_qlss_matrix(
    value: BlockEncodingProtocol | CKSSparseProtocol,
) -> BlockEncoding | SparseAccess:
    """按可用接口把矩阵输入转换成 ``BlockEncoding`` 或 ``SparseAccess``。

    满足 ``BlockEncodingProtocol`` 时优先返回块编码视图，否则按
    ``CKSSparseProtocol`` 取稀疏访问视图。

    Args:
        value: 矩阵输入对象。

    Returns:
        满足块编码协议时为 ``BlockEncoding``，否则为 ``SparseAccess``。

    Raises:
        ContractError: 两个协议都不满足，或取得的返回值类型不符。
    """
    return (
        as_block_encoding(value)
        if isinstance(value, BlockEncodingProtocol)
        else as_sparse_access(value)
    )


def operator_state_contract(
    name: str,
    *,
    matrix: str = "generator",
    state: str = "initial",
    composable: bool = True,
) -> AlgorithmContract:
    """当前 BE 型演化算法的需求；其他算法可定义完全不同的契约。

    Args:
        name: 契约名，用于报告与诊断。
        matrix: 块编码矩阵输入在契约中的路径名。
        state: 初态制备输入在契约中的路径名。
        composable: 初态是否要求厄米共轭与受控能力，供组合式算法使用。

    Returns:
        AlgorithmContract: 要求矩阵为可逆可控的块编码、初态为零输入且
        干净工作区的制备，并约束两者寄存器等宽的契约。
    """
    return AlgorithmContract(
        name,
        (
            InputRequirement(
                matrix,
                (BlockEncodingProtocol,),
                adapter=as_block_encoding,
                adjoint=True,
                controlled=True,
            ),
            InputRequirement(
                state,
                (StatePreparationProtocol,),
                adapter=as_state_preparation,
                adjoint=composable,
                controlled=composable,
                zero_input=True,
                clean_work=True,
            ),
        ),
        same_width=((matrix, state),),
    )

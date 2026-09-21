"算法库共享的 Python 结构协议；允许应用在自己的模块中继续定义协议。"

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

from pyqecclang.algorithms.input_model.contracts import (
    ContractError,
    InputRequirement,
    ProtocolContract,
    require_instance,
    requires,
)

if TYPE_CHECKING:
    from pyqecclang.algorithms.input_model.operators import BlockEncoding
    from pyqecclang.algorithms.input_model.oracles import (
        SparseAccess,
        StateOracle,
        StatePreparation,
        XorDatabase,
    )
    from pyqecclang.infrastructure.builder import Operation


@runtime_checkable
class StatePreparationProtocol(Protocol):
    """要求 ``state_preparation()`` 方法；满足即可作为算法的初态输入。"""

    def state_preparation(self) -> StatePreparation:
        """返回该输入对应的 ``StatePreparation`` 访问视图。"""
        ...


@runtime_checkable
class UnitaryProtocol(StatePreparationProtocol, Protocol):
    """在态制备协议之上增加 ``unitary()`` 的完整酉接口。"""

    def unitary(self) -> Operation:
        """返回完整寄存器空间上的酉 ``Operation``。"""
        ...


@runtime_checkable
class BlockEncodingProtocol(Protocol):
    """要求 ``block_encoding()`` 方法；满足即可作为块编码输入。"""

    def block_encoding(self) -> BlockEncoding:
        """返回该输入对应的 ``BlockEncoding`` 访问视图。"""
        ...


@runtime_checkable
class CKSSparseProtocol(Protocol):
    """要求 ``sparse_access()`` 方法；满足即可作为 CKS 稀疏矩阵输入。"""

    def sparse_access(self) -> SparseAccess:
        """返回 CKS 位置与元素操作组成的 ``SparseAccess``。"""
        ...


@runtime_checkable
class XorDatabaseProtocol(Protocol):
    """要求 ``xor_database()`` 方法；满足即可作为 XOR 数据库输入。"""

    def xor_database(self) -> XorDatabase:
        """返回该输入对应的 ``XorDatabase`` 访问视图。"""
        ...


@runtime_checkable
class StateOracleProtocol(Protocol):
    """要求 ``state_oracle()`` 方法；满足即可作为输出态 oracle。"""

    def state_oracle(self) -> StateOracle:
        """返回含成功信号的 ``StateOracle``。"""
        ...


def as_state_preparation(value):
    """把满足协议的输入转换成 ``StatePreparation``。

    Args:
        value: 满足 ``StatePreparationProtocol`` 的对象。

    Returns:
        StatePreparation: 调用 ``state_preparation()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``StatePreparation`` 时以 ``INPUT_TYPE`` 上报。
    """
    from pyqecclang.algorithms.input_model.oracles import StatePreparation

    result = requires(value, StatePreparationProtocol).state_preparation()
    require_instance(result, StatePreparation, "state_preparation.output")
    return result


def checked_state_preparation(value, *, adjoint=False, controlled=False, path="preparation"):
    """取得零输入、干净工作区的制备；由调用算法选择需要的变换能力。"""
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


def as_block_encoding(value):
    """把满足协议的输入转换成 ``BlockEncoding``。

    Args:
        value: 满足 ``BlockEncodingProtocol`` 的对象。

    Returns:
        BlockEncoding: 调用 ``block_encoding()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``BlockEncoding`` 时以 ``INPUT_TYPE`` 上报。
    """
    from pyqecclang.algorithms.input_model.operators import BlockEncoding

    result = requires(value, BlockEncodingProtocol).block_encoding()
    require_instance(result, BlockEncoding, "block_encoding.output")
    return result


def as_sparse_access(value):
    """把满足协议的输入转换成 ``SparseAccess``。

    Args:
        value: 满足 ``CKSSparseProtocol`` 的对象。

    Returns:
        SparseAccess: 调用 ``sparse_access()`` 得到的具体视图。

    Raises:
        ContractError: 缺少接口时以 ``INPUT_PROTOCOL`` 上报；返回值不是
            ``SparseAccess`` 时以 ``INPUT_TYPE`` 上报。
    """
    from pyqecclang.algorithms.input_model.oracles import SparseAccess

    result = requires(value, CKSSparseProtocol).sparse_access()
    require_instance(result, SparseAccess, "sparse_access.output")
    return result


def as_qlss_matrix(value):
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


def operator_state_contract(name, *, matrix="generator", state="initial", composable=True):
    """当前 BE 型演化算法的需求；其他算法可定义完全不同的契约。"""
    return ProtocolContract(
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

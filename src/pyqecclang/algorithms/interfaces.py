"算法库共享的 Python 结构协议；允许应用在自己的模块中继续定义协议。"

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

from pyqecclang.algorithms.contracts import (
    ContractError,
    InputRequirement,
    ProtocolContract,
    require_instance,
    requires,
)

if TYPE_CHECKING:
    from pyqecclang.algorithms.operators import BlockEncoding
    from pyqecclang.algorithms.oracles import (
        SparseAccess,
        StateOracle,
        StatePreparation,
        XorDatabase,
    )
    from pyqecclang.infrastructure.builder import Operation


@runtime_checkable
class StatePreparationProtocol(Protocol):
    def state_preparation(self) -> StatePreparation: ...


@runtime_checkable
class UnitaryProtocol(StatePreparationProtocol, Protocol):
    def unitary(self) -> Operation: ...


@runtime_checkable
class BlockEncodingProtocol(Protocol):
    def block_encoding(self) -> BlockEncoding: ...


@runtime_checkable
class CKSSparseProtocol(Protocol):
    def sparse_access(self) -> SparseAccess: ...


@runtime_checkable
class XorDatabaseProtocol(Protocol):
    def xor_database(self) -> XorDatabase: ...


@runtime_checkable
class StateOracleProtocol(Protocol):
    def state_oracle(self) -> StateOracle: ...


def as_state_preparation(value):
    from pyqecclang.algorithms.oracles import StatePreparation

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
    from pyqecclang.algorithms.operators import BlockEncoding

    result = requires(value, BlockEncodingProtocol).block_encoding()
    require_instance(result, BlockEncoding, "block_encoding.output")
    return result


def as_sparse_access(value):
    from pyqecclang.algorithms.oracles import SparseAccess

    result = requires(value, CKSSparseProtocol).sparse_access()
    require_instance(result, SparseAccess, "sparse_access.output")
    return result


def as_qlss_matrix(value):
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

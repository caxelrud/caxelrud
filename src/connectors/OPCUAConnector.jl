module OPCUAConnector

using PythonCall

import ..Config: OPCUAConfig

export OPCUAHandle, connect_opcua, disconnect_opcua!, read_values, read_value, write_value!

# Julia has no mature native OPC UA client, so this connector bridges to the
# well-established Python `asyncua` library (its synchronous wrapper) via
# PythonCall/CondaPkg. The dependency (and a private Python environment) is
# declared in CondaPkg.toml and provisioned automatically on first use.
const _asyncua_sync = Ref{Py}()

function __init__()
    _asyncua_sync[] = pyimport("asyncua.sync")
end

"""
    OPCUAHandle

Wraps a connected `asyncua.sync.Client` plus the node IDs configured for
polling.
"""
mutable struct OPCUAHandle
    client::Py
    node_ids::Vector{String}
end

"""
    connect_opcua(cfg::OPCUAConfig) -> OPCUAHandle

Connect to the OPC UA server at `cfg.endpoint`, optionally authenticating
with `cfg.username`/`cfg.password`.
"""
function connect_opcua(cfg::OPCUAConfig)
    client = _asyncua_sync[].Client(cfg.endpoint)

    if !isempty(cfg.username)
        client.set_user(cfg.username)
        isempty(cfg.password) || client.set_password(cfg.password)
    end

    client.connect()
    return OPCUAHandle(client, cfg.node_ids)
end

"""
    disconnect_opcua!(handle::OPCUAHandle)
"""
function disconnect_opcua!(handle::OPCUAHandle)
    handle.client.disconnect()
    return nothing
end

"""
    read_value(handle::OPCUAHandle, node_id) -> Any

Read a single node's current value, converted to a native Julia value.
"""
function read_value(handle::OPCUAHandle, node_id::AbstractString)
    node = handle.client.get_node(node_id)
    return pyconvert(Any, node.read_value())
end

"""
    read_values(handle::OPCUAHandle, node_ids = handle.node_ids) -> Dict{String,Any}

Read several nodes at once, keyed by node ID string.
"""
function read_values(handle::OPCUAHandle, node_ids::AbstractVector{<:AbstractString} = handle.node_ids)
    return Dict(nid => read_value(handle, nid) for nid in node_ids)
end

"""
    write_value!(handle::OPCUAHandle, node_id, value)

Write `value` to the given node (e.g. to change a setpoint from the
dashboard).
"""
function write_value!(handle::OPCUAHandle, node_id::AbstractString, value)
    node = handle.client.get_node(node_id)
    node.write_value(value)
    return nothing
end

end # module

module MQTTConnector

using Mosquitto

import ..Config: MQTTConfig

export MQTTHandle, connect_mqtt, disconnect_mqtt!, publish_value, latest_values, latest_value

"""
    MQTTHandle

Wraps a `Mosquitto.Client` together with a thread-safe cache of the latest
payload received on each subscribed topic, kept up to date by a background
task. Read it from the dashboard with [`latest_value`](@ref) /
[`latest_values`](@ref).
"""
mutable struct MQTTHandle
    client::Mosquitto.Client
    cache::Dict{String,String}
    lock::ReentrantLock
    running::Bool
end

"""
    connect_mqtt(cfg::MQTTConfig; qos = 1, poll_interval = 0.1) -> MQTTHandle

Connect to the MQTT broker described by `cfg`, subscribe to `cfg.topics`,
and start a background task that keeps the latest message per topic
available via [`latest_value`](@ref).

The network I/O itself runs on Mosquitto's own C thread (`loop_start`); the
Julia task only drains the incoming-message channel, so it never blocks the
Pluto notebook.
"""
function connect_mqtt(cfg::MQTTConfig; qos::Int = 1, poll_interval::Real = 0.1)
    client = Mosquitto.Client(; id = cfg.client_id)

    if cfg.use_tls && !isempty(cfg.cafile)
        Mosquitto.tls_set(client, cfg.cafile)
    end

    flag = Mosquitto.connect(client, cfg.host, cfg.port; username = cfg.username, password = cfg.password)
    flag == 0 || @warn "MQTT connect to $(cfg.host):$(cfg.port) returned error code $flag"

    Mosquitto.loop_start(client)

    for topic in cfg.topics
        Mosquitto.subscribe(client, topic; qos = qos)
    end

    handle = MQTTHandle(client, Dict{String,String}(), ReentrantLock(), true)
    @async _drain_loop(handle, poll_interval)
    return handle
end

function _drain_loop(handle::MQTTHandle, poll_interval::Real)
    msg_channel = Mosquitto.get_messages_channel(handle.client)
    while handle.running
        while isready(msg_channel)
            msg = take!(msg_channel)
            lock(handle.lock) do
                handle.cache[msg.topic] = String(msg.payload)
            end
        end
        sleep(poll_interval)
    end
end

"""
    disconnect_mqtt!(handle::MQTTHandle)

Stop the background drain task, disconnect from the broker and stop
Mosquitto's network thread.
"""
function disconnect_mqtt!(handle::MQTTHandle)
    handle.running = false
    Mosquitto.disconnect(handle.client)
    Mosquitto.loop_stop(handle.client)
    return nothing
end

"""
    publish_value(handle, topic, payload; retain = false, qos = 1)

Publish `payload` (converted with `string`) to `topic`.
"""
function publish_value(handle::MQTTHandle, topic::AbstractString, payload; retain::Bool = false, qos::Int = 1)
    return Mosquitto.publish(handle.client, topic, string(payload); qos = qos, retain = retain)
end

"""
    latest_values(handle::MQTTHandle) -> Dict{String,String}

A snapshot copy of the latest payload received per topic.
"""
function latest_values(handle::MQTTHandle)
    lock(handle.lock) do
        return copy(handle.cache)
    end
end

"""
    latest_value(handle::MQTTHandle, topic) -> Union{String,Missing}

The latest payload received on `topic`, or `missing` if none has arrived yet.
"""
function latest_value(handle::MQTTHandle, topic::AbstractString)
    lock(handle.lock) do
        return get(handle.cache, topic, missing)
    end
end

end # module

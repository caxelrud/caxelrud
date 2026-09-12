module Config

using TOML

export AppConfig, InfluxConfig, MQTTConfig, OPCUAConfig, load_config

"""
    InfluxConfig

Connection settings for an InfluxDB v2 instance (Flux query API).
"""
Base.@kwdef struct InfluxConfig
    url::String = "http://localhost:8086"
    org::String = ""
    bucket::String = ""
    token::String = ""
end

"""
    MQTTConfig

Connection settings for an MQTT broker, plus the list of topics to
subscribe to on connect.
"""
Base.@kwdef struct MQTTConfig
    host::String = "localhost"
    port::Int = 1883
    client_id::String = "hmi-dashboard"
    username::String = ""
    password::String = ""
    topics::Vector{String} = String[]
    use_tls::Bool = false
    cafile::String = ""
end

"""
    OPCUAConfig

Connection settings for an OPC UA server and the node IDs to poll.
"""
Base.@kwdef struct OPCUAConfig
    endpoint::String = "opc.tcp://localhost:4840/freeopcua/server/"
    username::String = ""
    password::String = ""
    node_ids::Vector{String} = String[]
end

"""
    AppConfig

Top-level configuration bundling all three data source configs.
"""
Base.@kwdef struct AppConfig
    influx::InfluxConfig = InfluxConfig()
    mqtt::MQTTConfig = MQTTConfig()
    opcua::OPCUAConfig = OPCUAConfig()
end

# Environment variables always win over the TOML file, so secrets can be
# injected at runtime (docker/CI) without editing/committing the file.
_env(key::AbstractString, default) = get(ENV, key, default)

_strvec(v) = v isa AbstractVector ? String.(v) : String[]

"""
    load_config(path = "config/config.toml") -> AppConfig

Load configuration from a TOML file, if present, then apply environment
variable overrides (`INFLUX_*`, `MQTT_*`, `OPCUA_*`). Missing values fall
back to sensible local-development defaults.
"""
function load_config(path::AbstractString = joinpath("config", "config.toml"))
    raw = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()

    influx_raw = get(raw, "influxdb", Dict{String,Any}())
    mqtt_raw = get(raw, "mqtt", Dict{String,Any}())
    opcua_raw = get(raw, "opcua", Dict{String,Any}())

    influx = InfluxConfig(
        url = _env("INFLUX_URL", get(influx_raw, "url", "http://localhost:8086")),
        org = _env("INFLUX_ORG", get(influx_raw, "org", "")),
        bucket = _env("INFLUX_BUCKET", get(influx_raw, "bucket", "")),
        token = _env("INFLUX_TOKEN", get(influx_raw, "token", "")),
    )

    mqtt = MQTTConfig(
        host = _env("MQTT_HOST", get(mqtt_raw, "host", "localhost")),
        port = parse(Int, string(_env("MQTT_PORT", get(mqtt_raw, "port", 1883)))),
        client_id = _env("MQTT_CLIENT_ID", get(mqtt_raw, "client_id", "hmi-dashboard")),
        username = _env("MQTT_USERNAME", get(mqtt_raw, "username", "")),
        password = _env("MQTT_PASSWORD", get(mqtt_raw, "password", "")),
        topics = _strvec(get(mqtt_raw, "topics", String[])),
        use_tls = get(mqtt_raw, "use_tls", false),
        cafile = _env("MQTT_CAFILE", get(mqtt_raw, "cafile", "")),
    )

    opcua = OPCUAConfig(
        endpoint = _env("OPCUA_ENDPOINT", get(opcua_raw, "endpoint", "opc.tcp://localhost:4840/freeopcua/server/")),
        username = _env("OPCUA_USERNAME", get(opcua_raw, "username", "")),
        password = _env("OPCUA_PASSWORD", get(opcua_raw, "password", "")),
        node_ids = _strvec(get(opcua_raw, "node_ids", String[])),
    )

    return AppConfig(influx = influx, mqtt = mqtt, opcua = opcua)
end

end # module

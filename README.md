# HMI Dashboard (Julia + Pluto)

A lightweight, reactive HMI/SCADA-style dashboard built with
[Pluto.jl](https://plutojl.org/) notebooks, pulling live and historical data
from three common industrial data sources:

- **InfluxDB** (v2, Flux) — historical trends, KPIs
- **MQTT** — live tag values pushed from PLCs/edge gateways/sensors
- **OPC UA** — direct polling/writing of PLC/controller nodes

The connectors are plain Julia modules (no Pluto dependency), so they can
also be reused from scripts, `Genie`/`HTTP.jl` services, or other notebooks.

## Project layout

```
Project.toml            # package dependencies
CondaPkg.toml            # Python env spec (asyncua, for the OPC UA bridge)
src/
  HMIDashboard.jl         # top-level module, ties the pieces together
  Config.jl               # TOML + environment variable configuration
  connectors/
    InfluxConnector.jl     # Flux queries + line-protocol writes
    MQTTConnector.jl        # Mosquitto.jl wrapper with a live tag cache
    OPCUAConnector.jl        # asyncua (Python) bridge via PythonCall
notebooks/
  dashboard.jl            # the actual Pluto dashboard
config/
  config.example.toml     # copy to config/config.toml and fill in
test/
  runtests.jl             # offline unit tests (no live servers required)
```

## Requirements

- Julia 1.9+
- A running InfluxDB v2, MQTT broker, and/or OPC UA server to connect to
  (all three are optional — the dashboard degrades gracefully if a source
  isn't configured)
- For OPC UA: Julia has no mature native OPC UA client, so this project
  bridges to Python's [`asyncua`](https://github.com/FreeOpcUa/opcua-asyncio)
  library via [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl).
  `CondaPkg.jl` provisions a private Python + `asyncua` environment
  automatically the first time it's needed — no manual Python setup required.

## Setup

```sh
julia --project=. -e 'import Pkg; Pkg.instantiate()'
cp config/config.example.toml config/config.toml
# edit config/config.toml with your InfluxDB/MQTT/OPC UA connection details
```

Secrets (tokens, passwords) can instead be supplied via environment
variables, which always override the TOML file and never need to be
committed:

| Variable          | Overrides             |
|-------------------|------------------------|
| `INFLUX_URL`      | `influxdb.url`         |
| `INFLUX_ORG`      | `influxdb.org`         |
| `INFLUX_BUCKET`   | `influxdb.bucket`      |
| `INFLUX_TOKEN`    | `influxdb.token`       |
| `MQTT_HOST`       | `mqtt.host`            |
| `MQTT_PORT`       | `mqtt.port`            |
| `MQTT_CLIENT_ID`  | `mqtt.client_id`       |
| `MQTT_USERNAME`   | `mqtt.username`        |
| `MQTT_PASSWORD`   | `mqtt.password`        |
| `MQTT_CAFILE`     | `mqtt.cafile`          |
| `OPCUA_ENDPOINT`  | `opcua.endpoint`       |
| `OPCUA_USERNAME`  | `opcua.username`       |
| `OPCUA_PASSWORD`  | `opcua.password`       |

`config/config.toml` is git-ignored so real credentials never end up in
version control — always use `config/config.example.toml` as the template.

## Running the dashboard

```sh
julia --project=. -e 'import Pkg; Pkg.add("Pluto"); using Pluto; Pluto.run(notebook = "notebooks/dashboard.jl")'
```

(Or, if you already have Pluto installed globally, just open
`notebooks/dashboard.jl` from the Pluto start page.) The notebook auto-loads
`config/config.toml`, connects to whichever sources you've configured, and
refreshes on a timer via `PlutoUI.Clock`.

## Using the connectors from your own code

```julia
using HMIDashboard

cfg = Config.load_config("config/config.toml")

# InfluxDB
trend_df = InfluxConnector.trend(cfg.influx, "temperature", "value"; range = "-1h")
InfluxConnector.write_record!(cfg.influx, "temperature", ["room" => "lab1"], ["value" => 21.5])

# MQTT
mqtt = MQTTConnector.connect_mqtt(cfg.mqtt)
MQTTConnector.latest_value(mqtt, "plant/line1/temperature")
MQTTConnector.publish_value(mqtt, "plant/line1/setpoint", 22.0)

# OPC UA
opcua = OPCUAConnector.connect_opcua(cfg.opcua)
OPCUAConnector.read_values(opcua)
OPCUAConnector.write_value!(opcua, "ns=2;i=3", 42.0)
```

## Testing

```sh
julia --project=. -e 'import Pkg; Pkg.test()'
```

Tests are self-contained (no live InfluxDB/MQTT/OPC UA server required) and
cover configuration loading, line-protocol encoding, and annotated-CSV
parsing.

## Notes on the OPC UA bridge

Native Julia OPC UA client packages don't currently exist in the General
registry. Rather than implementing the OPC UA binary protocol from scratch,
`OPCUAConnector.jl` bridges to `asyncua`'s synchronous client via
`PythonCall.jl`/`CondaPkg.jl`. This is a common, supported pattern in the
Julia ecosystem for interoperating with mature Python libraries, and keeps
the rest of the stack (InfluxDB, MQTT, the dashboard itself) pure Julia.

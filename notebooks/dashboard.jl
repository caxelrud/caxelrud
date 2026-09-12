### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 88a1b7c8-ae45-11f1-088b-0fda5f21c95c
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 88a1b8a6-ae45-11f1-2d16-23bb63d6e8db
begin
	using HMIDashboard
	using PlutoUI
	using DataFrames
	using Dates
end

# ╔═╡ 88a0fec8-ae45-11f1-3e9f-c99acd70081c
md"""
# 🏭 HMI Dashboard

Live and historical process data from **InfluxDB**, **MQTT**, and **OPC UA**,
in one reactive Pluto notebook.

Configure your data sources in `config/config.toml` (copy it from
`config/config.example.toml`), then just open this notebook — everything
below connects automatically and keeps itself up to date.
"""

# ╔═╡ 88a1b91c-ae45-11f1-0aea-05e06cd61b51
md"""
## Configuration

Loaded from `config/config.toml`, falling back to defaults / environment
variables (`INFLUX_*`, `MQTT_*`, `OPCUA_*`) when a value is missing — see
the README for the full list.
"""

# ╔═╡ 88a1b980-ae45-11f1-3c16-5da6783710f7
cfg = Config.load_config(joinpath(@__DIR__, "..", "config", "config.toml"))

# ╔═╡ 88a1b9f8-ae45-11f1-1d48-893b12b9890e
md"""
## Live refresh

Everything below re-runs automatically on this clock — adjust the interval
to taste (shorter for a faster-moving process, longer to be gentler on the
data sources).
"""

# ╔═╡ 88a1ba8e-ae45-11f1-114d-9f0ea22be44f
@bind tick PlutoUI.Clock(interval=2.0, fixed=false, start_running=true)

# ╔═╡ 88a1bb06-ae45-11f1-37fe-e9d3eec62567
md"""
## MQTT — live tags

Subscribes to the topics listed in `config/config.toml` (`mqtt.topics`) and
shows the latest payload received for each.
"""

# ╔═╡ 88a1bb6a-ae45-11f1-2423-fd1b0a35d242
mqtt_handle = try
	MQTTConnector.connect_mqtt(cfg.mqtt)
catch e
	@warn "Could not connect to the MQTT broker" exception = (e, catch_backtrace())
	nothing
end

# ╔═╡ 88a1bbc4-ae45-11f1-19f5-1f1e848a5961
begin
	tick # re-run this cell on every clock tick
	mqtt_tags = mqtt_handle === nothing ? Dict{String,String}() : MQTTConnector.latest_values(mqtt_handle)
end

# ╔═╡ 88a1bc1e-ae45-11f1-1b22-97bc031bbf07
if isempty(mqtt_tags)
	md"*No MQTT messages received yet — check the broker connection and `mqtt.topics` in `config/config.toml`.*"
else
	DataFrame(topic = collect(keys(mqtt_tags)), value = collect(values(mqtt_tags)))
end

# ╔═╡ 88a1bc8c-ae45-11f1-21c1-e7bcd414fb29
md"""
### Publish a value

Handy for sending setpoints or commands back out over MQTT.
"""

# ╔═╡ 88a1bcfa-ae45-11f1-06dd-27b9daecc554
md"""
Topic: $(@bind publish_topic PlutoUI.TextField(default="plant/line1/setpoint"))

Value: $(@bind publish_value_str PlutoUI.TextField(default="20.0"))

$(@bind do_publish PlutoUI.Button("Publish to MQTT"))
"""

# ╔═╡ 88a1bd5e-ae45-11f1-2517-3dff0fa89ac2
begin
	do_publish
	if mqtt_handle === nothing
		md"*MQTT not connected — nothing published.*"
	else
		MQTTConnector.publish_value(mqtt_handle, publish_topic, publish_value_str)
		md"Published **$(publish_value_str)** to `$(publish_topic)` at $(Dates.now())."
	end
end

# ╔═╡ 88a1bdcc-ae45-11f1-0f47-471280c25242
md"""
## OPC UA — live nodes

Polls the node IDs listed in `config/config.toml` (`opcua.node_ids`).
Connecting uses Python's `asyncua` library under the hood (via
`PythonCall`/`CondaPkg`) — the first run may take a little longer while the
Python environment is provisioned.
"""

# ╔═╡ 88a1be3a-ae45-11f1-04e9-8527c31ae358
opcua_handle = try
	OPCUAConnector.connect_opcua(cfg.opcua)
catch e
	@warn "Could not connect to the OPC UA server" exception = (e, catch_backtrace())
	nothing
end

# ╔═╡ 88a1be9e-ae45-11f1-2d4d-c9f42be349ac
begin
	tick # re-run this cell on every clock tick
	opcua_values = if opcua_handle === nothing
		Dict{String,Any}()
	else
		try
			OPCUAConnector.read_values(opcua_handle)
		catch e
			@warn "OPC UA read failed" exception = (e, catch_backtrace())
			Dict{String,Any}()
		end
	end
end

# ╔═╡ 88a1bf18-ae45-11f1-2878-516c26cf2fe3
if isempty(opcua_values)
	md"*No OPC UA values available — check the server connection and `opcua.node_ids` in `config/config.toml`.*"
else
	DataFrame(node_id = collect(keys(opcua_values)), value = collect(values(opcua_values)))
end

# ╔═╡ 88a1c0ce-ae45-11f1-0eac-87583938847c
md"""
### Write a node value
"""

# ╔═╡ 88a1c182-ae45-11f1-1b50-5d592325d21f
md"""
Node ID: $(@bind opcua_write_node PlutoUI.TextField(default = isempty(cfg.opcua.node_ids) ? "" : cfg.opcua.node_ids[1]))

Value: $(@bind opcua_write_value PlutoUI.NumberField(-1000:0.1:1000; default=0.0))

$(@bind do_write_opcua PlutoUI.Button("Write to OPC UA"))
"""

# ╔═╡ 88a1c1dc-ae45-11f1-1f31-a7315c2e1b84
begin
	do_write_opcua
	if opcua_handle === nothing
		md"*OPC UA not connected — nothing written.*"
	else
		OPCUAConnector.write_value!(opcua_handle, opcua_write_node, opcua_write_value)
		md"Wrote **$(opcua_write_value)** to `$(opcua_write_node)` at $(Dates.now())."
	end
end

# ╔═╡ 88a1c236-ae45-11f1-1219-7779d3784632
md"""
## InfluxDB — historical trend

Queries InfluxDB (Flux) for a `measurement`/`field` over a time range and
plots it as a simple trend line.
"""

# ╔═╡ 88a1c290-ae45-11f1-139b-df4f7f5158fd
md"""
Measurement: $(@bind trend_measurement PlutoUI.TextField(default="temperature"))

Field: $(@bind trend_field PlutoUI.TextField(default="value"))

Range: $(@bind trend_range PlutoUI.Select(["-15m" => "Last 15 min", "-1h" => "Last hour", "-6h" => "Last 6 hours", "-24h" => "Last 24 hours"]; default="-1h"))
"""

# ╔═╡ 88a1c2fe-ae45-11f1-361c-37aa7cef15f2
begin
	tick # re-run this cell on every clock tick
	trend_df = try
		InfluxConnector.trend(cfg.influx, trend_measurement, trend_field; range = trend_range)
	catch e
		@warn "InfluxDB query failed" exception = (e, catch_backtrace())
		DataFrame()
	end
end

# ╔═╡ 88a1c36c-ae45-11f1-311c-2553028b4fcb
function sparkline_svg(df::DataFrame; width::Int = 640, height::Int = 180, pad::Int = 28, color::String = "#2563eb")
	if !("_value" in names(df)) || nrow(df) < 2
		return md"*No data returned for this range — check the measurement/field names and that InfluxDB has recent data.*"
	end

	ys = Float64.(df[!, "_value"])
	n = length(ys)
	ymin, ymax = extrema(ys)
	yspan = ymax == ymin ? one(ymax) : (ymax - ymin)

	px(i, y) = (
		pad + (i - 1) / (n - 1) * (width - 2pad),
		height - pad - (y - ymin) / yspan * (height - 2pad),
	)

	pts = join((join(px(i, y), ",") for (i, y) in enumerate(ys)), " ")
	latest = round(ys[end]; digits = 3)

	HTML("""
	<svg width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)" style="max-width:100%">
		<polyline fill="none" stroke="$(color)" stroke-width="2" points="$(pts)" />
		<text x="$(width - pad)" y="16" text-anchor="end" font-size="13" fill="currentColor">latest: $(latest)</text>
	</svg>
	""")
end

# ╔═╡ 88a1c3d0-ae45-11f1-2863-6365630f4cb0
sparkline_svg(trend_df)

# ╔═╡ 88a1c42a-ae45-11f1-11a6-81e1d4df52af
trend_df

# ╔═╡ 88a1c486-ae45-11f1-0dfb-219a7dd0f016
md"""
---

Built with Pluto.jl, Mosquitto.jl, HTTP.jl and PythonCall.jl.
See the connector source in `src/connectors/` to extend this dashboard.
"""

# ╔═╡ Cell order:
# ╠═88a0fec8-ae45-11f1-3e9f-c99acd70081c
# ╠═88a1b7c8-ae45-11f1-088b-0fda5f21c95c
# ╠═88a1b8a6-ae45-11f1-2d16-23bb63d6e8db
# ╠═88a1b91c-ae45-11f1-0aea-05e06cd61b51
# ╠═88a1b980-ae45-11f1-3c16-5da6783710f7
# ╠═88a1b9f8-ae45-11f1-1d48-893b12b9890e
# ╠═88a1ba8e-ae45-11f1-114d-9f0ea22be44f
# ╠═88a1bb06-ae45-11f1-37fe-e9d3eec62567
# ╠═88a1bb6a-ae45-11f1-2423-fd1b0a35d242
# ╠═88a1bbc4-ae45-11f1-19f5-1f1e848a5961
# ╠═88a1bc1e-ae45-11f1-1b22-97bc031bbf07
# ╠═88a1bc8c-ae45-11f1-21c1-e7bcd414fb29
# ╠═88a1bcfa-ae45-11f1-06dd-27b9daecc554
# ╠═88a1bd5e-ae45-11f1-2517-3dff0fa89ac2
# ╠═88a1bdcc-ae45-11f1-0f47-471280c25242
# ╠═88a1be3a-ae45-11f1-04e9-8527c31ae358
# ╠═88a1be9e-ae45-11f1-2d4d-c9f42be349ac
# ╠═88a1bf18-ae45-11f1-2878-516c26cf2fe3
# ╠═88a1c0ce-ae45-11f1-0eac-87583938847c
# ╠═88a1c182-ae45-11f1-1b50-5d592325d21f
# ╠═88a1c1dc-ae45-11f1-1f31-a7315c2e1b84
# ╠═88a1c236-ae45-11f1-1219-7779d3784632
# ╠═88a1c290-ae45-11f1-139b-df4f7f5158fd
# ╠═88a1c2fe-ae45-11f1-361c-37aa7cef15f2
# ╠═88a1c36c-ae45-11f1-311c-2553028b4fcb
# ╠═88a1c3d0-ae45-11f1-2863-6365630f4cb0
# ╠═88a1c42a-ae45-11f1-11a6-81e1d4df52af
# ╠═88a1c486-ae45-11f1-0dfb-219a7dd0f016

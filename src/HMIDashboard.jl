module HMIDashboard

include("Config.jl")
using .Config
export Config, AppConfig, InfluxConfig, MQTTConfig, OPCUAConfig, load_config

include(joinpath("connectors", "InfluxConnector.jl"))
using .InfluxConnector
export InfluxConnector

include(joinpath("connectors", "MQTTConnector.jl"))
using .MQTTConnector
export MQTTConnector

include(joinpath("connectors", "OPCUAConnector.jl"))
using .OPCUAConnector
export OPCUAConnector

end # module

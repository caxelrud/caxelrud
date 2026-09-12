using Test
using HMIDashboard
using HMIDashboard.Config
using HMIDashboard.InfluxConnector
using DataFrames
using Dates

@testset "HMIDashboard" begin

    @testset "Config" begin
        cfg = Config.load_config("does/not/exist.toml")
        @test cfg.influx.url == "http://localhost:8086"
        @test cfg.mqtt.port == 1883
        @test cfg.opcua.node_ids == String[]

        mktempdir() do dir
            path = joinpath(dir, "config.toml")
            write(path, """
            [influxdb]
            url = "http://influx.example.com:8086"
            org = "acme"
            bucket = "plant"
            token = "abc123"

            [mqtt]
            host = "broker.example.com"
            port = 8883
            topics = ["a/b", "c/d"]

            [opcua]
            endpoint = "opc.tcp://plc.example.com:4840/"
            node_ids = ["ns=2;i=2"]
            """)
            cfg2 = Config.load_config(path)
            @test cfg2.influx.url == "http://influx.example.com:8086"
            @test cfg2.influx.bucket == "plant"
            @test cfg2.mqtt.port == 8883
            @test cfg2.mqtt.topics == ["a/b", "c/d"]
            @test cfg2.opcua.node_ids == ["ns=2;i=2"]
        end

        withenv("INFLUX_TOKEN" => "env-token", "MQTT_PORT" => "1884") do
            cfg3 = Config.load_config("does/not/exist.toml")
            @test cfg3.influx.token == "env-token"
            @test cfg3.mqtt.port == 1884
        end
    end

    @testset "InfluxConnector line protocol" begin
        line = InfluxConnector.to_line_protocol(
            "temperature",
            ["room" => "lab 1", "sensor" => "A"],
            ["value" => 21.5, "ok" => true, "count" => 3],
        )
        @test startswith(line, "temperature,room=lab\\ 1,sensor=A ")
        @test occursin("value=21.5", line)
        @test occursin("ok=true", line)
        @test occursin("count=3i", line)

        ts = DateTime(2024, 1, 1, 0, 0, 0)
        line_ts = InfluxConnector.to_line_protocol("m", [], ["v" => 1]; timestamp = ts)
        @test occursin(r"v=1i \d+$", line_ts)
    end

    @testset "InfluxConnector annotated CSV parsing" begin
        csv = """
        #group,false,false,true,true,false,false,true,true
        #datatype,string,long,dateTime:RFC3339,dateTime:RFC3339,dateTime:RFC3339,double,string,string
        #default,_result,,,,,,,
        ,result,table,_start,_stop,_time,_value,_field,_measurement
        ,_result,0,2024-01-01T00:00:00Z,2024-01-01T01:00:00Z,2024-01-01T00:00:10Z,21.5,value,temperature
        ,_result,0,2024-01-01T00:00:00Z,2024-01-01T01:00:00Z,2024-01-01T00:00:20Z,21.7,value,temperature
        """
        df = HMIDashboard.InfluxConnector._parse_annotated_csv(csv)
        @test nrow(df) == 2
        @test "_value" in names(df)
        @test df[1, "_value"] == 21.5
    end

end

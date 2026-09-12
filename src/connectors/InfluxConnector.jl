module InfluxConnector

using HTTP
using CSV
using DataFrames
using Dates

import ..Config: InfluxConfig

export query, write_record!, write_records!, trend, latest_value

# ---------------------------------------------------------------------------
# Querying (Flux, InfluxDB v2 /api/v2/query)
# ---------------------------------------------------------------------------

"""
    query(cfg::InfluxConfig, flux::AbstractString) -> DataFrame

Run a Flux query against InfluxDB and return the result as a `DataFrame`.
InfluxDB's "annotated CSV" response can contain several tables (e.g. one per
series); when they share the same columns they are stacked into a single
`DataFrame`, otherwise a `Vector{DataFrame}` is returned (one per table).
"""
function query(cfg::InfluxConfig, flux::AbstractString)
    url = rstrip(cfg.url, '/') * "/api/v2/query?org=" * HTTP.escapeuri(cfg.org)
    body = """{"query": $(repr(flux)), "type": "flux"}"""

    resp = HTTP.post(
        url,
        [
            "Authorization" => "Token $(cfg.token)",
            "Content-Type" => "application/json",
            "Accept" => "application/csv",
        ],
        body;
        status_exception = false,
    )

    if resp.status >= 300
        error("InfluxDB query failed with status $(resp.status): $(String(resp.body))")
    end

    return _parse_annotated_csv(String(resp.body))
end

"""
    _parse_annotated_csv(text) -> DataFrame or Vector{DataFrame}

Parse InfluxDB's "annotated CSV" format: one or more CSV tables separated by
blank lines, each preceded by `#`-prefixed annotation rows (datatype, group,
default) that we simply strip.
"""
function _parse_annotated_csv(text::AbstractString)
    blocks = split(text, r"\r?\n\r?\n")
    tables = DataFrame[]

    for block in blocks
        lines = split(block, r"\r?\n")
        data_lines = [l for l in lines if !isempty(l) && !startswith(l, "#")]
        isempty(data_lines) && continue

        csv_text = join(data_lines, "\n")
        df = CSV.read(IOBuffer(csv_text), DataFrame; missingstring = "")
        # InfluxDB annotated CSV has a leading empty-named column (row index marker).
        if ncol(df) > 0 && (names(df)[1] == "Column1" || names(df)[1] == "")
            select!(df, Not(1))
        end
        push!(tables, df)
    end

    isempty(tables) && return DataFrame()
    length(tables) == 1 && return only(tables)

    # Stack tables that share a schema, otherwise return them separately.
    try
        return reduce(vcat, tables; cols = :union)
    catch
        return tables
    end
end

# ---------------------------------------------------------------------------
# Writing (Line Protocol, InfluxDB v2 /api/v2/write)
# ---------------------------------------------------------------------------

_escape_measurement(s) = replace(string(s), "," => "\\,", " " => "\\ ")
_escape_tag(s) = replace(string(s), "," => "\\,", "=" => "\\=", " " => "\\ ")

function _escape_field_value(v)
    if v isa AbstractString
        return "\"" * replace(v, "\"" => "\\\"") * "\""
    elseif v isa Bool
        return v ? "true" : "false"
    elseif v isa Integer
        return string(v) * "i"
    else
        return string(v)
    end
end

_kvpairs(x) = x
_kvpairs(x::NamedTuple) = pairs(x)

"""
    to_line_protocol(measurement, tags, fields; timestamp=nothing) -> String

Build one InfluxDB line-protocol row. `tags` and `fields` are any
iterable of `key => value` pairs (e.g. `Dict` or `NamedTuple`). `timestamp`,
if given, should be a `DateTime`/`ZonedDateTime` and is sent as nanoseconds
since the epoch.

Prefer a `NamedTuple` (e.g. `(value = 21.5, ok = true)`) over an `Array` or
`Dict` literal when mixing field types: array/dict literals promote all
values to a common numeric type (e.g. `true` silently becomes `1.0`), while
each `NamedTuple` field keeps its own type.
"""
function to_line_protocol(measurement, tags, fields; timestamp = nothing)
    isempty(fields) && error("at least one field is required")

    tag_str = join((string(_escape_tag(k), "=", _escape_tag(v)) for (k, v) in _kvpairs(tags)), ",")
    field_str = join((string(_escape_tag(k), "=", _escape_field_value(v)) for (k, v) in _kvpairs(fields)), ",")

    line = _escape_measurement(measurement)
    isempty(tag_str) || (line *= "," * tag_str)
    line *= " " * field_str

    if timestamp !== nothing
        ns = round(Int64, datetime2unix(timestamp) * 1e9)
        line *= " " * string(ns)
    end

    return line
end

"""
    write_records!(cfg::InfluxConfig, lines::Vector{<:AbstractString})

Write raw line-protocol rows to InfluxDB.
"""
function write_records!(cfg::InfluxConfig, lines::AbstractVector{<:AbstractString})
    url = rstrip(cfg.url, '/') * "/api/v2/write?org=" * HTTP.escapeuri(cfg.org) *
          "&bucket=" * HTTP.escapeuri(cfg.bucket) * "&precision=ns"

    resp = HTTP.post(
        url,
        ["Authorization" => "Token $(cfg.token)", "Content-Type" => "text/plain; charset=utf-8"],
        join(lines, "\n");
        status_exception = false,
    )

    if resp.status >= 300
        error("InfluxDB write failed with status $(resp.status): $(String(resp.body))")
    end

    return nothing
end

"""
    write_record!(cfg, measurement, tags, fields; timestamp=nothing)

Write a single measurement point to InfluxDB.
"""
function write_record!(cfg::InfluxConfig, measurement, tags, fields; timestamp = nothing)
    write_records!(cfg, [to_line_protocol(measurement, tags, fields; timestamp = timestamp)])
end

# ---------------------------------------------------------------------------
# Dashboard convenience helpers
# ---------------------------------------------------------------------------

"""
    trend(cfg, measurement, field; bucket=cfg.bucket, range="-1h", every=nothing) -> DataFrame

Fetch a time series for `field` of `measurement` over the last `range`
(a Flux relative duration such as `"-1h"`, `"-15m"`), optionally downsampled
with `aggregateWindow` when `every` (e.g. `"10s"`) is given. Handy for
feeding a trend chart on the dashboard.
"""
function trend(cfg::InfluxConfig, measurement, field; bucket = cfg.bucket, range = "-1h", every = nothing)
    agg = every === nothing ? "" : """
      |> aggregateWindow(every: $(every), fn: mean, createEmpty: false)
    """
    flux = """
    from(bucket: "$(bucket)")
      |> range(start: $(range))
      |> filter(fn: (r) => r._measurement == "$(measurement)" and r._field == "$(field)")
      $(agg)
      |> keep(columns: ["_time", "_value"])
      |> sort(columns: ["_time"])
    """
    return query(cfg, flux)
end

"""
    latest_value(cfg, measurement, field; bucket=cfg.bucket, lookback="-5m")

Return the most recent value of `field` in `measurement`, or `missing` if
none was found within `lookback`.
"""
function latest_value(cfg::InfluxConfig, measurement, field; bucket = cfg.bucket, lookback = "-5m")
    flux = """
    from(bucket: "$(bucket)")
      |> range(start: $(lookback))
      |> filter(fn: (r) => r._measurement == "$(measurement)" and r._field == "$(field)")
      |> last()
    """
    df = query(cfg, flux)
    (df isa DataFrame && "_value" in names(df) && nrow(df) > 0) || return missing
    return df[end, "_value"]
end

end # module

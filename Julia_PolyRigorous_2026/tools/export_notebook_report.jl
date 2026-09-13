"""
Render `notebooks/ThermoExplorer.jl` to a static HTML report — every
markdown cell and every plot, with each `@bind` slider/dropdown at its
documented default — without needing Pluto's own JS frontend (which loads
its bundle from a CDN at runtime; not usable in a network-sandboxed
environment). Used to generate `docs/ThermoExplorer_notebook_printout.pdf`.

Usage: `julia --project=<repo root> tools/export_notebook_report.jl`
Output: `docs/_notebook_report/report.html` (+ `imgs/`), self-contained
(no external CSS/JS/fonts), so any browser's own "Print to PDF" — or
headless Chrome/Chromium's `--print-to-pdf` — turns it into a PDF.

How it works: this is not a Pluto client. It parses the notebook's own
`# ╔═╡ <uuid>` cell markers (the same format checked for header/footer
UUID consistency whenever the notebook is edited), evaluates each cell's
source directly in `Main` in the order given by the `# ╔═╡ Cell order:`
footer, and renders whatever each cell evaluates to: a `Markdown.MD` cell
becomes HTML, a `Plots.Plot` cell gets `savefig`'d to a PNG and embedded,
everything else is dropped from the report (raw values are visible in the
surrounding markdown's own interpolated text anyway).

Two real bugs surfaced (and fixed here, not worked around) the first time
this ran:

1. `Plots` and `PolyRigorous` both export a function called `density`
   (Plots: KDE plots; PolyRigorous: Sanchez-Lacombe PVT) — referencing the
   bare name is genuinely ambiguous once both are `using`'d in the same
   plain `Main` scope (Pluto's own per-cell scoping avoids this inside the
   real notebook). Resolved by binding `density` explicitly to
   `PolyRigorous.density` below, matching the notebook's own intent.
2. Evaluating a markdown cell directly through the `md` string macro
   throws a spurious internal ParseError on this Julia version whenever a
   string interpolation contains a dotted Unicode field access (e.g.
   `crit.χ_c`) alongside other interpolations in the same cell —
   confirmed in isolation, not an artifact of the cell-splitting here.
   Sidestepped by parsing the cell body as a *plain* interpolated string
   first (Julia's ordinary tokenizer handles that part correctly) and
   handing the fully-interpolated text to `Markdown.parse` (the function,
   not the macro) separately.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Plots
using Markdown
using PolyRigorous
using PlutoUI

gr()
ENV["GKSwstype"] = "100"  # headless GR, no X server needed

const density = PolyRigorous.density  # see bug (1) above

const NOTEBOOK_PATH = joinpath(@__DIR__, "..", "notebooks", "ThermoExplorer.jl")
const REPORT_DIR = joinpath(@__DIR__, "..", "docs", "_notebook_report")
const IMG_DIR = joinpath(REPORT_DIR, "imgs")
mkpath(IMG_DIR)

# Real @bind: pull the actual default out of the PlutoUI element (Slider,
# Select, Scrubbable, CheckBox all carry a `.default` field) instead of the
# notebook's own mock (which falls back to `missing` outside real Pluto) --
# so the report reflects the notebook's own documented default slider
# positions, not hand-transcribed guesses.
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = try
            el.default
        catch
            missing
        end
        el
    end
end

text = read(NOTEBOOK_PATH, String)

cell_re = r"# ╔═╡ ([0-9a-f\-]{36})\n"
markers = collect(eachmatch(cell_re, text))

footer_start = findfirst("# ╔═╡ Cell order:", text)
footer = text[first(footer_start):end]
order_uuids = [m.match for m in eachmatch(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", footer)]

cell_src = Dict{String,String}()
for (i, m) in enumerate(markers)
    uuid = m.captures[1]
    start_pos = m.offset + ncodeunits(m.match)
    end_pos = i < length(markers) ? markers[i+1].offset - 1 : first(footer_start) - 1
    cell_src[uuid] = strip(text[start_pos:end_pos])
end

println("Parsed $(length(cell_src)) cells; $(length(order_uuids)) in Cell order footer.")
@assert Set(keys(cell_src)) == Set(order_uuids) "cell UUIDs and footer UUIDs must match exactly"

# The very first notebook cell (Pkg.activate + using) is already covered by
# this script's own preamble above -- skip it rather than re-eval it (its
# `@__DIR__` would resolve to this script's directory, not the notebook's,
# if evaluated via Meta.parse/eval rather than a real `include`).
const SETUP_CELL_UUID = "1d69aa37-9d2a-4af7-a1d3-9c791a983c12"
@assert haskey(cell_src, SETUP_CELL_UUID)

struct ReportItem
    kind::Symbol   # :markdown, :image, :error
    content::Any
end

report_items = ReportItem[]
img_counter = 0

for uuid in order_uuids
    uuid == SETUP_CELL_UUID && continue
    src = cell_src[uuid]
    isempty(src) && continue
    stripped = strip(src)

    if startswith(stripped, "@bind")
        # Evaluate it (so the bound variable gets its default value for
        # later cells to use) but don't render the Slider/Select widget
        # itself in a static report.
        try
            Base.eval(Main, Meta.parse(src))
        catch e
            push!(report_items, ReportItem(:error, "bind cell $uuid: $e"))
        end
        continue
    end

    if startswith(stripped, "md\"\"\"")
        # See bug (2) in the module docstring above.
        plain_src = replace(src, r"^md" => ""; count=1)
        result = try
            interpolated = Base.eval(Main, Meta.parse(plain_src))
            Markdown.parse(interpolated)
        catch e
            push!(report_items, ReportItem(:error, "markdown cell $uuid: $e"))
            continue
        end
        push!(report_items, ReportItem(:markdown, result))
        continue
    end

    expr = try
        Meta.parse(src)
    catch e
        push!(report_items, ReportItem(:error, "parse error in cell $uuid: $e"))
        continue
    end

    result = try
        Base.eval(Main, expr)
    catch e
        push!(report_items, ReportItem(:error, "runtime error in cell $uuid ($(first(split(src, '\n')))...): $e"))
        continue
    end

    if result isa Markdown.MD
        push!(report_items, ReportItem(:markdown, result))
    elseif result isa Plots.Plot
        global img_counter += 1
        fname = joinpath(IMG_DIR, "plot_$(lpad(img_counter, 3, '0')).png")
        savefig(result, fname)
        push!(report_items, ReportItem(:image, fname))
    end
end

n_md = count(i -> i.kind == :markdown, report_items)
n_img = count(i -> i.kind == :image, report_items)
n_err = count(i -> i.kind == :error, report_items)
println("Executed. markdown=$n_md images=$n_img errors=$n_err")
for item in report_items
    item.kind == :error && println("  ERROR: ", item.content)
end
n_err > 0 && @warn "Some cells failed -- the report is incomplete. See ERROR lines above."

# ---- Assemble one HTML report -----------------------------------------

io = IOBuffer()
println(io, """
<!doctype html><html><head><meta charset="utf-8">
<title>PolyRigorous — ThermoExplorer notebook printout</title>
<style>
body { font-family: -apple-system, Helvetica, Arial, sans-serif; max-width: 860px; margin: 2em auto; padding: 0 1em; line-height: 1.5; color: #222; }
h1,h2,h3 { color: #111; }
h2 { border-top: 2px solid #ddd; padding-top: 1.2em; margin-top: 2em; }
code, pre { background: #f5f5f5; border-radius: 4px; }
pre { padding: 0.7em; overflow-x: auto; font-size: 0.85em; }
code { padding: 0.1em 0.3em; }
img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; margin: 0.5em 0; }
hr { border: none; border-top: 1px solid #ddd; }
</style></head><body>
<h1>PolyRigorous — ThermoExplorer notebook printout</h1>
<p><em>Generated from notebooks/ThermoExplorer.jl, evaluated with each slider/dropdown at its documented default. Every plot below was produced by actually running the notebook's own code in Julia, not recreated by hand.</em></p>
""")

for item in report_items
    if item.kind == :markdown
        println(io, Markdown.html(item.content))
    elseif item.kind == :image
        println(io, "<img src=\"imgs/$(basename(item.content))\">")
    elseif item.kind == :error
        println(io, "<div style=\"color:#a00;font-family:monospace;\">", replace(string(item.content), "<" => "&lt;"), "</div>")
    end
end

println(io, "</body></html>")

report_path = joinpath(REPORT_DIR, "report.html")
open(report_path, "w") do f
    write(f, String(take!(io)))
end
println("REPORT HTML WRITTEN to ", report_path)
println("To get a PDF: open that file in a browser and Print > Save as PDF, or")
println("  <chromium> --headless --disable-gpu --print-to-pdf=out.pdf --print-to-pdf-no-header \\")
println("    --virtual-time-budget=10000 --run-all-compositor-stages-before-draw file://" * report_path)

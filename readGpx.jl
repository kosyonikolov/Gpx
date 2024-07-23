using EzXML
using Dates
using Interpolations
using Plots
using Distances

include("kernelRegression.jl")

function readGpx(fileName::AbstractString)
    lines = readlines(fileName)

    # Remove namespace stuff - breaks findall
    lines[2] = "<gpx>"
    txt = join(lines, "\n")

    # Silence warnings about missing namespace
    txt = replace(txt, "gpxtpx:" => "")

    doc = parsexml(txt)
    segs = findall("//trkseg", doc)
    if lastindex(segs) < 1
        throw("No track segments in $fileName")
    end
    println("Number of segments = $(lastindex(segs))")

    lat = Vector{Float64}()
    lon = Vector{Float64}()
    elev = Vector{Float64}()
    timeStr = Vector{String}()

    # Read only the first segment
    seg = segs[1]
    for node in eachelement(seg)
        if nodename(node) != "trkpt"
            continue
        end
        push!(lat, parse(Float64, node["lat"]))
        push!(lon, parse(Float64, node["lon"]))
        push!(elev, parse(Float64, nodecontent(findfirst("ele", node))))
        push!(timeStr, nodecontent(findfirst("time", node)))
    end

    return timeStr, lat, lon, elev
end

# Matrix with cols = [relTime(s) lat lon elev]
function readGpx2(fileName::AbstractString)
    t, lat, lon, elev = readGpx(fileName)
    df = DateFormat("yyyy-mm-ddTHH:MM:SSZ")
    t0 = DateTime(t[1], df)

    function seconds(tStr::AbstractString)
        v = DateTime(tStr, df) - t0
        return v.value / 1000
    end

    n = lastindex(t)
    result = zeros(n, 4)
    result[:,1] .= seconds.(t)
    result[:,2] .= lat
    result[:,3] .= lon
    result[:,4] .= elev

    return result
end

function calculateVam(t::AbstractVector{<:Number}, elevation::AbstractVector{<:Number};
                      avgTime::Number = 120, sampleDuration::Number = 10)
    elev = linear_interpolation(t, elevation)
    t0 = t[1] + avgTime / 2
    t1 = t[end] - avgTime / 2
    n = floor(Int64, (t1 - t0) / sampleDuration)

    samplePts = range(t0, t1, n)
    vam = zeros(n)
    for i = 1:n
        e0 = elev(samplePts[i] - avgTime / 2)
        e1 = elev(samplePts[i] + avgTime / 2)
        vam[i] = 3600 * (e1 - e0) / avgTime
    end

    return samplePts, vam
end

function plotVam(t::AbstractVector{<:Number}, vam::AbstractVector{<:Number})
    # Convert time to minutes
    t = t ./ 60

    # Pretty ranges
    tMax = 10 * (ceil(Int64, t[end] / 10))
    vamMax = 100 * (ceil(Int64, maximum(vam) / 100))

    p = plot(t, vam, yticks=0:100:vamMax, gridalpha=0.7, gridwidth=1,
             size=(1600,900), minorgrid=true, minorgridalpha=0.4, minorgridwidth = 1,
             xticks=0:10:tMax, label="", title="VAM", linewidth=2,
             ylims = (0, vamMax + 20))

    return p
end

function plotGpxVam(fileName::AbstractString)
    m = readGpx2(fileName)
    tVam, vam = calculateVam(m[:,1], m[:,4])
    p = plotVam(tVam, vam)
    
    id, _ = splitext(fileName)
    outName = "$id.svg"
    println(outName)
    savefig(p, outName)
end

function calculateSpeedFiltered(gpx::AbstractMatrix{<:Number}, σ::Number = 5)
    n, m = size(gpx)
    @assert(m >= 3)

    R = 6372.8 * 1000 # meters

    # Filter latitude and longitude
    t, lat = kernelRegression(gpx[:,1], gpx[:,2], 1, σ)
    t, lon = kernelRegression(gpx[:,1], gpx[:,3], 1, σ)

    # Calculate speed on filtered data
    n = lastindex(t)
    result = zeros(n - 1, 2)
    for i = 1:(n - 1)
        result[i, 1] = t[i] + 0.5
        # !!! Need to swap lat / lon for correct distance !!!
        p0 = (lon[i, 1], lat[i, 1])
        p1 = (lon[i + 1, 1], lat[i + 1, 1])
        d = haversine(p0, p1, R) # meters
        result[i, 2] = 3.6 * d # km/h
    end

    return result
end

function filterGpx(gpx::AbstractMatrix{<:Number}, σ::Number)
    m = size(gpx)[2]
    t, lat = kernelRegression(gpx[:,1], gpx[:,2], 1, σ)
    n = lastindex(t)
    result = zeros(n, m)
    result[:,1] .= t
    result[:,2] .= lat[:,1]
    for i = 3:m
        _, val = kernelRegression(gpx[:,1], gpx[:,i], 1, σ)
        result[:,i] .= val[:,1]
    end

    return result
end

function plotSpeedAndAltitude(gpx::AbstractMatrix{<:Number}, σ::Number = 12)
    minH = minimum(gpx[:,4])
    n = size(gpx)[1]
    
    spd = calculateSpeedFiltered(gpx, σ)

    # Convert time to minutes
    t = spd[:,1] ./ 60

    # Pretty ranges
    tMax = 10 * (ceil(Int64, t[end] / 10))
    vMax = ceil(Int64, maximum(spd[:,2]))

    p = plot(t, spd[:,2], yticks=0:1:vMax, gridalpha=0.5, gridwidth=0.5,
             size=(1600,900), xticks=0:1:tMax, label = "",
             title = "Speed")
    hline!(p, 0:1:vMax, label = "", color = :black, alpha = 0.5, linewidth=0.5)

    p2 = twinx(p)
    plot!(p2, gpx[:,1] ./ 60, minH .* ones(n), fillrange = gpx[:,4], 
          color = :gray, label = "", alpha=0.3)

    return p
end

function plotSpeedAndAltitude(fileName::AbstractString, σ::Number = 12)
    gpx = readGpx2(fileName)
    p = plotSpeedAndAltitude(gpx, σ)

    id, _ = splitext(fileName)
    outName = "$(id)_speed.svg"
    println(outName)
    savefig(p, outName)
end

function plotGpxVam2(gpx::AbstractMatrix{<:Number}; h::Number = 5, Δ::Number = 60, σ::Number = 15)
    @assert(size(gpx)[2] >= 4)
    t0 = gpx[1,1] + Δ/2
    t1 = gpx[end,1] - Δ/2
    ts = range(t0, t1, ceil(Int64, (t1 - t0) / h))
    n = lastindex(ts)
    vam = zeros(n)
    for i = 1:n
        h0 = kernelRegressionPoint(gpx[:,1], gpx[:,4], ts[i] - Δ/2, σ)[1]
        h1 = kernelRegressionPoint(gpx[:,1], gpx[:,4], ts[i] + Δ/2, σ)[1]
        vam[i] = max(0, 3600 * (h1 - h0) / Δ)
    end

    vamMax = 100 * ceil(Int64, maximum(vam) / 100)
    vamTicks = 0:100:vamMax
    tMax = 10 * ceil(Int64, (ts[end] ./ 60) / 10)
    timeTicks = 0:10:tMax

    p = plot(ts ./ 60, vam, yticks = vamTicks, ylims=(0,vamMax),
             label = "", xticks = timeTicks, linewidth=1.2,
             size=(1600,900))
    hline!(p, vamTicks, color = :black, alpha = 0.7, linewidth=0.8, label="")
    timeTicks = collect(timeTicks)
    vline!(p, timeTicks, color = :black, alpha = 0.7, linewidth=0.8, label="")

    subTimeTicks = 0:2:tMax
    vline!(p, subTimeTicks, color = :black, alpha = 0.5, linewidth=0.5, label="")

    p2 = twinx(p)
    minH = minimum(gpx[:,4])
    maxH = maximum(gpx[:,4])
    plot!(p2, gpx[:,1] ./ 60, minH .* ones(size(gpx)[1]), fillrange = gpx[:,4], 
          color = :gray, label = "", alpha=0.3, ylims=(minH,maxH))
    return p
end

function plotGpxVam2(fileName::AbstractString; h::Number = 5, Δ::Number = 60, σ::Number = 15)
    gpx = readGpx2(fileName)
    p = plotGpxVam2(gpx, h = h, Δ = Δ, σ = σ)

    id, _ = splitext(fileName)
    outName = "$(id)_vam.svg"
    println(outName)
    savefig(p, outName)
end

function plotRunGpx(fileName::AbstractString)
    plotSpeedAndAltitude(fileName)
    plotGpxVam2(fileName)
end
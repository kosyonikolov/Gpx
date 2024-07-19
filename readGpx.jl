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
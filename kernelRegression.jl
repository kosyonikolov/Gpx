
function fitPolyWeighted(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, w::AbstractVector{<:Number}, deg::Integer)
    n = lastindex(x)
    @assert(n == lastindex(y) && n == lastindex(w))
    @assert(n >= deg + 1)

    m = zeros(deg + 1, deg + 1)
    b = zeros(deg + 1)
    for i = 1:n
        φ = ones(deg + 1)
        for j = 2:(deg + 1)
            φ[j] = φ[j - 1] * x[i]
        end

        localM = φ * φ'
        localB = φ .* y[i]

        m += w[i] * localM
        b += w[i] * localB
    end

    q = m \ b
    return q
end

function kernelRegression(x::AbstractVector{<:Number}, y::AbstractVector{<:Number},
                          h::Number, sigma::Number; deg::Integer = 1, cutoffSigma::Number = 3)
    @assert(issorted(x))
    n = lastindex(x)
    @assert(n == lastindex(y))

    m = ceil(Int64, (x[end] - x[1]) / h)
    q = range(x[1], x[end], m)
    cutoff = sigma * cutoffSigma

    iLow = 1 # Smallest index such that q - x[i] < cutoff
    result = zeros(m, deg + 1)
    for i = 1:m
        xs = Vector{Float32}()
        ys = Vector{Float32}()
        ws = Vector{Float32}()

        # Find first point
        while iLow <= n && q[i] - x[iLow] > cutoff
            iLow += 1
        end

        j = iLow
        while j <= n && x[j] - q[i] < cutoff
            d = x[j] - q[i]
            push!(xs, d)
            push!(ys, y[j])
            push!(ws, exp(-d^2 / sigma^2))
            j += 1
        end

        poly = fitPolyWeighted(xs, ys, ws, deg)
        #display(poly)
        result[i, :] = poly
    end

    return q, result
end
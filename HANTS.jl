#!/usr/bin/env julia
# -*- coding: utf-8 -*-
"""
module HANTS

Converted from MATLAB

https://mabouali.wordpress.com/projects/harmonic-analysis-of-time-series-hants/

Wout Verhoef
NLR, Remote Sensing Dept.
June 1998

Converted to MATLAB:
Mohammad Abouali (2011) (BSD 2-Clause)

## Modified:

Apply suppression of high amplitudes for near-singular case by
adding a number δ to the diagonal elements of matrix A,
except element (1,1), because the average should not be affected.

Output of reconstructed time series in array yr.

Change call and input arguments to accommodate a base period length (nb).
All frequencies from 1 (base period) until nf are included.

Converted from MATLAB to Julia:
Shen Ruoque (2019) (MIT License)

"""
module HANTS
using LinearAlgebra, Base.Cartesian

export hants, reconstruct

"""
    hants(y, fet, dod, δ; nb, nf, validrange, ts, outlier)

Apply the HANTS process on the series `y`.

## Arguments:

- `y`   : array of input sample values (e.g. NDVI values)
- `fet` : fit error tolerance (points deviating more than fet from curve fit are rejected)
- `dod` : degree of overdeterminedness (iteration stops if number of points reaches the
          minimum required for curve fitting, plus dod). This is a safety measure
- `δ`   : small positive number (e.g. 0.1) to suppress high amplitudes
- `nb`  : length of the base period, measured in virtual samples
          (days, dekads, months, etc.)
- `nf`  : number of frequencies to be considered above the zero frequency
- `validrange` : tuple of valid range minimum and maximum
                 (values outside the valid range are rejeced right away)
- `ts`  : array same size as `y` of time sample indicators (indicates virtual sample
          number relative to the base period); numbers in array `ts` maybe greater
          than nb. Set to 1:length(ts) by default
- `outlier` : 2-character string indicating rejection of high or low outliers.
              Set to `nothing` by default

## Outputs:

- `amp` : returned array of amplitudes, first element is the average of the curve
- `φ`   : returned array of phases, first element is zero
- `yr`  : array holding reconstructed time series
"""
function hants(
    y::AbstractVector{T}, fet, dod, δ;
    nb=length(y), nf=3, validrange=extrema(y[.!isnan.(y)]),
    ts=1:length(y), outlier=nothing
) where {T<:AbstractFloat}

    ni = length(y)
    nr = min(2nf+1, ni)
    mat = Array{T, 2}(undef, nr, ni)
    amp = Vector{T}(undef, nf+1)
    φ = Vector{T}(undef, nf+1)
    local zr, yr

    soutlier = outlier in ["Hi", "High"] ? -1 : outlier in ["Lo", "Low"] ? 1 : 0

    low, high = validrange

    noutmax = ni - nr - dod
    mat[1, :] .= 1

    ang = 2 * (0:nb-1) / nb
    cs = cospi.(ang); sn = sinpi.(ang)
    for i = 1:nf
        for j = 1:ni
            index = 1 + mod(i * (ts[j] - 1), nb)
            mat[2i  , j] = cs[index]
            mat[2i+1, j] = sn[index]
        end
    end

    y_in = replace(y, NaN=>low-1)
    p = ones(ni)
    p[(y_in .< low) .| (y_in .> high)] .= 0
    nout = sum(p .== 0)

    if nout > noutmax return end

    ready = false; nloop = 0; nloopmax = ni

    while (!ready) && (nloop < nloopmax)
        nloop += 1
        za = mat * (p .* y_in)

        A = mat * diagm(0=>p) * mat'
        A += diagm(0=>ones(nr)) * δ
        A[1, 1] -= δ
        zr = A \ za

        yr = mat' * zr
        diffvec = soutlier * (yr - y_in)
        err = p .* diffvec

        rankvec = sortperm(err)

        maxerr = diffvec[Int(rankvec[ni])]
        ready = maxerr ≤ fet || nout == noutmax

        if !ready
            i = ni; j = rankvec[i]
            while p[j] * diffvec[j] > 0.5maxerr && nout < noutmax
                p[j] = 0; nout += 1; i -= 1; j = rank(i)
            end
        end
    end

    amp[1] = zr[1]
    φ[1] = 0

    push!(zr, 0)

    i = 2:2:nr
    ifr = Int.((i .+ 2) ./ 2)
    ra = zr[i]; rb = zr[i.+1]
    amp[ifr] = .√(ra.*ra .+ rb.*rb)
    phase = atand.(rb, ra)
    [phase[j] += 360 for j in eachindex(phase) if phase[j] < 0]
    φ[ifr] = phase

    amp, φ, yr
end

hants(
    y::AbstractVector{<:Integer}, fet, dod, δ; nb, nf, validrange, ts, outlier
) = hants(
    convert(Vector{Float64}, y), fet, dod, δ;
    nb=nb, nf=nf, validrange=validrange,
    ts=ts, outlier=outlier
)


"""
    hants(A, fet, dod, δ; dims::Integer, nb, nf, validrange, ts, outlier)

Apply the HANTS process on the multidimensional array `A`
along the given dimension `dims`.
"""
@generated function hants(
    A::AbstractArray{T,N}, fet, dod, δ;
    dims::Integer, nb=size(A)[dims], nf=3, validrange=extrema(A[.!isnan.(A)]),
    ts=1:size(A)[dims], outlier=nothing
) where{T, N}
    quote
        if dims > N error("dims must less than N") end

        Asize = size(A)
        ampsize = copy(collect(Asize))
        ampsize[dims] = nf + 1
        Ar = Array{Float64, N}(undef, Asize...)
        amp = Array{Float64, N}(undef, ampsize...)
        φ = Array{Float64, N}(undef, ampsize...)

        @nloops $N i A (
            d -> d == dims ? (j_d = (:); if i_d > 1 break end) : j_d = i_d
        ) begin
            y_ = @nref $N A j
            _amp, _φ, _yr = hants(
                y_, fet, dod, δ; nb=nb, nf=nf,
                validrange=validrange, ts=ts, outlier=outlier
            )
            setindex!(amp, _amp, (@ntuple $N j)...)
            setindex!(φ, _φ, (@ntuple $N j)...)
            setindex!(Ar, _yr, (@ntuple $N j)...)
        end

        amp, φ, Ar
    end
end


"""
    reconstruct(amp, φ, nb)

Comput reconstructed series.

## Arguments:

- `amp` : array of amplitudes, first element is the average of the curve
- `φ`   : array of phases, first element is zero
- `nb`  : length of the base period, measured in virtual samples
          (days, dekads, months, etc.)

## Outputs:

- `y`   : reconstructed array of sample values
"""
function reconstruct(
    amp::AbstractVector{T}, φ::AbstractVector{T}, nb::Integer
) where {T<:AbstractFloat}
    nf = length(amp)
    y = zeros(T, nb)
    a_coef = @. amp * cosd(φ)
    b_coef = @. amp * sind(φ)
    for i = 1:nf
        tt = @. (i - 1) * 2 * (0:nb-1) / nb
        y .+= a_coef[i] .* cospi.(tt) .+ b_coef[i] .* sinpi.(tt)
    end
    y
end


"""
    reconstruct(amp, φ, nb; dims::Integer)

Comput reconstructed multidimensional array.
"""
@generated function reconstruct(
    amp::AbstractArray{T,N}, φ::AbstractArray{T,N}, nb::Integer; dims::Integer
) where{T<:AbstractFloat, N}
    quote
        if dims > N error("dims must less than N") end

        ysize = collect(size(amp))
        ysize[dims] = nb
        y = Array{T}(undef, ysize...)

        @nloops $N i y (
            d -> d == dims ? (j_d = (:); if i_d > 1 break end) : j_d = i_d
        ) begin
            _amp = @nref $N amp j
            _φ = @nref $N φ j
            _y = reconstruct(_amp, _φ, nb)
            setindex!(y, _y, (@ntuple $N j)...)
        end

        y
    end
end


end # module HANTS

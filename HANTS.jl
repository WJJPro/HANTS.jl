#!/usr/bin/env julia
# -*- coding: utf-8 -*-
"""
    module HANTS

A Julia module of the Harmonic ANalysis of Time Series (HANTS)

Wout Verhoef
NLR, Remote Sensing Dept.
June 1998

Convert to MATLAB:
Mohammad Abouali (2011) (BSD 2-Clause)

https://mabouali.wordpress.com/projects/harmonic-analysis-of-time-series-hants/

Convert to Julia:
Shen Ruoque (2019) (MIT License)

https://github.com/shenrq/HANTS.jl

## Modified

Use long name instead of acronym.

Add support to multidimensional arrays.
"""
module HANTS
using LinearAlgebra

export hants, reconstruct

"""
    hants(y, fet, dod, δ; nbase, nfreq, validrange, tseries, outlier)

Apply the HANTS process on the series `y`.

## Arguments:

- `y`   : array of input sample values (e.g. NDVI values)
- `fet` : fit error tolerance (points deviating more than fet from curve fit are rejected)
- `dod` : degree of overdeterminedness (iteration stops if number of points reaches the
          minimum required for curve fitting, plus dod). This is a safety measure
- `δ`   : small positive number (e.g. 0.1) to suppress high amplitudes
- `nbase` : length of the base period, measured in virtual samples
            (days, dekads, months, etc.)
- `nfreq` : number of frequencies to be considered above the zero frequency.
            Set to 3 by default
- `validrange` : tuple of valid range minimum and maximum
                 (values outside the valid range are rejeced right away)
                 Set to the extremas by default.
- `tseries` : array same size as `y` of time sample indicators (indicates virtual sample
              number relative to the base period); numbers in array `tseries` maybe greater
              than nbase. Set to 1:length(tseries) by default
- `outlier` : a symbol `:Hi` or `:Lo` indicating rejection of high or low outliers.
              Set to `nothing` by default

## Outputs:

- `A` : returned array of amplitudes, first element is the average of the curve
- `φ` : returned array of phases, first element is zero
- `yrec` : array holding reconstructed time series
"""
function hants(
    y::AbstractVector{T}, fet, dod, δ;
    nbase=length(y), nfreq=3, validrange=extrema(y[.!isnan.(y)]),
    tseries=1:length(y), outlier=nothing
) where {T<:AbstractFloat}

    ny = length(y)
    nr = min(2nfreq+1, ny)
    matirx = Matrix{T}(undef, nr, ny)
    A = Vector{T}(undef, nfreq+1)
    φ = Vector{T}(undef, nfreq+1)
    local zr, yrec

    soutlier = outlier in [:Hi, :High] ? -1 : outlier in [:Lo, :Low] ? 1 : 0

    low, high = validrange

    noutmax = ny - nr - dod
    matirx[1, :] .= 1

    ang = 2 * (0:nbase-1) / nbase
    cs = cospi.(ang); sn = sinpi.(ang)
    for i = 1:nfreq
        j = 1:ny
        index = @. 1 + mod(i * (tseries[j] - 1), nbase)
        matirx[2i  , j] = cs[index]
        matirx[2i+1, j] = sn[index]
    end

    y_in = replace(y, NaN=>low-1)
    p = ones(ny)
    p[(y_in .< low) .| (y_in .> high)] .= 0
    nout = sum(p .== 0)

    if nout > noutmax return end

    ready = false; nloop = 0; nloopmax = ny

    while (!ready) && (nloop < nloopmax)
        nloop += 1
        za = matirx * (p .* y_in)

        arr = matirx * diagm(0=>p) * matirx'
        arr += diagm(0=>ones(nr)) * δ
        arr[1, 1] -= δ
        zr = arr \ za

        yrec = matirx' * zr
        diffvec = soutlier * (yrec - y_in)
        err = p .* diffvec

        rankvec = sortperm(err)

        maxerr = diffvec[Int(rankvec[ny])]
        ready = maxerr ≤ fet || nout == noutmax

        if !ready
            i = ny; j = rankvec[i]
            while p[j] * diffvec[j] > 0.5maxerr && nout < noutmax
                p[j] = 0; nout += 1; i -= 1; j = rank(i)
            end
        end
    end

    A[1] = zr[1]
    φ[1] = 0

    push!(zr, 0)

    i = 2:2:nr
    ifr = 2:1:nr÷2+1
    ra = zr[i]; rb = zr[i.+1]
    A[ifr] = @. √(ra^2 + rb^2)
    phase = atand.(rb, ra)
    replace!(x -> x < 0 ? x + 360 : x, phase)
    φ[ifr] = phase

    A, φ, yrec
end

hants(
    y::AbstractVector{<:Integer}, fet, dod, δ; nbase, nfreq, validrange, tseries, outlier
) = hants(
    convert(Vector{Float64}, y), fet, dod, δ;
    nbase=nbase, nfreq=nfreq, validrange=validrange,
    tseries=tseries, outlier=outlier
)


"""
    hants(arr, fet, dod, δ; dims::Integer, nbase, nfreq, validrange, tseries, outlier)

Apply the HANTS process on the multidimensional array `arr`
along the given dimension `dims`.
"""
function hants(
    arr::AbstractArray{T,N}, fet, dod, δ;
    dims::Integer, nbase=size(arr)[dims], nfreq=3, validrange=extrema(arr[.!isnan.(arr)]),
    tseries=1:size(arr)[dims], outlier=nothing
) where{T, N}
    if dims > N error("dims must less than N") end

    arrsize = size(arr)
    Asize = copy(collect(arrsize))
    Asize[dims] = nfreq + 1

    ny = arrsize[dims]
    pdims = (dims, setdiff(1:ndims(arr), dims)...)
    arr_rec = Array{Float64, N}(undef, arrsize...)
    A = Array{Float64, N}(undef, Asize...)
    φ = Array{Float64, N}(undef, Asize...)
    arrv = permutedims(arr, pdims) |> vec
    arr_recv = PermutedDimsArray(arr_rec, pdims) |> vec
    Av = PermutedDimsArray(A, pdims) |> vec
    φv = PermutedDimsArray(φ, pdims) |> vec

    for i = 1:prod(arrsize)÷ny
        yid = (i-1)*ny+1:i*ny
        Aid = (i-1)*(nfreq+1)+1:i*(nfreq+1)
        Av[Aid], φv[Aid], arr_recv[yid] = hants(
            arrv[yid], fet, dod, δ; nbase=nbase, nfreq=nfreq,
            validrange=validrange, tseries=tseries, outlier=outlier
        )
    end

    A, φ, arr_rec
end


"""
    reconstruct(A, φ, nbase)

Comput reconstructed series.

## Arguments:

- `A` : array of amplitudes, first element is the average of the curve
- `φ` : array of phases, first element is zero
- `nbase` : length of the base period, measured in virtual samples
            (days, dekads, months, etc.)

## Outputs:

- `y` : reconstructed array of sample values
"""
function reconstruct(
    A::AbstractVector{T}, φ, nbase::Integer
) where {T<:AbstractFloat}
    nfreq = length(A)
    y = zeros(T, nbase)
    a_coef = @. A * cosd(φ)
    b_coef = @. A * sind(φ)
    for i = 1:nfreq
        tt = @. (i - 1) * 2 * (0:nbase-1) / nbase
        @. y += a_coef[i] * cospi(tt) + b_coef[i] * sinpi(tt)
    end
    y
end


"""
    reconstruct(A, φ, nbase; dims::Integer)

Comput reconstructed multidimensional array.
"""
function reconstruct(
    A::AbstractArray{T,N}, φ, nbase::Integer; dims::Integer
) where{T<:AbstractFloat, N}
    if dims > N error("dims must less than N") end

    Asize = collect(size(A))
    arrsize = copy(Asize)
    arrsize[dims] = nbase

    nA = Asize[dims]
    pdims = (dims, setdiff(1:ndims(A), dims)...)
    arr = Array{T}(undef, arrsize...)
    arrv = PermutedDimsArray(arr, pdims) |> vec
    Av = permutedims(A, pdims) |> vec
    φv = permutedims(φ, pdims) |> vec

    for i = 1:prod(Asize)÷nA
        yid = (i-1)*nbase+1:i*nbase
        Aid = (i-1)*nA+1:i*nA
        arrv[yid] = reconstruct(Av[Aid], φv[Aid], nbase)
    end

    arr
end


end # module HANTS

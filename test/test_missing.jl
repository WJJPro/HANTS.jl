#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
A missing data set to NaN.
Missing data can also set to a value below the low value or above the high value.
=#
cd(dirname(@__FILE__))
include("../HANTS.jl")
using .HANTS, Plots

y = [
    0.055952, 0.037081, 0.062657, 0.120110, 0.178219, 0.244443, 0.213190, 0.494648,
    0.328767, 0.561776, 0.380300, 0.641233, NaN     , 0.827757, 0.654052, 0.816695,
    0.611424, 0.661557, 0.398067, 0.339881, 0.113957, 0.086790, 0.058600
]

nb = 365
nf = 3
ts = 1:16:365
validrange = (-1, 1)
fet = 0.02
dod = 1
δ = 0.1

amp, φ, yr = hants(
    y, fet, dod, δ; nb=nb, nf=nf, validrange=validrange, ts=ts, outlier=nothing
)
yre = reconstruct(amp, φ, nb)

amp_Hi, φ_Hi, yr_Hi = hants(
    y, fet, dod, δ; nb=nb, nf=nf, validrange=validrange, ts=ts, outlier="Hi"
)
yre_Hi = reconstruct(amp_Hi, φ_Hi, nb)

amp_Lo, φ_Lo, yr_Lo = hants(
    y, fet, dod, δ; nb=nb, nf=nf, validrange=validrange, ts=ts, outlier="Lo"
)
yre_Lo = reconstruct(amp_Lo, φ_Lo, nb)

fig1 = plot(ts, y, lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(ts, yr, lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(ts, yr_Hi, lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(ts, yr_Lo, lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

fig2 = plot(ts, y, lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(1:365, yre, lc=:blue, label="HANTS - None")
plot!(1:365, yre_Hi, lc=:red, label="HANTS - Hi")
plot!(1:365, yre_Lo, lc=:green, label="HANTS - Lo")
title!("Testing HANTS Reconstruct")

plot(fig1, fig2; layout=(1, 2), size=(800,300))

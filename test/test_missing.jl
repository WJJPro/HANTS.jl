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

nbase = 365
nfreq = 3
tseries = 1:16:365
validrange = (-1, 1)
fet = 0.02
dod = 1
δ = 0.1

A, φ, yrec = hants(
    y, fet, dod, δ; nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=nothing
)
y_reconstruct = reconstruct(A, φ, nbase)

A_Hi, φ_Hi, yrec_Hi = hants(
    y, fet, dod, δ; nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Hi
)
y_reconstruct_Hi = reconstruct(A_Hi, φ_Hi, nbase)

A_Lo, φ_Lo, yrec_Lo = hants(
    y, fet, dod, δ; nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Lo
)
y_reconstruct_Lo = reconstruct(A_Lo, φ_Lo, nbase)

fig1 = plot(tseries, y, lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(tseries, yrec, lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(tseries, yrec_Hi, lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(tseries, yrec_Lo, lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

fig2 = plot(tseries, y, lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(1:365, y_reconstruct, lc=:blue, label="HANTS - None")
plot!(1:365, y_reconstruct_Hi, lc=:red, label="HANTS - Hi")
plot!(1:365, y_reconstruct_Lo, lc=:green, label="HANTS - Lo")
title!("Testing HANTS Reconstruct")

plot(fig1, fig2; layout=(1, 2), size=(800,300))

#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Apply HANTS process on an integer array
=#
cd(dirname(@__FILE__))
include("../HANTS.jl")
using .HANTS, Plots

y = [
     559,  370,  626, 1201, 1782, 2444, 2131, 4946,
    3287, 5617, 3803, 6412, 5324, 8277, 6540, 8166,
    6114, 6615, 3980, 3398, 1139,  867,  586
]

nbase = 365
nfreq = 3
tseries = 1:16:365
validrange = (-10000, 10000)
fet = 0.02
dod = 5
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

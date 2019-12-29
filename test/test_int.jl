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

nb = 365
nf = 3
ts = 1:16:365
validrange = (-10000, 10000)
fet = 0.02
dod = 5
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

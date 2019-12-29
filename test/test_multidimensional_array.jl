#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Apply HANTS process on a multidimensional array
=#
cd(dirname(@__FILE__))
include("../HANTS.jl")
using .HANTS, Plots, DelimitedFiles

y = readdlm("2-d_array.csv", ',', Float64, '\n')

nb = 365
nf = 3
ts = 1:16:365
validrange = (-1, 1)
fet = 0.02
dod = 1
δ = 0.1

amp, φ, yr = hants(
    y, fet, dod, δ; dims=2, nb=nb, nf=nf, validrange=validrange, ts=ts, outlier=nothing
)
yre = reconstruct(amp, φ, nb; dims=2)

amp_Hi, φ_Hi, yr_Hi = hants(
    y, fet, dod, δ; dims=2, nb=nb, nf=nf, validrange=validrange, ts=ts, outlier="Hi"
)
yre_Hi = reconstruct(amp_Hi, φ_Hi, nb; dims=2)

amp_Lo, φ_Lo, yr_Lo = hants(
    y, fet, dod, δ; dims=2, nb=nb, nf=nf, validrange=validrange, ts=ts, outlier="Lo"
)
yre_Lo = reconstruct(amp_Lo, φ_Lo, nb; dims=2)

fig1 = plot(ts, y[18, :], lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(ts, yr[18, :], lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(ts, yr_Hi[18, :], lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(ts, yr_Lo[18, :], lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

fig2 = plot(ts, y[18, :], lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(1:365, yre[18, :], lc=:blue, label="HANTS - None")
plot!(1:365, yre_Hi[18, :], lc=:red, label="HANTS - Hi")
plot!(1:365, yre_Lo[18, :], lc=:green, label="HANTS - Lo")
title!("Testing HANTS Reconstruct")

fig = plot(fig1, fig2; layout=(1, 2), size=(800,300))

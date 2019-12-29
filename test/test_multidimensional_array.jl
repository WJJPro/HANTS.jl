#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Apply HANTS process on a multidimensional array
=#
cd(dirname(@__FILE__))
include("../HANTS.jl")
using .HANTS, Plots, DelimitedFiles

y = readdlm("2-d_array.csv", ',', Float64, '\n')

nbase = 365
nfreq = 3
tseries = 1:16:365
validrange = (-1, 1)
fet = 0.02
dod = 1
δ = 0.1

A, φ, yrec = hants(
    y, fet, dod, δ; dims=2, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=nothing
)
y_reconstruct = reconstruct(A, φ, nbase; dims=2)

A_Hi, φ_Hi, yrec_Hi = hants(
    y, fet, dod, δ; dims=2, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Hi
)
y_reconstruct_Hi = reconstruct(A_Hi, φ_Hi, nbase; dims=2)

A_Lo, φ_Lo, yrec_Lo = hants(
    y, fet, dod, δ; dims=2, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Lo
)
y_reconstruct_Lo = reconstruct(A_Lo, φ_Lo, nbase; dims=2)

fig1 = plot(tseries, y[18, :], lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(tseries, yrec[18, :], lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(tseries, yrec_Hi[18, :], lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(tseries, yrec_Lo[18, :], lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

fig2 = plot(tseries, y[18, :], lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(1:365, y_reconstruct[18, :], lc=:blue, label="HANTS - None")
plot!(1:365, y_reconstruct_Hi[18, :], lc=:red, label="HANTS - Hi")
plot!(1:365, y_reconstruct_Lo[18, :], lc=:green, label="HANTS - Lo")
title!("Testing HANTS Reconstruct")

plot(fig1, fig2; layout=(1, 2), size=(800,300))

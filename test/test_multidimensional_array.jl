#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#=
Apply HANTS process on a multidimensional array
=#
cd(dirname(@__FILE__))
include("../HANTS.jl")
using .HANTS, Plots, DelimitedFiles

arr = readdlm("test_data.csv", ',', Float64, '\n')
arr = reshape(arr, 4, 5, 23) # 3d array, the 3rd dimension to apply HANTS

nbase = 365
nfreq = 3
tseries = 1:16:365
validrange = (-1, 1)
fet = 0.02
dod = 1
δ = 0.1

A, φ, arr_rec = hants(
    arr, fet, dod, δ; dims=3, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=nothing
)
arr_reconstruct = reconstruct(A, φ, nbase; dims=3)

A_Hi, φ_Hi, arr_rec_Hi = hants(
    arr, fet, dod, δ; dims=3, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Hi
)
arr_reconstruct_Hi = reconstruct(A_Hi, φ_Hi, nbase; dims=3)

A_Lo, φ_Lo, arr_rec_Lo = hants(
    arr, fet, dod, δ; dims=3, nbase=nbase, nfreq=nfreq, validrange=validrange, tseries=tseries, outlier=:Lo
)
arr_reconstruct_Lo = reconstruct(A_Lo, φ_Lo, nbase; dims=3)

fig1 = plot(tseries, arr[2, 4, :]; lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(tseries, arr_rec[2, 4, :]; lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(tseries, arr_rec_Hi[2, 4, :]; lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(tseries, arr_rec_Lo[2, 4, :]; lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

fig2 = plot(tseries, arr[2, 4, :]; lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(1:365, arr_reconstruct[2, 4, :]; lc=:blue, label="HANTS - None")
plot!(1:365, arr_reconstruct_Hi[2, 4, :]; lc=:red, label="HANTS - Hi")
plot!(1:365, arr_reconstruct_Lo[2, 4, :]; lc=:green, label="HANTS - Lo")
title!("Testing HANTS Reconstruct")

plot(fig1, fig2; layout=(1, 2), size=(800,300))

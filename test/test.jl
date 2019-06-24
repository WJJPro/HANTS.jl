#!/usr/bin/env julia
# -*- coding: utf-8 -*-
include("../HANTS.jl")
using .HANTS, Plots

y = [
    0.055952, 0.037081, 0.062657, 0.120110, 0.178219, 0.244443, 0.213190, 0.494648,
    0.328767, 0.561776, 0.380300, 0.641233, 0.532486, 0.827757, 0.654052, 0.816695,
    0.611424, 0.661557, 0.398067, 0.339881, 0.113957, 0.086790, 0.058600
]

nb = 365
nf = 3
ts = 1:16:365
low = -1
high = 1
fet = 0.02
dod = 1
δ = 0.1

amp, φ, yr = hants(nb, nf, y, ts, nothing, low, high, fet, dod, δ)

amp_Hi, φ_Hi, yr_Hi = hants(nb, nf, y, ts, "Hi", low, high, fet, dod, δ)

amp_lo, φ_lo, yr_Lo = hants(nb, nf, y, ts, "Lo", low, high, fet, dod, δ)

fig = plot(ts, y, lc=:black, label="Original Data", shape=:circle, mc=:black)
plot!(ts, yr, lc=:blue, label="HANTS - None", shape=:circle, mc=:blue)
plot!(ts, yr_Hi, lc=:red, label="HANTS - Hi", shape=:circle, mc=:red)
plot!(ts, yr_Lo, lc=:green, label="HANTS - Lo", shape=:circle, mc=:green)
title!("Testing HANTS Algorithm")

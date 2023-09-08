# DroughtIndices

[![Build Status](https://github.com/semran9/DroughtIndices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/semran9/DroughtIndices.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Current Development Notes
This package in its current state calculates only monthly SPEI based on thornthwaite's equation. More functionality, such as user inputted start times and scales will hopefully be coming soon. Note that in its current, somewhat complete state, this package does not have any safeguards against N/A values, etc. Another note is that in the project's current state, the calculate_SPI() function works (in theory, it has not been tested fully).

## Overview
This package calculates drought indices such as SPEI and SPI. The package implements probability weighted moments and the L-moments procedure to estimate the parameters of the log-logistic distrbution before fitting data to the distribution and standardizing it. In testing, execution of this code is >100x faster than the standard R package.  More on the SPEI index and its calculation can be found at https://spei.csic.es/home.html.

## Main Functions

calculate_SPEI(Precip, Tave, lat, start_month = 1)

Accepts the arguments of an average precipitation vector, an average temperature vector, a latitude and a start month. Returns the the Standarized Precipitation Evapotrasnpiration Index for the provided values (in vector form).

calculate_SPI(Precip)

Accepts average precipitation as an argument and fits it to the log-logistic distribution. Returns the standardized values as a vector.


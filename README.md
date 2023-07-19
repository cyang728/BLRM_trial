# BLRM Dose-Finding Design

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Table of Contents
- [Introduction](#introduction)
- [Functions](#functions)
- [Getting Started](#getting-started)
- [Installation](#installation)

## Introduction

The BLRM Dose-Finding Design is a collection of R functions for implementing the Bayesian Logistic Regression Model (BLRM) in dose-finding clinical trials. The function (BLRM_design) provides the results of the Maximum Tolerated Dose (MTD) and the summary data. 

## Functions

- **interval_prob**: Calculates the probabilities of different intervals based on given probability intervals.
- **checkstop_TI**: Checks stopping conditions for the BLRM dose-finding trial.
- **action_BLRM**: Determines the action to take in the BLRM dose-finding trial.
- **BLRM_design**: Implements the BLRM dose-finding trials.

## Getting Started

To run the BLRM Dose-Finding Design, ensure you have R and JAGS (Just Another Gibbs Sampler) installed on your system. You can install the required packages using the following command:

```R
install.packages("rjags")
```

## Installation
```R
install.packages("devtools")
devtools::install_github("cyang728/BLRM_trial")
```

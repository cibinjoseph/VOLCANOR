#!/bin/bash
# Code to save parameters from current case
# Redirect output to required file

# Parameters from config file and rotor files
cat <(echo '==> config.in <==') config.in <(echo)\
  <(echo '==> rotor??.in <==') rotor??.in

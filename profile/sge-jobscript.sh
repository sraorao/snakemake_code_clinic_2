#!/bin/bash
# properties = {properties}

# exit on first error
set -o errexit
module load Anaconda3/2020.07

{exec_job}

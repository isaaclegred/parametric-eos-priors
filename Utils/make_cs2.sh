#!/bin/bash
cd Utils/Plotting/
mkdir run_cs2
condor_submit cs2_submit
cd ../..

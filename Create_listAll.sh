#!/bin/bash
ls $SPS_LINK/rootDataRun317/streams/$2/$1/$3/CalibratedData/CalibratedData_*  > List/$1
ls $SPS_LINK/rootDataRun317/streams/$2/$1/$3/ProcessedData_*  > List/$1_processed

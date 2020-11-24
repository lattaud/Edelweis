#!/bin/bash
ls /sps/edelweis/rootDataRun317/streams/$2/$1/$3/CalibratedData/CalibratedData_*  > List/$1
ls /sps/edelweis/rootDataRun317/streams/$2/$1/$3/ProcessedData_*  > List/$1_processed

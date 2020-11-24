#!/bin/bash
ls /sps/edelweis/rootDataRun317/streams/$2/$1/$3/CalibratedData/CalibratedData_*_S*[13579]_*_* > List/$1
ls /sps/edelweis/rootDataRun317/streams/$2/$1/$3/ProcessedData_*_S*[13579]_*_*  > List/$1_processed

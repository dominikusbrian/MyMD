#!/bin/bash
TRAJ_STATE="PI"
REC_STATE="PI"

for j in {4100..5000..100}
do
CASE="${j}"
ONE=$((CASE+1))
TWO=$((CASE+2))
THREE=$((CASE+3))
FOUR=$((CASE+4))
cd /xspace/db4271/Triad/FL/RESULT_NVT
sed -i '100001,100002d' triad_thf_${REC_STATE}_TRAJ_${TRAJ_STATE}_NVT_${ONE}.dat
sed -i '100001,100002d' triad_thf_${REC_STATE}_TRAJ_${TRAJ_STATE}_NVT_${TWO}.dat
sed -i '100001,100002d' triad_thf_${REC_STATE}_TRAJ_${TRAJ_STATE}_NVT_${THREE}.dat
sed -i '100001,100002d' triad_thf_${REC_STATE}_TRAJ_${TRAJ_STATE}_NVT_${FOUR}.dat

done


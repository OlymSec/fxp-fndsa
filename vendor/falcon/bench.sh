#!/bin/bash

if [ -f /sys/devices/system/cpu/intel_pstate/no_turbo ]; then
	echo "0" | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo
fi
echo "performance" | sudo tee /sys/devices/system/cpu/cpufreq/policy*/scaling_governor

make clean && \
    CFLAGS="-DFALCON_FPR_MODE=0 -DTEST_REFERENCE" make -j build/speed && \
    mv -f build/speed build/speed_ref_fpnative && \
    echo "################################################################################" > bench.log && \
    echo "#                                 Native (ref)                                 #" >> bench.log && \
    echo "################################################################################" >> bench.log && \
    taskset -c 1 ./build/speed_ref_fpnative 5 >> bench.log
echo "" >>  bench.log

make clean && \
    CFLAGS="-DFALCON_FPR_MODE=1 -DTEST_REFERENCE" make -j build/speed && \
    mv -f build/speed build/speed_ref_fpemu && \
    echo "################################################################################" >> bench.log && \
    echo "#                            Emulated floats (ref)                             #" >> bench.log && \
    echo "################################################################################" >> bench.log && \
    taskset -c 1 ./build/speed_ref_fpemu 5 >> bench.log
echo "" >>  bench.log

make clean && \
    CFLAGS="-DFALCON_FPR_MODE=2" make -j build/speed && \
    mv -f build/speed build/speed_new_fxp && \
    echo "################################################################################" >> bench.log && \
    echo "#                              Fixed point (new)                               #" >> bench.log && \
    echo "################################################################################" >> bench.log && \
    taskset -c 1 ./build/speed_new_fxp 5 >> bench.log
echo "" >>  bench.log

make clean && \
    CFLAGS="-DFALCON_FPR_MODE=0" make -j build/speed && \
    mv -f build/speed build/speed_new_fpnative && \
    echo "################################################################################" >> bench.log && \
    echo "#                                 Native (new)                                 #" >> bench.log && \
    echo "################################################################################" >> bench.log && \
    taskset -c 1 ./build/speed_new_fpnative 5 >> bench.log
echo "" >>  bench.log


if [ -f /sys/devices/system/cpu/intel_pstate/no_turbo ]; then
	echo "1" | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo
fi
echo "powersave" | sudo tee /sys/devices/system/cpu/cpufreq/policy*/scaling_governor

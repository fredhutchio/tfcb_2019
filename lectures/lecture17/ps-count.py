#!/usr/bin/env python

import subprocess

process_count = 0

for process in subprocess.check_output("ps aux", shell=True).split():
    if "sshd" in process:
        process_count += 1

print(process_count)

#!/bin/sh

set -eu

ps aux | grep sshd | wc -l

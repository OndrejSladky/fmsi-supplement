#! /usr/bin/env bash

set -e
set -o pipefail
set -u
#set -f
#PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x


cp input/symlinks/* .

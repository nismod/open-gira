#!/bin/bash

# Given two image files as input, compare them with imagemagick's identfy
# If the 'visual hashes' are equal, exit successfully, otherwise exit 1

set -e

# Imagemagick v6 offers the 'identify' command
# Imagemagick v7+ changes to the 'magick identify' command

if [ -n "$(which magick)" ]; then
    CMD="magick identify"
else
    if [ -n "$(which identify)" ]; then
        CMD="identify"
    else
        echo "Cannot find 'magick' or 'identify' executable, please install imagemagick"
        exit 1
    fi
fi

if [ -z "$1" ]; then
    echo "Require path to image file as argument 1"
    exit 1
fi

if [ -z "$2" ]; then
    echo "Require path to image file as argument 2"
    exit 1
fi

# See https://imagemagick.org/script/identify.php for more information
if [[ $($CMD -quiet -format "%#" "$1") != $($CMD -quiet -format "%#" "$2") ]]; then
    echo "Files $1 and $2 are not a match"
    echo "This could be because they are different formats, or a (subtle) visual difference"
    echo "It should not be because of (embedded) file creation times"
    exit 1
fi

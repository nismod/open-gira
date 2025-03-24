#! /bin/bash

# given two image files as input, compare them with imagemagick's identfy
# if the 'visual hashes' are equal, exit successfully, otherwise exit 1

set -e

if [ -z $(which identify) ]; then
    echo "require executable 'identify', please install imagemagick package"
    exit 1
fi

if [ -z $1 ]; then
    echo "require path to image file as argument 1"
    exit 1
fi

if [ -z $2 ]; then
    echo "require path to image file as argument 2"
    exit 1
fi

# see https://imagemagick.org/script/identify.php for more information
if [[ $(identify -quiet -format "%#" $1) != $(identify -quiet -format "%#" $2) ]]; then
    echo "files $1 and $2 are not a match"
    echo "this could be because they are different formats, or a (subtle) visual difference"
    echo "it should not be because of (embedded) file creation times"
    exit 1
fi

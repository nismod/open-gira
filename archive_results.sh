#/bin/bash

# archive contents of open-gira results directory
# sample usage:
# ./archive_results.sh open-gira/results $DATA/archive

set -e

##### ARGUMENT HANDLING #####

# location of results to be archived, typically open-gira/results
RESULTS_DIR=$1
# location to store the archive (do not include a timestamp)
ARCHIVE_ROOT=$2

usage () {
    echo "sample usage:"
    echo "$0 <results_dir_path> <archive_dir_path>"
}

if [ -z $RESULTS_DIR ]; then
    usage
    exit 1
else
    RESULTS_DIR=$(readlink -f $RESULTS_DIR)
fi
# readlink to get full path however the path was supplied, e.g. . or ../results
# dirname to get the parent directory of that path
OPEN_GIRA_DIR=$(dirname $(readlink -f $RESULTS_DIR))

if [ -z $ARCHIVE_ROOT ]; then
    usage
    exit 1
else
    ARCHIVE_ROOT=$(readlink -f $ARCHIVE_ROOT)
    if [ ! -d $ARCHIVE_ROOT ]; then
        echo "suggested location for archive: $ARCHIVE_ROOT, does not exist, quitting!"
        exit 1
    fi
fi

##### TAR BALLING & COMPRESSING #####

# full date, time, time zone
TIME=$(date +%FT%H%M%S%z)
ARCHIVE_DIR=$ARCHIVE_ROOT/$TIME
if [ -d $ARCHIVE_DIR ]; then
    echo "archive already exists, skipping!"
    exit 1
fi
mkdir $ARCHIVE_DIR

ARCHIVE_FILE="$ARCHIVE_DIR/archive.tar.gz"

# do not back up the input directory (files download from web)
ITEMS_TO_ARCHIVE=$(ls -1 $RESULTS_DIR | grep -vx ^input$)

# if pigz is available, use that for multicore compression
if [ -z $(which pigz) ]; then
    COMPRESSOR=gzip
else
    COMPRESSOR=pigz
fi

# move to results directory (ITEMS_TO_ARCHIVE is relative to this dir)
cd $RESULTS_DIR
echo "archiving results..."
tar cfv - $ITEMS_TO_ARCHIVE | $COMPRESSOR -v > $ARCHIVE_FILE

##### REPORTING #####

cd $OPEN_GIRA_DIR
ARCHIVE_README="$ARCHIVE_DIR/README.md"

# use 'here document' for multiline string while respecting whitespace
cat > $ARCHIVE_README << EOF
open-gira repository state at time of archiving: $TIME
N.B. Not necessarily the same state as when the results were created!

Commit: $(git rev-parse HEAD)

Branch: $(git branch --show-current)

Status:
$(git status --porcelain)

Diff:
$(git diff)
EOF

echo "open-gira results directory (excl. input) archived to $ARCHIVE_DIR"

##### DONE #####

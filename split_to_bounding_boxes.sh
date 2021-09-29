# Takes a OpenStreetmap dataset and slices it into smaller rectangle
# areas.
# Usage: bash split_to_bounding_boxes [data] [ratio]
# Example: bash split_to_bounding_boxes france-latest.osm.pbf 3
# Would generate 9 dataset slices of size L/3 * K/3 where L * K is the
# area of the original bounding box of france-latest.osm.pbf

datafile=$1
ratio=$2

bbox=$(osmium fileinfo $datafile | grep ' \+(-\?[0-9]\+.[0-9]\+,')
xmin0=$(echo $bbox | sed -e 's/[()]//g' | cut -d, -f1)
ymin0=$(echo $bbox | sed -e 's/[()]//g' | cut -d, -f2)
xmax0=$(echo $bbox | sed -e 's/[()]//g' | cut -d, -f3)
ymax0=$(echo $bbox | sed -e 's/[()]//g' | cut -d, -f4)

# Compute longitudinal and latitudinal size of slices
dlong=$(echo "($xmax0 - $xmin0) / $ratio" | bc -l)
dlat=$(echo "($ymax0 - $ymin0) / $ratio" | bc -l)

for i in $(seq 1 $ratio); do
    xmin=$(echo "$xmin0 + ($i - 1) * $dlong" | bc -l)
    xmax=$(echo "$xmin0 + $i * $dlong" | bc -l)
    for j in $(seq 1 $ratio); do
        ymin=$(echo "$ymin0 + ($j - 1)* $dlat" | bc -l)
        ymax=$(echo "$ymin0 + $j * $dlat" | bc -l)

        slug=${datafile%.osm.pbf}
        n=$((j + ratio * $((i - 1))))
        outputfile=${slug}-slice${n}.osm.pbf
        echo "Slicing slice $n"
        osmium extract --no-progress --bbox \
            $xmin,$ymin,$xmax,$ymax $datafile -o $outputfile
    done
done

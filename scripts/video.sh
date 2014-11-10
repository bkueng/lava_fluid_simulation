#! /bin/sh

script_dir="$(dirname $(readlink -f $0))"

if [ "$#" -lt 1 -o "$#" -gt 3 ]
then
    echo "Usage: $0 <rendering_directory> [fps]"
	echo "       rendering_directory eg: output/rendering/sample"
	exit
fi

rendering_dir="$1"
FPS=25
[ "$#" -gt 1 ] && FPS=$2

output_dir="$script_dir/../output/videos"
mkdir "$output_dir" &>/dev/null

suffix=.avi
file_name="$(basename "$rendering_dir")"
# get a non-existing file name
idx=
while [ -f "$output_dir/$file_name$idx$suffix" ]; do
	let "idx=$idx+1"
done
file_name="$file_name$idx"


# -crf 22: high quality constant bitrate
ffmpeg -i "$rendering_dir"/out_%06d.tif -r "$FPS" -crf 22 \
	-vcodec libx264 -threads 6 "$output_dir/$file_name$suffix"

echo "Generated video: $output_dir/$file_name$suffix"


#! /bin/sh

pushd "`dirname $0`" > /dev/null
script_dir="`pwd`"
popd > /dev/null

if [ "$#" -lt 1 -o "$#" -gt 3 ]
then
    echo "Usage: $0 <rendering_directory> [fps [start_frame]]"
	echo "       <rendering_directory> eg: output/rendering/sample"
	echo "       fps: input frame rate (default = 25)"
	echo "       start_frame: first input frame (default = 0)"
	exit
fi

rendering_dir="$1"
FPS=25
[ "$#" -gt 1 ] && FPS=$2
start_frame_params=""
[ "$#" -gt 2 ] && start_frame_params="-start_number $3"

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
ffmpeg -framerate "$FPS" $start_frame_params -i "$rendering_dir"/out_%06d.tif \
	-vf fps=30 -crf 22 \
	-vcodec libx264 -threads 6 "$output_dir/$file_name$suffix"

echo "Generated video: $output_dir/$file_name$suffix"


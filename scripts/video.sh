#! /bin/sh

pushd "`dirname $0`" > /dev/null
script_dir="`pwd`"
popd > /dev/null

usage() {
    echo "Usage: $0 <rendering_directory> [--fps <fps>] [--start <start_frame>]"
	echo "              [--prefix <prefix>]"
	echo "       <rendering_directory> eg: output/rendering/sample"
	echo "       fps: input frame rate (default = 25)"
	echo "       start_frame: first input frame (default = 0)"
	echo "       prefix: input image file name prefix (eg tblur for applied filter)"
	exit
}

if [ "$#" -lt 1 -o "$1" = "--help" ]; then
	usage
fi

rendering_dir="$1"
FPS=25
file_prefix=""
start_frame_params=""
# read parameters
shift
while [ $# -gt 0 ]
do
	key="$1"
	shift
	case $key in
		-h|--help)
			usage
			;;
		--fps)
			FPS="$1"
			shift
			;;
		--start)
			start_frame_params="-start_number $1"
			shift
			;;
		--prefix)
			file_prefix="$1_"
			shift
			;;
		*)
			echo "Error: invalid parameter $key"
			usage
			;;
	esac
done


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
ffmpeg -framerate "$FPS" $start_frame_params -i "$rendering_dir"/"$file_prefix"out_%06d.tif \
	-vf fps=30 -crf 22 \
	-vcodec libx264 -threads 6 "$output_dir/$file_name$suffix"

echo "Generated video: $output_dir/$file_name$suffix"


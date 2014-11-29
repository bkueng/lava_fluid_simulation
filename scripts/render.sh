#! /bin/sh
# this script is used to render the output from the simulation into separate
# image files using a scene template. it uses the same xml config file as the
# simulator

old_pwd="`pwd`"
config_file="$(realpath $1 2>/dev/null)"
pushd "`dirname $0`" > /dev/null
script_name="$0"
#we want to work in the project root directory
cd ..


usage() {
	echo "usage: $script_name <config> [-f <frames>] [-r <renderer>]"
	echo "       <config>             xml config file (under config/)"
	echo "       -f,--frames <frames>"
	echo "                            select which frames to render"
	echo "                            frames is a list of frame indexes, eg:"
	echo "                            0"
	echo "                            \"0 1 2\""
	echo "                            \"\`seq 3 8\`\""
	echo "       -r,--renderer <renderer_binary>"
	echo "                            select renderer (default=$renderer)"
	echo "       -t,--timeout <timeout>"
	echo "                            set rendering timeout in sec (for stuck rendering)"
	echo "       --threads <threads>  maximum number of threads for rendering & image blurring"
	exit -1

}

# parsing xml with xmllint:
# newer versions:
# - attribute:
#   a="$(xmllint --xpath 'string(/config/output/rendering/@constantwidth)' config/sample.xml)"
# - node:
#   a="$(xmllint --xpath '/config/output/text()' config/sample.xml)"
#
# older versions (have no --xpath):
# - attribute:
#   echo "$(echo "cat /config/output/rendering/@constantwidth" | xmllint --nocdata --shell config/sample.xml | sed '1d;$d')"
#   But: this inludes the attribute name!
# - node:
#   a="$(echo "cat /config/output/text()" | xmllint --nocdata --shell config/sample.xml | sed '1d;$d')"

# $1 = xml file, $2 = xpath expression
getxmlattr() {
    echo "cat $2" | xmllint --shell $1 |\
    	sed -n 's/[^\"]*\"\([^\"]*\)\"[^\"]*/\1/gp'
}

# render single frame. $1 = scene file, $2 = data rib file
render() {
	scene_file_tmp="$1"
	file="$2"
	file_name="$(basename -s .rib "$file")"
	frame_nr="${file_name:6}"

	sed -i.tmp "s/frame_[0-9]\{6\}.rib/frame_${frame_nr}.rib/g" "$scene_file_tmp"
	sed -i.tmp "s/out_[0-9]\{6\}.tif/out_${frame_nr}.tif/g" "$scene_file_tmp"
	echo "Rendering $file"
	runs=4 #max retries
	while [ $runs -gt 0 ]; do
		$timeout_prefix $renderer $threads_renderer_params "$scene_file_tmp"
		ret=$?
		[[ $ret -eq 0 ]] && runs=0
		let runs="$runs-1"
	done
}

frames=""
timeout_prefix=""
renderer=rndr #pixie renderer binary
threads_renderer_params=""
[ ! -f "$config_file" ] && echo "Error: no/invalid config file given" && usage
#parse parameters
shift
while [ $# -gt 0 ]
do
	key="$1"
	shift
	case $key in
		-h|--help)
			usage
			;;
		-f|--frames)
			frames="$1"
			shift
			;;
		-r|--renderer)
			renderer="$1"
			shift
			;;
		-t|--timeout)
			timeout_prefix="timeout $1"
			shift
			;;
		--threads)
			threads_renderer_params="-t:$1"
			export OMP_NUM_THREADS=$1
			shift
			;;
		*)
			echo "Error: invalid parameter $key"
			usage
			;;
	esac
done

command -v "$renderer" >/dev/null 2>&1 || \
	{ echo >&2 "Error: renderer \"$renderer\" not found. Is Pixie installed??"; exit 1; }

# data input & rendering directory
base_name="$(basename -s .xml "$config_file")"
data_dir="output/simulation/$base_name"
rendering_dir="output/rendering/$base_name"
mkdir -p "$rendering_dir" &>/dev/null

# additional shaders path
export SHADERS=data/shaders

# parse config file
render_passes=()
for i in `seq 0 15`; do
	scene_file="$(getxmlattr "$config_file" "/config/output/rendering/@pass$i")"
	[ "$scene_file" != "" ] && render_passes+=( "$scene_file" )
done

point_width="$(getxmlattr "$config_file" "/config/output/rendering/@constantwidth")"
format="$(getxmlattr "$config_file" "/config/output/@format")"
height_field_file="$(getxmlattr "$config_file" "/config/simulation/heightfield/@file")"
height_scaling="$(getxmlattr "$config_file" "/config/simulation/heightfield/@scaling")"
[ "$height_scaling" = "" ] && height_scaling=1
height_scaling_larger=`echo "$height_scaling+0.1" | bc`
#height_scaling_larger=`echo "$height_scaling 0.1" | awk '{printf "%f", $1 + $2}'`

echo "Render Passes: ${render_passes[@]}"

# init tmp rendering files: copy the main scene file(s) & apply config variables
scene_files_tmp=()
for render_pass in "${render_passes[@]}"; do
	case "$render_pass" in
		blur*|temporalblur*)
			scene_files_tmp+=( "" )
			;;
		*.rib)
			scene_file_tmp="$(mktemp -t "render_${base_name}.XXXXXXXXXX.rib")"
			scene_files_tmp+=( "$scene_file_tmp" )
			cp "$render_pass" "$scene_file_tmp"

			# iterate 'remove_line_from_scene' & apply to scene file
			i=1
			line="x"
			while [[ "$line" != "" ]]; do
				line="$(getxmlattr "$config_file" "/config/output/rendering/remove_line_from_scene[$i]/@match")"
				if [[ "$line" != "" ]]; then
					line_esc="$(echo "$line" | sed -e 's/[\/&]/\\&/g')"
					sed -i.tmp "/$line_esc/d" "$scene_file_tmp"
				fi
				let i="$i+1"
			done

			# replace directories
			sed -i.tmp "s/OUTPUT_DIRECTORY/${base_name}/g" "$scene_file_tmp"
			# escape the path (-> or better use awk?)
			height_field_file_esc="$(echo "$height_field_file" | sed -e 's/[\/&]/\\&/g')"
			sed -i.tmp "s/HEIGHT_FIELD_FILE/${height_field_file_esc}/g" "$scene_file_tmp"
			if [ "$format" == "point" ]; then
				sed -i.tmp 's/\"constantwidth\"[ ]* \[[ ]* [0-9\.]*/\"constantwidth\" [ '${point_width}'/g' "$scene_file_tmp"
			else
				# delete the line
				sed -i.tmp '/\"constantwidth\"[ ]* \[[ ]* [0-9\.]*/d' "$scene_file_tmp"
			fi
			sed -i.tmp 's/amplitude\"[ ]* \[[ ]* [0-9\.-]*/amplitude\" [ -'${height_scaling}'/g' "$scene_file_tmp"
			sed -i.tmp 's/sphere\"[ ]* \[[ ]* [0-9\.-]*/sphere\" [ '${height_scaling_larger}'/g' "$scene_file_tmp"
			;;
		*)
			echo "Error: invalid render config $render_pass"
			exit
			;;
	esac
done

temporalBlur() {
	args="$1"
	frame_nr="$2"
	file_prefix=""
	filter_size="${args% *}"
	if [[ "$args" == *" "* ]]; then
		file_prefix="${args#* }_"
	fi

	file_out="$rendering_dir/tblur_${file_prefix}out_${frame_nr}.tif"
	file_in_cur="$rendering_dir/${file_prefix}out_${frame_nr}.tif"
	case "$filter_size" in
	1)
		filter_idx=( -1 0 1 )
		filter_factor=( 0.25 0.5 0.25 )
		;;
	2)
		filter_idx=( -2 -1 0 1 2 )
		filter_factor=( 0.15 0.2 0.3 0.2 0.15 )
		;;
	3)
		filter_idx=( -3 -2 -1 0 1 2 3 )
		filter_factor=(0.1052 0.13157 0.15789 0.21052 0.15789 0.13157 0.1052 )
		;;
	*)
		echo "Error: unsupported temporal filter size $filter_size"
		exit -1
		;;
	esac
	convert_args=""

	idx=0
	file_not_found=0
	while [[ "$idx" -lt "${#filter_idx[@]}" ]]; do
		add=${filter_idx[$idx]}
		factor=${filter_factor[$idx]}
		frame_nr_t="$(printf "%06i" $((10#$frame_nr + $add)))"
		file_in="$rendering_dir/${file_prefix}out_${frame_nr_t}.tif"
		if [[ ! -f "$file_in" ]]; then
			echo "Info: $file_in does not exist. no temporal blur is applied"
			file_not_found=1
		fi
		convert_args="$convert_args ( $file_in -evaluate Multiply $factor )"
		let idx="$idx+1"
	done

	if [[ $file_not_found = 1 ]]; then
		# subsequent render passes expect the file to exist
		cp "$file_in_cur" "$file_out"
	else
		#FIXME: alpha channel is handled in a weird way here (by ImageMagick),
		#but it seems to work anyway
		convert $convert_args -evaluate-sequence sum -channel a -negate "$file_out"
	fi
}

renderPass() {
	file="$1" #frame to render
	i="$2" #render pass
	render_pass="${render_passes[$i]}"
	file_name="$(basename -s .rib "$file")"
	frame_nr="${file_name:6}"

	case "$render_pass" in
		temporalblur*)
			echo "Temporal blurring $frame_nr"
			temporalBlur "${render_pass#* }" "$frame_nr"
			;;
		blur*)
			file_base=${render_pass##* }
			file_blur="$rendering_dir/${file_base}_out_${frame_nr}.tif"
			blur="${render_pass% *}"
			echo "Blurring $file_blur"
			convert "$file_blur" -$blur "$file_blur"
			;;
		*.rib)
			scene_file_tmp="${scene_files_tmp[$i]}"
			render "$scene_file_tmp" "$file" "$i"
			;;
	esac
}

# iterate frames & render
i=0
while [[ "$i" -lt "${#render_passes[@]}" ]]; do
	echo "Rendering Pass $i"
	if [ "$frames" = "" ]; then
		for file in $data_dir/frame_*.rib; do
			renderPass "$file" $i
		done
	else
		for frame in $frames; do
			file="$(printf "$data_dir/frame_%06i.rib" $frame)"
			renderPass "$file" $i
		done
	fi
	let i="$i+1"
done

for scene_file_tmp in "${scene_files_tmp[@]}"; do
	[ "$scene_file_tmp" != "" ] && \
		rm "$scene_file_tmp" "$scene_file_tmp".tmp
done


rm -rf images/* &&
make &&
mpiexec -n 3 a.out  &&
./tools/draw_plots_dump_files.py &&
cd images &&
ffmpeg -loglevel fatal -f image2 -pattern_type glob -framerate 12 -i '*.png' -s 400x400 foo.mp4 &&
ffmpeg -loglevel fatal -i foo.mp4 -pix_fmt rgb24 -s qcif -loop 0 -s 400x400 output.gif

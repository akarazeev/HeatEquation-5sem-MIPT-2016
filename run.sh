# rm -rf images/* &&
# make &&
# ./a.out &&
# ./tools/draw_plots_dump_files.py &&
cd images &&
# ffmpeg -loglevel fatal -f image2 -pattern_type glob -framerate 12 -i '*.png' -s 222x620 foo.mp4 &&
ffmpeg -r 12 -y -i "kek_0%2d.png" output.mp4
# ffmpeg -loglevel fatal -i foo.mp4 -pix_fmt rgb24 -s qcif -loop 0 -s 620x222 output.gif

rm M*.png
mv final_simulation.mp4 final_simulation-prev.mp4
gnuplot plotchannel_jb.gplot
avconv -r $1 -i M%05d.png -c:v libx264 final_simulation.mp4
#totem final_simulation.mp4

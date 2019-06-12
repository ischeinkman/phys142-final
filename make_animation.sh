rm out/*
python3 animation.py 
ffmpeg -loglevel warning -r 8 -f image2 -s 1920x1080 -i out/%03d.png -vcodec libx264 -crf 15 -pix_fmt yuv420p -y anim_a.mp4

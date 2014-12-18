BEGIN{save=0;savetime=0}{if (save > 0 && $6 < 0) {print (savetime + ($1 - savetime) * save / (save - $6))}; save = $6; savetime = $1}

D:
cd D:\Program Files\FlightGear

SET FG_ROOT=D:\Program Files\FlightGear\\data
.\\bin\\win32\\fgfs --aircraft=f16 --fdm=network,localhost,5501,5502,5503 --fog-fastest --disable-clouds --start-date-lat=2004:06:01:09:00:00 --disable-sound --in-air --enable-freeze --airport=KSFO --runway=10L --altitude=1000 --heading=113 --offset-distance=4.72 --offset-azimuth=0

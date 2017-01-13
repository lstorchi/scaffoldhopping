
README (Loriano Storchi) loriano@storchi.org
---------------------------------------------------------------------

Simil scaffold hopping test source

Example: 

$ python scaffhop.py "7;11;10;6;5" sc1.sdf

will produce all possibile structures using atoms 7, 11, 10, 6 and 5:

$ babel --gen2D -isdf unique.sdf -osdf unique2d.sdf

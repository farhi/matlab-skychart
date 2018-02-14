# matlab-skychart
A skychart for Matlab: display the sky view

![Image of SkyChart](https://github.com/farhi/matlab-skychart/blob/master/doc/SkyChart.png)

Purpose
-------

**SKYCHART** a class to plot a sky chart with stars/objects
 
This class comptes and plots the sky seen at given location and time. About
43000 stars and 13000 deep sky objects are considered, as well as the Sun, the 
Moon and 7 planets. The actual number of rendered objects depends on the zomm 
level in the sky chart.

You may zoom the plot using the Zoom tool (in the Toolbar). You may as well 
use the drag tool to move the visible area. Right-click shows a contextual
menu with the under-lying object properties (coordinates, type, ...).

Usage
-----

To use this code, type:

```matlab
>> sc = skychart
```

displays the view at the current UTC and location.

**Methods:**

- skychart:   create the view
- date:       set/get the date (UTC)
- load:       load the catalogs. Done at start.
- getplace:   get the current GPS location from the network.
- plot:       plot/replot the view.

You may force a re-computation and replot of the sky view with:

```matlab
>> compute(sc, 'force')
>> plot(sc, 'force)
```

Requirements/Installation
-------------------------
Matlab, no external toolbox.

Just copy the files and go into the src directory. Then type commands above.

Credits
-------

- Local Time to UTC from https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
- Parse JSON from https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
- Amazing work from Eran O. Ofek (MAAT). URL : http://weizmann.ac.il/home/eofek/matlab/
- Stars (~46000) data base from http://astrosci.scimuze.com/stellar_data.htm
- deep sky objects (~13000) from http://klima-luft.de/steinicke/ngcic/ngcic_e.htm


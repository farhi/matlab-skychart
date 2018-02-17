# matlab-skychart
A skychart for Matlab: display the sky view

![Image of SkyChart](https://github.com/farhi/matlab-skychart/blob/master/doc/SkyChart.png)

Purpose
-------

**SKYCHART** a class to plot a sky chart with stars/objects
 
This class computes and plots the sky seen at given location and time. About
43000 stars and 13000 deep sky objects are considered, as well as the Sun, the 
Moon and 7 planets. The actual number of rendered objects depends on the zoom 
level in the sky chart.

You may zoom the plot using the Zoom tool (in the Toolbar). You may as well 
use the drag/pan tool to move the visible area. Right-click shows a contextual
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
- getplace:   get the current GPS location from the network.
- plot:       plot/replot the view.
- connect:    connect to a scope controler
- goto:       send connected scope to selected location
- findobj:    search for a named object and select it

You may force a re-computation and replot of the sky view with:

```matlab
>> compute(sc, 'force')
>> plot(sc, 1)
```

Connecting to a Scope
---------------------

You may connect to a telescope mount using e.g.

```matlab
>> connect(sc, scope)
```

where 'scope' should be an object with methods:

- getstatus: read the mount status, and update the scope properties: scope.ra.h, scope.ra.min, scope.dec.deg, scope.dec.min
- gotoradec(RA,DEC): send the mount to location (RA,DEC)

when scope is omitted, a connection with a Vixen StarBook is attempted. This
controler can be set in 'simulate' mode.

Requirements/Installation
-------------------------

Matlab, no external toolbox.
Just copy the files and go into the directory. Then type commands above.

Credits
-------

- Local Time to UTC from https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
- Parse JSON from https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
- Amazing work from Eran O. Ofek (MAAT). URL : http://weizmann.ac.il/home/eofek/matlab/
- Stars (~46000) data base from http://astrosci.scimuze.com/stellar_data.htm
- deep sky objects (~13000) from http://klima-luft.de/steinicke/ngcic/ngcic_e.htm
- https://fr.mathworks.com/matlabcentral/fileexchange/65944-vixen-starbook-control


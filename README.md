# matlab-skychart
A skychart for Matlab: display the sky view.
Version: 19.08. Download at https://github.com/farhi/matlab-skychart.

![Image of SkyChart](https://github.com/farhi/matlab-skychart/blob/master/%40skychart/doc/SkyChart.png)

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

displays the view at the current UTC and location. The location is obtained from the closes network router (requires internet connection). Else defaults to Grenoble. You can change the location by giving its GPS coordinates with e.g. (in [deg]):

```matlab
>> sc.place = [5.7 45.2]
>> compute(sc,'force'); plot(sc,1);
```

**Methods:**

- **skychart**:   create the view
- **date**:       set/get the date (UTC)
- **getplace**:   get the current GPS location from the network.
- **plot**:       plot/replot the view.
- **connect**:    connect to a scope controler
- **goto**:       send connected scope to selected location
- **findobj**:    search for a named object and select it

You may force a re-computation and replot of the sky view with:

```matlab
>> compute(sc, 'force')
>> plot(sc, 1)
```

or use the "SkyChart/Replot" menu item.

The Chart view
--------------

As seen above, the Chart view represent the visible sky, with stars and deep sky objects. The actual number of visible objects adapts to the zoom level (use the Zoom +/- tool from the tool bar, or the mouse wheel).

The color of stars (filled circles) depend on their spectral type:

- <span style="color:red">red: M K</span>
- <span style="color:green">green: F G</span>
- <span style="color:blue">blue: O B A</span>

Deep sky objects are shown as empty circles:

- <span style="color:red">red: galaxies</span>
- <span style="color:green">green: nebulae</span>
- <span style="color:blue">blue: star clusters</span>

The intensity of the color matches the actual visible magnitude.

The plot can be closed without loosing its information, and re-opened/updated with:

```matlab
>> plot(sc)
```

The plot auto-updates on zoom and pan tools, as well as with Selections (see below).

Connecting to a Scope
---------------------

You may connect to a telescope mount using e.g.

```matlab
>> connect(sc, scope)
```

where 'scope' should be an object with properties/methods:

- scope.getstatus: read the mount status, and update the scope properties: 
- scope.gotoradec(RA,DEC): send the mount to location (RA,DEC)
- properties: scope.ra.h, scope.ra.min, scope.dec.deg, scope.dec.min

We recommend for instance the StarBook and StarGo controllers.

Get the **Vixen StarBook** controller for Matlab at https://github.com/farhi/matlab-starbook. It even works without the physical mount, and then sets itself in 'simulate' mode. A red pointer will then show its current location.

Get the **Avalon StarGo** controller at https://github.com/farhi/matlab-stargo. It also has a simulate mode, and will show as a red pointer as well. 

Selecting objects and Planning observations
-------------------------------------------

Clicking on any visible object on the chart selects it. Right-click shows a contextual menu allowing, when connected, to send the mount immediately to this object (requires a mount controller such as the StarBook - see above). All coordinates (RA, DEC, Alt, Az) are given in degrees.

You can also select objects from the "Planning/Find object..." menu item, or with the 'findobj' method:

```matlab
>> findobj(sc, 'M 42')
```

will search for the given named object, and select it. Valid names include usual names (such as Betelgeuse, Rigel, Capella, Orion nebula, ...), as well as the major Catalogs such as:

- Proper Name
- StarID
- HD (Henry Draper)
- HR (Harvard Revised)
- Gliese (Gliese Catalog of Nearby Stars)
- BayerFlamsteed denomination (Fifth Edition of the Yale Bright Star Catalog)
- M (Messier)
- NGC (New General Catalog)
- IC (Index Catalog)
- CGCG (Zwicky's Catalogue of Galaxies and of Clusters of Galaxies)
- ESO (ESO/Uppsala Survey of the ESO(B)-Atlas)
- IRAS (IRAS catalogue of infrared sources)
- MCG (Morphological Catalogue of Galaxies)
- PGC (Catalog of Principal Galaxies and LEDA)
- UGC (Uppsala General Catalogue of Galaxies)
- ...

You may alternatively add the target to the Planning list. This list can be editted, and you may remove some of its elements, or completely clear it (see the Planning menu). The current selection, as well as the list content are shown on the chart as **white** crosses. 

You can add objects to the list with the "Planning/Add Selected Object", or alternatively:

```matlab
>> findobj(sc, 'M 42');
>> listAdd(sc);
```

When a target is selected, you may create a rectangular grid around it, in order to plan a larger panorama view. The grid size can be customized with the number of steps, and the angular step between lines (for DEC and RA). This can be done, after selecting an object, with:

```matlab
>> findobj(sc, 'M 42');
>> listGrid(sc);
```

which uses a 3x3 grid with angular step 0.75 deg. To change the grid size, use:

```matlab
>> listGrid(sc, sc.selected, n, da);
>> listGrid(sc, sc.findobj('M 42'), n, da);
```

where **n** is the grid binning, e.g. [3 4], and **da** is the angular step between DEC and RA lines, e.g. [0.75 1.22] in [deg]. This step should be small enough to allow partial overlap of the mosaic grid images. When using a focal length F with a camera sensor size S, the field of view (FOV) is da=S/F*57.3 [deg], where S and F should e.g. be in [mm].

With a 1200 mm focal length and an APS-C sensor 23.5x15.6, the FOV is:

- da = 0.74 and 1.12 [deg]
    
With a 400 mm focal length and similar sensor:

- da = 2.23 and 3.36 [deg]

Once you are satisfied with the list, you can set the time of observation between mount moves as the list Period from the "Planning/Set Period" menu entry (time given in [s]), and the method:

```matlab
>>listPeriod(sc, 1800);
```

Then start the Planning with the "Planning/Start" menu entry, which can also be used to Pause/Stop. The same can be done with commands:


```matlab
>>listRun(sc);
```

which can be doubled to pause the execution.

Requirements/Installation
-------------------------

Matlab, no external toolbox.
Just copy the files and go into the directory. Then type commands above.

When behind a proxy/firewall, in order to retrieve your initial GPS location, you may need to define the proxy host and port with e.g.:

```matlab
ProxyHost='proxy.ill.fr'; % Proxy address if you are behind a proxy [e.g. myproxy.mycompany.com or empty]
ProxyPort=8888;           % Proxy port if you are behind a proxy [8888 or 0 or empty]

java.lang.System.setProperty('http.proxyHost', ProxyHost); 
com.mathworks.mlwidgets.html.HTMLPrefs.setUseProxy(true);
com.mathworks.mlwidgets.html.HTMLPrefs.setProxyHost(ProxyHost);
java.lang.System.setProperty('http.proxyPort', num2str(ProxyPort));
com.mathworks.mlwidgets.html.HTMLPrefs.setProxyPort(num2str(ProxyPort));
```
The default GPS location will be used otherwise, and can be changed from the SkyChart menu.

Credits
-------

- Local Time to UTC from https://fr.mathworks.com/matlabcentral/fileexchange/22295-local-time-to-utc
- Parse JSON from https://fr.mathworks.com/matlabcentral/fileexchange/23393--another--json-parser
- Amazing work from Eran O. Ofek (MAAT). URL : http://weizmann.ac.il/home/eofek/matlab/
- Stars (~46000) data base from http://astrosci.scimuze.com/stellar_data.htm
- deep sky objects (~13000) from http://klima-luft.de/steinicke/ngcic/ngcic_e.htm
- Vixen StarBook controller for Matlab at https://github.com/farhi/matlab-starbook.
- Avalon StarGo controller at https://github.com/farhi/matlab-stargo.


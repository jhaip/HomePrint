HomePrint - MIT Mediated Matter
=========

HomePrint is a tool used to break down a .STL model of a building to output code for building-scale 3D printing. Models are sliced up into layers and parameters are adjusted to output the toolpaths used by a KUKA robotic arm for 3D printing.

This project is a derivative/expansion of the [wxblackcat project by Zhigang Liu](http://code.google.com/p/wxblackcat/) which provided the basis for the GUI and slicing of models. I updated the visualization of the models and added additional features needed for 3D printing such as more user friendly parameter adjustment, and predicted print information.

**Requirements:**

- Python 2.7
- wxPython (http://www.wxpython.org)
- pyopengl

**To run:**

1. Navigate to HomePrint folder containing homeprint.py
2. type:

    python homeprint.py

for Mac users:

    python-32 homeprint.py

	

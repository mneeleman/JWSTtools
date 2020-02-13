# create JWST GUI for high-z QSOs
import sys
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import ascii


class ApplicationWindow(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'JWST High-z QSO GUI'
        self.left = 100
        self.top = 100
        self.width = 1200
        self.height = 700
        self.initGraph()
        self.initUI()

    def initGraph(self):
        self.range1 = [[0.5, 5], [0.0, 0.6]]
        self.range2 = [[5, 30], [0, 0.6]]
        self.linecolor = ['black', 'gray']
        self.linescale = [3., 100.]
        self.linestyle = ['-', '-']
        self.linewidth = [2, 2]
        self.linealpha = [1.0, 1.0]
        self.z = [7.54, 7.54]
        self.NWfilters = ['F150W2', 'F322W2', 'F070W', 'F090W', 'F115W',
                          'F150W', 'F200W', 'F277W', 'F356W', 'F444W']
        self.NMfilters = ['F140M', 'F162M', 'F182M', 'F210M', 'F250M',
                          'F300M', 'F335M', 'F360M', 'F410M', 'F430M',
                          'F460M', 'F480M']
        self.NNfilters = ['F164N', 'F187N', 'F212N', 'F323N', 'F405N',
                          'F466N', 'F470N']
        self.NWfcolors = [0.2222, 0.6044, 0.0444, 0.0889, 0.1444,
                          0.2222, 0.3333, 0.5044, 0.6800, 0.8756]
        self.NMfcolors = [0.2000, 0.2489, 0.2933, 0.3556, 0.4444,
                          0.5556, 0.6333, 0.6889, 0.8000, 0.8444,
                          0.9111, 0.9556]
        self.NNfcolors = [0.2533, 0.3044, 0.3600, 0.6067, 0.7889,
                          0.9244, 0.9333]
        self.Mfilters = ['F560W', 'F770W', 'F1000W', 'F1130W', 'F1280W',
                         'F1500W', 'F1800W', 'F2100W', 'F2550W']
        self.Mfcolors = [0.024, 0.108, 0.200, 0.252, 0.312, 0.400,
                         0.520, 0.640, 0.820]
        self.Nsmodes = ['F070LP/G140', 'F100LP/G140', 'F170LP/G235',
                        'F290LP/G395']
        self.Nscolors = [0.108, 0.200, 0.422, 0.800]
        self.Nsdata = [[0.7, 1.27], [0.97, 1.89], [1.66, 3.17], [2.87, 5.27]]
        self.Msmodes = ['SHORT1', 'SHORT2', 'SHORT3', 'SHORT4',
                        'MEDIUM1', 'MEDIUM2', 'MEDIUM3', 'MEDIUM4',
                        'LONG1', 'LONG2', 'LONG3', 'LONG4']
        self.Mscolors = [0.012, 0.124, 0.3, 0.572, 0.044, 0.176, 0.38, 0.7,
                         0.084, 0.228, 0.468, 0.848]
        self.Msdata = [[4.89, 5.75], [7.49, 8.78], [11.53, 13.48],
                       [17.66, 20.92], [5.65, 6.64], [8.65, 10.14],
                       [13.37, 15.63], [20.54, 24.40], [6.52, 7.66],
                       [9.99, 11.71], [15.44, 18.05], [23.95, 28.45]]
        self.filteralpha = 0.3

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # get the individual panels
        self.get_canvas1_layout()
        self.get_canvas2_layout()
        self.get_template_layout()
        self.get_filter_layout()

        # add to layout
        mainLayout = QtWidgets.QGridLayout()
        mainLayout.addWidget(self.canvas1, 0, 0, 2, 3)
        mainLayout.addWidget(self.canvas2, 2, 0, 2, 3)
        mainLayout.addWidget(self.template, 0, 3, 1, 1)
        mainLayout.addWidget(self.filter, 1, 3, 3, 1)
        self.setLayout(mainLayout)

        # draw the intial figure
        self.update_figures()

    def get_canvas1_layout(self):
        self.canvas1 = QtWidgets.QGroupBox('Canvas 1')

        # canvas 1
        self.fig1 = FigureCanvas(Figure(figsize=(1, 1.5)))
        self._ax1 = self.fig1.figure.subplots()

        # range of canvas 1
        self.r1 = []
        for idx, label in enumerate(['X range: ', 'Y range: ']):
            self.r1.append(QtWidgets.QLabel(label, self))
            self.r1.append(QtWidgets.QLineEdit(
                    '{:4.3f}'.format(self.range1[idx][0]), self))
            self.r1.append(QtWidgets.QLabel(' to ', self))
            self.r1.append(QtWidgets.QLineEdit(
                    '{:4.3f}'.format(self.range1[idx][1]), self))
        for idx in np.arange(1, 9, 2):
            self.r1[idx].returnPressed.connect(self.update_range1)

        # go button
        self.updaterange1button = QtWidgets.QPushButton('Update Range', self)
        self.updaterange1button.clicked.connect(self.update_range1)

        # the layout
        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.fig1, 0, 0, 1, 9)
        for idx in np.arange(8):
            layout.addWidget(self.r1[idx], 1, idx)
        layout.addWidget(self.updaterange1button, 1, 8)
        self.canvas1.setLayout(layout)

    def get_canvas2_layout(self):
        self.canvas2 = QtWidgets.QGroupBox('Canvas 2')

        # canvas 2
        self.fig2 = FigureCanvas(Figure(figsize=(1, 1.5)))
        self._ax2 = self.fig2.figure.subplots()

        # range of canvas 2
        self.r2 = []
        for idx, label in enumerate(['X range: ', 'Y range: ']):
            self.r2.append(QtWidgets.QLabel(label, self))
            self.r2.append(QtWidgets.QLineEdit(
                    '{:4.3f}'.format(self.range2[idx][0]), self))
            self.r2.append(QtWidgets.QLabel(' to ', self))
            self.r2.append(QtWidgets.QLineEdit(
                    '{:4.3f}'.format(self.range2[idx][1]), self))
        for idx in np.arange(1, 9, 2):
            self.r2[idx].returnPressed.connect(self.update_range2)

        # go button
        self.updaterange2button = QtWidgets.QPushButton('Update Range', self)
        self.updaterange2button.clicked.connect(self.update_range2)

        # the layout
        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.fig2, 0, 0, 1, 9)
        for idx in np.arange(8):
            layout.addWidget(self.r2[idx], 1, idx)
        layout.addWidget(self.updaterange2button, 1, 8)
        self.canvas2.setLayout(layout)

    def get_template_layout(self):
        self.template = QtWidgets.QTabWidget()

        self.get_QSOtemplate_layout()
        self.get_Galtemplate_layout()

        self.template.addTab(self.QSOtemplate, 'QSO')
        self.template.addTab(self.Galtemplate, 'Galaxy')

    def get_QSOtemplate_layout(self):
        self.QSOtemplate = QtWidgets.QGroupBox('QSO Template')
        self.QSOplot = QtWidgets.QCheckBox('Plot template', self)
        self.QSOz_label = QtWidgets.QLabel('z = ', self)
        self.QSOz = QtWidgets.QLineEdit('{:5.4f}'.format(self.z[0]), self)
        self.QSOcolor_label = QtWidgets.QLabel('Color: ', self)
        self.QSOcolor = QtWidgets.QLineEdit(self.linecolor[0], self)
        self.QSOscale_label = QtWidgets.QLabel('Scale factor: ', self)
        self.QSOscale = QtWidgets.QLineEdit(
                '{:5.2f}'.format(self.linescale[0], self))
        self.QSOlabels = QtWidgets.QCheckBox('Plot labels', self)

        self.QSOplot.setChecked(True)
        self.QSOplot.stateChanged.connect(self.update_figures)
        self.QSOz.returnPressed.connect(self.update_redshifts)
        self.QSOcolor.returnPressed.connect(self.update_linecolors)
        self.QSOscale.returnPressed.connect(self.update_linescales)
        self.QSOlabels.stateChanged.connect(self.update_figures)

        # the layout
        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.QSOplot, 0, 0, 1, 2)
        layout.addWidget(self.QSOlabels, 0, 2, 1, 2)
        layout.addWidget(self.QSOz_label, 1, 0, 1, 1)
        layout.addWidget(self.QSOz, 1, 1, 1, 1)
        layout.addWidget(self.QSOscale_label, 1, 2, 1, 1)
        layout.addWidget(self.QSOscale, 1, 3, 1, 1)
        layout.addWidget(self.QSOcolor_label, 2, 0, 1, 1)
        layout.addWidget(self.QSOcolor, 2, 1, 1, 1)
        self.QSOtemplate.setLayout(layout)

    def get_Galtemplate_layout(self):
        self.Galtemplate = QtWidgets.QGroupBox('Galaxy Template')
        self.Galplot = QtWidgets.QCheckBox('Plot', self)
        self.Galz_label = QtWidgets.QLabel('z = ', self)
        self.Galz = QtWidgets.QLineEdit('{:5.4f}'.format(self.z[1]), self)
        self.Galcolor_label = QtWidgets.QLabel('Color: ', self)
        self.Galcolor = QtWidgets.QLineEdit(self.linecolor[1], self)
        self.Galscale_label = QtWidgets.QLabel('Scale factor: ', self)
        self.Galscale = QtWidgets.QLineEdit(
                '{:5.2f}'.format(self.linescale[1], self))
        self.Gallabels = QtWidgets.QCheckBox('Plot line labels', self)

        self.Galplot.stateChanged.connect(self.update_figures)
        self.Galz.returnPressed.connect(self.update_redshifts)
        self.Galcolor.returnPressed.connect(self.update_linecolors)
        self.Galscale.returnPressed.connect(self.update_linescales)
        self.Gallabels.stateChanged.connect(self.update_figures)

        # the layout
        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.Galplot, 0, 0, 1, 2)
        layout.addWidget(self.Gallabels, 0, 2, 1, 2)
        layout.addWidget(self.Galz_label, 1, 0, 1, 1)
        layout.addWidget(self.Galz, 1, 1, 1, 1)
        layout.addWidget(self.Galscale_label, 1, 2, 1, 1)
        layout.addWidget(self.Galscale, 1, 3, 1, 1)
        layout.addWidget(self.Galcolor_label, 2, 0, 1, 1)
        layout.addWidget(self.Galcolor, 2, 1, 1, 1)
        self.Galtemplate.setLayout(layout)

    def get_filter_layout(self):
        self.filter = QtWidgets.QTabWidget()
        self.filter1 = QtWidgets.QTabWidget()
        self.filter2 = QtWidgets.QTabWidget()
        self.filter3 = QtWidgets.QTabWidget()

        self.get_NWfilters_layout()
        self.get_NMfilters_layout()
        self.get_NNfilters_layout()
        self.get_MIRIfilters_layout()
        self.get_MIRIspec_layout()
        self.get_NIRSpec_layout()

        self.filter.addTab(self.filter1, 'NIRCam')
        self.filter.addTab(self.filter2, 'MIRI')
        self.filter.addTab(self.filter3, 'NIRSpec')

        self.filter1.addTab(self.NIRCamWfilter, 'Wide')
        self.filter1.addTab(self.NIRCamMfilter, 'Medium')
        self.filter1.addTab(self.NIRCamNfilter, 'Narrow')
        self.filter2.addTab(self.MIRIfilter, 'Filters')
        self.filter2.addTab(self.MIRIspec, 'Spectrograph')
        self.filter3.addTab(self.NIRSpec, 'NIRSpec')

    def get_NWfilters_layout(self):
        self.NIRCamWfilter = QtWidgets.QGroupBox('NIRCam wide filters')
        self.NWfilter = []
        for idx, filt in enumerate(self.NWfilters):
            self.NWfilter.append(QtWidgets.QCheckBox(filt, self))
            self.NWfilter.append(QtWidgets.QLabel('Color: ', self))
            self.NWfilter.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.NWfcolors[idx]), self))
        for idx in np.arange(len(self.NWfilters)):
            self.NWfilter[3 * idx].stateChanged.connect(self.update_figures)
            self.NWfilter[3 * idx + 2].returnPressed.connect(
                    self.update_NWfiltercolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.NWfilters)):
            for idx2 in np.arange(3):
                layout.addWidget(self.NWfilter[3 * idx + idx2], idx, idx2)
        self.NIRCamWfilter.setLayout(layout)

    def get_NMfilters_layout(self):
        self.NIRCamMfilter = QtWidgets.QGroupBox('NIRCam medium filters')
        self.NMfilter = []
        for idx, filt in enumerate(self.NMfilters):
            self.NMfilter.append(QtWidgets.QCheckBox(filt, self))
            self.NMfilter.append(QtWidgets.QLabel('Color: ', self))
            self.NMfilter.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.NMfcolors[idx]), self))
        for idx in np.arange(len(self.NMfilters)):
            self.NMfilter[3 * idx].stateChanged.connect(self.update_figures)
            self.NMfilter[3 * idx + 2].returnPressed.connect(
                    self.update_NMfiltercolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.NMfilters)):
            for idx2 in np.arange(3):
                layout.addWidget(self.NMfilter[3 * idx + idx2], idx, idx2)
        self.NIRCamMfilter.setLayout(layout)

    def get_NNfilters_layout(self):
        self.NIRCamNfilter = QtWidgets.QGroupBox('NIRCam narrow filters')
        self.NNfilter = []
        for idx, filt in enumerate(self.NNfilters):
            self.NNfilter.append(QtWidgets.QCheckBox(filt, self))
            self.NNfilter.append(QtWidgets.QLabel('Color: ', self))
            self.NNfilter.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.NNfcolors[idx]), self))
        for idx in np.arange(len(self.NNfilters)):
            self.NNfilter[3 * idx].stateChanged.connect(self.update_figures)
            self.NNfilter[3 * idx + 2].returnPressed.connect(
                    self.update_NNfiltercolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.NNfilters)):
            for idx2 in np.arange(3):
                layout.addWidget(self.NNfilter[3 * idx + idx2], idx, idx2)
        self.NIRCamNfilter.setLayout(layout)

    def get_MIRIfilters_layout(self):
        self.MIRIfilter = QtWidgets.QGroupBox('MIRI filters')
        self.Mfilter = []
        for idx, filt in enumerate(self.Mfilters):
            self.Mfilter.append(QtWidgets.QCheckBox(filt, self))
            self.Mfilter.append(QtWidgets.QLabel('Color: ', self))
            self.Mfilter.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.Mfcolors[idx]), self))
        for idx in np.arange(len(self.Mfilters)):
            self.Mfilter[3 * idx].stateChanged.connect(self.update_figures)
            self.Mfilter[3 * idx + 2].returnPressed.connect(
                    self.update_Mfiltercolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.Mfilters)):
            for idx2 in np.arange(3):
                layout.addWidget(self.Mfilter[3 * idx + idx2], idx, idx2)
        self.MIRIfilter.setLayout(layout)

    def get_MIRIspec_layout(self):
        self.MIRIspec = QtWidgets.QGroupBox('MIRI spectral ranges')
        self.Mspec = []
        for idx, filt in enumerate(self.Msmodes):
            self.Mspec.append(QtWidgets.QCheckBox(filt, self))
            self.Mspec.append(QtWidgets.QLabel('Color: ', self))
            self.Mspec.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.Mscolors[idx]), self))
        for idx in np.arange(len(self.Msmodes)):
            self.Mspec[3 * idx].stateChanged.connect(self.update_figures)
            self.Mspec[3 * idx + 2].returnPressed.connect(
                    self.update_Mspeccolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.Msmodes)):
            for idx2 in np.arange(3):
                layout.addWidget(self.Mspec[3 * idx + idx2], idx, idx2)
        self.MIRIspec.setLayout(layout)

    def get_NIRSpec_layout(self):
        self.NIRSpec = QtWidgets.QGroupBox('NIRSpec spectral ranges')
        self.Nspec = []
        for idx, filt in enumerate(self.Nsmodes):
            self.Nspec.append(QtWidgets.QCheckBox(filt, self))
            self.Nspec.append(QtWidgets.QLabel('Color: ', self))
            self.Nspec.append(QtWidgets.QLineEdit(
                    '{:6.4f}'.format(self.Nscolors[idx]), self))
        for idx in np.arange(len(self.Nsmodes)):
            self.Nspec[3 * idx].stateChanged.connect(self.update_figures)
            self.Nspec[3 * idx + 2].returnPressed.connect(
                    self.update_Nspeccolor)

        # the layout
        layout = QtWidgets.QGridLayout()
        for idx in np.arange(len(self.Nsmodes)):
            for idx2 in np.arange(3):
                layout.addWidget(self.Nspec[3 * idx + idx2], idx, idx2)
        self.NIRSpec.setLayout(layout)

    def update_redshifts(self):
        try:
            newvalue = np.float(self.QSOz.text())
            self.z[0] = newvalue
        except ValueError:
            print(self.QSOz.text() + 'Not a valid value!')
            self.QSOz.setText(str(self.z[0]))
        try:
            newvalue = np.float(self.Galz.text())
            self.z[1] = newvalue
        except ValueError:
            print(self.Galz.text() + 'Not a valid value!')
            self.Galz.setText(str(self.z[1]))
        self.update_figures()

    def update_linecolors(self):
        if self.QSOcolor.text() in mpl.colors.cnames.keys():
            self.linecolor[0] = self.QSOcolor.text()
            self.update_figures()
        else:
            print(self.QSOcolor.text() + ' is not a valid color!')
            self.QSOcolor.setText(self.linecolor[0])

        if self.Galcolor.text() in mpl.colors.cnames.keys():
            self.linecolor[1] = self.Galcolor.text()
            self.update_figures()
        else:
            print(self.Galcolor.text() + ' is not a valid color!')
            self.Galcolor.setText(self.linecolor[1])
        self.update_figures()

    def update_linescales(self):
        try:
            newvalue = np.float(self.QSOscale.text())
            self.linescale[0] = newvalue
        except ValueError:
            print(self.QSOscale.text() + 'Not a valid value!')
            self.QSOscale.setText(str(self.linescale[0]))
        try:
            newvalue = np.float(self.Galscale.text())
            self.linescale[1] = newvalue
        except ValueError:
            print(self.Galscale.text() + 'Not a valid value!')
            self.Galscale.setText(str(self.linescale[1]))
        self.update_figures()

    def update_range1(self):
        try:
            newvalue = np.float(self.r1[1].text())
            self.range1[0][0] = newvalue
        except ValueError:
            print(self.r1[1].text() + 'Not a valid value!')
            self.r1[1].setText(str(self.range1[0][0]))
        try:
            newvalue = np.float(self.r1[3].text())
            self.range1[0][1] = newvalue
        except ValueError:
            print(self.r1[3].text() + 'Not a valid value!')
            self.r1[3].setText(str(self.range1[0][1]))
        try:
            newvalue = np.float(self.r1[5].text())
            self.range1[1][0] = newvalue
        except ValueError:
            print(self.r1[5].text() + 'Not a valid value!')
            self.r1[5].setText(str(self.range1[1][0]))
        try:
            newvalue = np.float(self.r1[7].text())
            self.range1[1][1] = newvalue
        except ValueError:
            print(self.r1[7].text() + 'Not a valid value!')
            self.r1[7].setText(str(self.range1[1][1]))
        self.update_figures()

    def update_range2(self):
        try:
            newvalue = np.float(self.r2[1].text())
            self.range2[0][0] = newvalue
        except ValueError:
            print(self.r2[1].text() + 'Not a valid value!')
            self.r2[1].setText(str(self.range2[0][0]))
        try:
            newvalue = np.float(self.r2[3].text())
            self.range2[0][1] = newvalue
        except ValueError:
            print(self.r2[3].text() + 'Not a valid value!')
            self.r2[3].setText(str(self.range2[0][1]))
        try:
            newvalue = np.float(self.r2[5].text())
            self.range2[1][0] = newvalue
        except ValueError:
            print(self.r2[5].text() + 'Not a valid value!')
            self.r2[5].setText(str(self.range2[1][0]))
        try:
            newvalue = np.float(self.r2[7].text())
            self.range2[1][1] = newvalue
        except ValueError:
            print(self.r2[7].text() + 'Not a valid value!')
            self.r2[7].setText(str(self.range2[1][1]))
        self.update_figures()

    def update_NWfiltercolor(self):
        for idx in np.arange(len(self.NWfilters)):
            try:
                newvalue = np.float(self.NWfilter[3 * idx + 2].text())
                self.NWfcolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.NWfilter[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.NWfcolors[idx] = self.NWfilter[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.NWfilter[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.NWfilter[3 * idx+2].setText('black')
                    self.update_figures()

    def update_NMfiltercolor(self):
        for idx in np.arange(len(self.NMfilters)):
            try:
                newvalue = np.float(self.NMfilter[3 * idx + 2].text())
                self.NMfcolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.NMfilter[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.NMfcolors[idx] = self.NMfilter[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.NMfilter[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.NMfilter[3 * idx+2].setText('black')
                    self.update_figures()

    def update_NNfiltercolor(self):
        for idx in np.arange(len(self.NNfilters)):
            try:
                newvalue = np.float(self.NNfilter[3 * idx + 2].text())
                self.NNfcolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.NNfilter[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.NNfcolors[idx] = self.NNfilter[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.NNfilter[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.NNfilter[3 * idx+2].setText('black')
                    self.update_figures()

    def update_Mfiltercolor(self):
        for idx in np.arange(len(self.Mfilters)):
            try:
                newvalue = np.float(self.Mfilter[3 * idx + 2].text())
                self.Mfcolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.Mfilter[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.Mfcolors[idx] = self.Mfilter[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.Mfilter[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.Mfilter[3 * idx+2].setText('black')
                    self.update_figures()

    def update_Mspeccolor(self):
        for idx in np.arange(len(self.Msmodes)):
            try:
                newvalue = np.float(self.Mspec[3 * idx + 2].text())
                self.Mscolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.Mspec[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.Mscolors[idx] = self.Mspec[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.Mspec[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.Mspec[3 * idx+2].setText('black')
                    self.update_figures()

    def update_Nspeccolor(self):
        for idx in np.arange(len(self.Nsmodes)):
            try:
                newvalue = np.float(self.Nspec[3 * idx + 2].text())
                self.Nscolors[idx] = newvalue
                self.update_figures()
            except ValueError:
                if self.Nspec[3 * idx+2].text() in mpl.colors.cnames.keys():
                    self.Nscolors[idx] = self.Nspec[3 * idx+2].text()
                    self.update_figures()
                else:
                    print(self.Nspec[3 * idx+2].text() +
                          ' is not a valid color!')
                    self.Nspec[3 * idx+2].setText('black')
                    self.update_figures()

    def update_figures(self):

        # clear figures
        self._ax1.clear()
        self._ax2.clear()

        # plot templates
        if self.QSOplot.isChecked():
            self.plot_qso()
        if self.Galplot.isChecked():
            self.plot_galaxy()

        # plot filters
        self.plot_filters()

        # define range, axis-label and draw
        self._ax1.set_xlim(self.range1[0][0], self.range1[0][1])
        self._ax1.set_ylim(self.range1[1][0], self.range1[1][1])
        self._ax1.figure.canvas.draw()

        self._ax2.set_xlim(self.range2[0][0], self.range2[0][1])
        self._ax2.set_ylim(self.range2[1][0], self.range2[1][1])
        self._ax2.figure.canvas.draw()

    def plot_qso(self):
        qso = ascii.read("templates/z75qso.dat")
        wave_mu = qso['lambda_mu'] / (1 + 7.5)
        wave_mu *= (1 + self.z[0])
        flux = qso['flux_mjy']

        self._ax1.plot(wave_mu, flux * self.linescale[0], ls=self.linestyle[0],
                       color=self.linecolor[0], alpha=self.linealpha[0],
                       lw=self.linewidth[0])
        self._ax2.plot(wave_mu, flux * self.linescale[0], ls=self.linestyle[0],
                       color=self.linecolor[0], alpha=self.linealpha[0],
                       lw=self.linewidth[0])

        if self.QSOlabels.isChecked():
            self.add_lines(self.z[0], 0)

    def plot_galaxy(self):
        s99 = ascii.read("templates/starburst99_host_models_z75_lines.txt")
        wave_mu = s99['wave'] / (1 + 7.5)
        wave_mu *= (1 + self.z[1])
        flux = s99["300Myr_0.2"]

        self._ax1.plot(wave_mu, flux * self.linescale[1], ls=self.linestyle[1],
                       color=self.linecolor[1], alpha=self.linealpha[1],
                       lw=self.linewidth[1])
        self._ax2.plot(wave_mu, flux * self.linescale[1], ls=self.linestyle[1],
                       color=self.linecolor[1], alpha=self.linealpha[1],
                       lw=self.linewidth[1])

        if self.Gallabels.isChecked():
            self.add_lines(self.z[1], 1)

    def plot_filters(self):

        # NIRCam-Wide
        for idx, filt in enumerate(self.NWfilters):
            if self.NWfilter[3 * idx].isChecked():
                # read the data
                Folder = ('./filters/jwst/nircam_throughputs/modAB_mean/' +
                          'nrc_plus_ote/')
                Ext = '_NRC_and_OTE_ModAB_mean.txt'
                data = ascii.read(Folder+filt+Ext)
                if type(self.NWfcolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.NWfcolors[idx])
                else:
                    color = self.NWfcolors[idx]
                self._ax1.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)
                self._ax2.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)

        # NIRCam-Medium
        for idx, filt in enumerate(self.NMfilters):
            if self.NMfilter[3 * idx].isChecked():
                # read the data
                Folder = ('./filters/jwst/nircam_throughputs/modAB_mean/' +
                          'nrc_plus_ote/')
                Ext = '_NRC_and_OTE_ModAB_mean.txt'
                data = ascii.read(Folder+filt+Ext)
                if type(self.NMfcolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.NMfcolors[idx])
                else:
                    color = self.NMfcolors[idx]
                self._ax1.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)
                self._ax2.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)

        # NIRCam-Narrow
        for idx, filt in enumerate(self.NNfilters):
            if self.NNfilter[3 * idx].isChecked():
                # read the data
                Folder = ('./filters/jwst/nircam_throughputs/modAB_mean/' +
                          'nrc_plus_ote/')
                Ext = '_NRC_and_OTE_ModAB_mean.txt'
                data = ascii.read(Folder+filt+Ext)
                if type(self.NNfcolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.NNfcolors[idx])
                else:
                    color = self.NNfcolors[idx]
                self._ax1.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)
                self._ax2.fill_between(data['microns'], data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)

        # MIRI
        for idx, filt in enumerate(self.Mfilters):
            if self.Mfilter[3 * idx].isChecked():
                # read the data
                Folder = ('./filters/jwst/miri/JWST-MIRI.')
                Ext = '.dat'
                data = ascii.read(Folder+filt+Ext)
                if type(self.Mfcolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.Mfcolors[idx])
                else:
                    color = self.Mfcolors[idx]
                self._ax1.fill_between(data['angstrom'] / 1e4,
                                       data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)
                self._ax2.fill_between(data['angstrom'] / 1e4,
                                       data['throughput'],
                                       edgecolor=color, facecolor=color,
                                       alpha=self.filteralpha)

        # MIRI - spec
        for idx, filt in enumerate(self.Msmodes):
            if self.Mspec[3 * idx].isChecked():
                if type(self.Mscolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.Mscolors[idx])
                else:
                    color = self.Mscolors[idx]
                self._ax1.axvspan(self.Msdata[idx][0], self.Msdata[idx][1],
                                  ymin=(0.1 + 0.07 * np.mod(idx, 4)),
                                  ymax=(0.15 + 0.07 * np.mod(idx, 4)),
                                  edgecolor=color, facecolor=color,
                                  alpha=self.filteralpha)
                self._ax2.axvspan(self.Msdata[idx][0], self.Msdata[idx][1],
                                  ymin=(0.1 + 0.05 * np.mod(idx, 4)),
                                  ymax=(0.15 + 0.07 * np.mod(idx, 4)),
                                  edgecolor=color, facecolor=color,
                                  alpha=self.filteralpha)

        # NIRSPEC
        for idx, filt in enumerate(self.Nsmodes):
            if self.Nspec[3 * idx].isChecked():
                if type(self.Nscolors[idx]) is np.float:
                    color = plt.cm.Spectral_r(self.Nscolors[idx])
                else:
                    color = self.Nscolors[idx]
                self._ax1.axvspan(self.Nsdata[idx][0], self.Nsdata[idx][1],
                                  ymin=(0.1 + 0.07 * idx),
                                  ymax=(0.15 + 0.07 * idx),
                                  edgecolor=color, facecolor=color,
                                  alpha=self.filteralpha)
                self._ax2.axvspan(self.Nsdata[idx][0], self.Nsdata[idx][1],
                                  ymin=(0.1 + 0.05 * idx),
                                  ymax=(0.15 + 0.07 * idx),
                                  edgecolor=color, facecolor=color,
                                  alpha=self.filteralpha)

    def add_lines(self, z, idx):
        line_label = [r'Ly$\alpha$+N$\,$V', r'C$\,$IV', r'Mg$\,$II',
                      r'C$\,$III] + Fe$\,$III',
                      r'H$\,\beta$' + r'+ [O$\,$III]', r'H$\, \alpha$',
                      r'Pa$\, \alpha$', r'Pa$\, \beta$',
                      r'Pa$\, \gamma+\delta$']
        line_wave = [1270 / 1e4 - 0.01, 1546.15 / 1e4, 2800.26 / 1e4,
                     1905.97 / 1e4, 4835 / 1e4 + 0.01, 6564.93 / 1e4,
                     1.87,  1.28, 1.06]
        for label, wave in zip(line_label, line_wave):
            if (wave * (1 + z) > self.range1[0][0] and
                    wave * (1 + z) < self.range1[0][1]):
                self._ax1.text(wave * (1 + z), 0.9 * self.range1[1][1],
                               label, color=self.linecolor[idx], fontsize=8,
                               ha='center')
            if (wave * (1 + z) > self.range2[0][0] and
                    wave * (1 + z) < self.range2[0][1]):
                self._ax2.text(wave * (1 + z), 0.9 * self.range1[1][1],
                               label, color=self.linecolor[idx], fontsize=8,
                               ha='center')


def main():
    app = QtWidgets.QApplication(sys.argv)
    if len(sys.argv) != 2:
        print(sys.argv)
        print('Please supply the fits file for the data')
        sys.exit()
    else:
        AppWin = ApplicationWindow(sys.argv[1])
        AppWin.show()
        sys.exit(app.exec_())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    AppWin = ApplicationWindow()
    AppWin.show()
    sys.exit(app.exec_())
